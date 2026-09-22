"""Mathematical regression tests; run with PYTHONPATH=src sage -python -m unittest discover -s tests -v."""

import unittest
from unittest.mock import patch
from collections import Counter,defaultdict
from sage.all import AA,QQ,ZZ,RealBallField,PolynomialRing,NumberField,matrix,vector,identity_matrix

from rigorous_cnt import (NumberFieldContext,RayConditions,RandomBits,Limits,
                          ResourceLimit,EmptySupport,UnresolvedBound,WalkParameters,Algorithm1,algorithm2)
from rigorous_cnt.numerics import endpoints,nearest,bernoulli_exp_pi,norm_one_exp_approx
from rigorous_cnt.gaussian import discrete_gaussian
from rigorous_cnt.lattice import gram_schmidt,lll_transform,hkz_transform,shortest_vector_coefficients
from rigorous_cnt.geometry import ProductRegion,UniformLatticeSampler
from rigorous_cnt.primes import PrimeIdealSampler
from rigorous_cnt.parameters import distortion_grid,corollary65_steps


R=PolynomialRing(QQ,'t')
t=R.gen()


class NumericalTests(unittest.TestCase):
    def test_endpoints_are_directed_exact_rationals(self):
        lo,hi=endpoints(RealBallField(128)(2).sqrt())
        self.assertLess(lo*lo,2)
        self.assertGreater(hi*hi,2)
        self.assertTrue(ZZ(lo.denominator()).is_power_of(2))

    def test_rounding_commutes_with_integer_translation_at_ties(self):
        for x in (QQ(-3)/2,QQ(-1)/2,QQ(1)/2,QQ(3)/2):
            for k in range(-3,4):
                self.assertEqual(nearest(x+k),nearest(x)+k)

    def test_exact_bernoulli_decisions(self):
        class Zero:
            def getrandbits(self,n): return 0
        class One:
            def getrandbits(self,n): return (1<<n)-1
        self.assertTrue(bernoulli_exp_pi(1,RandomBits(source=Zero())))
        self.assertFalse(bernoulli_exp_pi(1,RandomBits(source=One())))
        self.assertTrue(bernoulli_exp_pi(0,RandomBits(source=One())))

    def test_norm_one_distortion_real_and_complex(self):
        for weights,a in [((1,1),(QQ(1)/7,-QQ(1)/7)),
                          ((1,2),(QQ(1)/8,-QQ(1)/8)),
                          ((2,2),(QQ(1)/3,-QQ(1)/3))]:
            eta=QQ(1)/100000
            values=norm_one_exp_approx(a,weights,eta)
            product=QQ(1)
            for v,w in zip(values,weights): product*=v**w
            self.assertEqual(product,1)
            RB=RealBallField(1024)
            for v,ai,w in zip(values,a,weights):
                self.assertLess(endpoints(abs(RB(v)/(RB(ai)/w).exp()-1))[1],eta)

    def test_distortion_grid_satisfies_paper_bound(self):
        eps=QQ(1)/100
        delta=distortion_grid(3,108,eps)
        RB=RealBallField(512)
        bound=(RB(eps/40)**5/RB(9))/(RB(90).exp()*108*(3*RB(240/eps).log()).sqrt())
        self.assertLessEqual(delta,endpoints(bound)[0])

    def test_walk_steps_use_a_certified_volume_upper_bound(self):
        N=corollary65_steps(3,108,16,QQ(1)/100)
        RB=RealBallField(512)
        target=23+2*RB(200).log()+RB(108*16).log()
        self.assertGreaterEqual(N,endpoints(target)[1])
        self.assertLessEqual(N,endpoints(target)[0]+1)
        self.assertGreaterEqual(corollary65_steps(3,108,16,QQ(1)/10000),N)


class GaussianTests(unittest.TestCase):
    def test_membership_tail_and_replay_on_nonorthogonal_basis(self):
        B=matrix(QQ,[[1,0],[QQ(1)/2,1]])
        center=vector(QQ,[QQ(1)/3,-QQ(2)/7])
        for seed in range(6):
            a=discrete_gaussian(B,4,QQ(1)/100,center,RandomBits(seed))
            b=discrete_gaussian(B,4,QQ(1)/100,center,RandomBits(seed))
            self.assertEqual(a.value,b.value)
            self.assertEqual(vector(ZZ,a.coefficients)*B,a.value)
            squared=(a.value-center).dot_product(a.value-center)
            lower=endpoints(32*RealBallField(256)(800).log())[0]
            self.assertLessEqual(squared,lower)

    def test_rectangular_and_zero_rank(self):
        B=matrix(QQ,[[1,-1,0],[0,1,-1]])/100
        a=discrete_gaussian(B,1,QQ(1)/100,rng=RandomBits(1))
        self.assertEqual(sum(a.value),0)
        z=discrete_gaussian(matrix(QQ,0,1),1,QQ(1)/10)
        self.assertEqual(tuple(z.value),(0,))

    def test_reject_invalid_hypotheses(self):
        with self.assertRaises(ValueError):
            discrete_gaussian(identity_matrix(QQ,2),QQ(1)/100,QQ(1)/100)
        with self.assertRaises(ValueError):
            discrete_gaussian(matrix(QQ,[[1,1],[2,2]]),10,QQ(1)/100)
        with self.assertRaises(ValueError):
            discrete_gaussian(matrix(QQ,[[1,-1]]),10,QQ(1)/100,center=[1,1])


class LatticeTests(unittest.TestCase):
    def test_lll_unimodular_and_lovasz(self):
        G=matrix(AA,[[100,99],[99,100]])
        U=lll_transform(G)
        self.assertEqual(abs(U.det()),1)
        mu,norms=gram_schmidt(U*G*U.transpose())
        self.assertLessEqual(abs(mu[1,0]),QQ(1)/2)
        self.assertGreaterEqual(norms[1],(QQ(3)/4-mu[1,0]**2)*norms[0])

    def test_svp_and_hkz_against_independent_small_enumeration(self):
        G=matrix(AA,[[9,8,1],[8,9,2],[1,2,5]])
        candidates=[]
        for i in range(-3,4):
            for j in range(-3,4):
                for k in range(-3,4):
                    if i or j or k:
                        v=vector(ZZ,[i,j,k]); candidates.append(v*G*v)
        c,length=shortest_vector_coefficients(G)
        self.assertEqual(length,min(candidates))
        U=hkz_transform(G)
        H=U*G*U.transpose()
        self.assertEqual(abs(U.det()),1)
        self.assertEqual(H[0,0],length)
        for start in range(3):
            # Check the entire recursive HKZ condition independently through GSO.
            mu,norms=gram_schmidt(H)
            P=matrix(AA,3-start,3-start,lambda i,j:sum(
                mu[start+i,q]*mu[start+j,q]*norms[q]
                for q in range(start,min(start+i,start+j)+1)))
            _,shortest=shortest_vector_coefficients(P)
            self.assertEqual(shortest,P[0,0])


class GeometryTests(unittest.TestCase):
    def test_exact_uniform_mass_including_boundaries(self):
        # Exhaust ALL proposal grid points, rather than a statistical frequency test.
        sampler=UniformLatticeSampler(matrix(AA,[[2]]),ProductRegion([30],[]),
                                      [QQ(3)/2],limits=Limits(attempts=1))
        N=sampler.certificate.denominator
        e=sampler.epsilon
        outer=ZZ((N*(1+4*e)*30).floor())
        counts=Counter()
        class FixedProposal:
            def __init__(self,k): self.k=k
            def integer(self,lo,hi):
                assert lo<=self.k<=hi
                return self.k
        for k in range(-outer,outer+1):
            try:
                v,_=sampler.sample(FixedProposal(k))
                counts[ZZ(v[0])]+=1
            except ResourceLimit: pass
        expected={k for k in range(-20,21) if abs(2*k-QQ(3)/2)<=30}
        self.assertEqual(set(counts),expected)
        self.assertEqual(set(counts.values()),{2*N})

    def test_complex_disks_and_geometric_hypotheses(self):
        region=ProductRegion([],[2])
        self.assertTrue(region.contains([2,0]))
        self.assertFalse(region.contains([2,QQ(1)/100]))
        for seed in range(10):
            self.assertTrue(region.contains(region.sample_grid(4,RandomBits(seed))))
        with self.assertRaises(ValueError):
            UniformLatticeSampler(identity_matrix(AA,2),ProductRegion([1,1],[]),[0,0])


class NumberFieldTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.K=NumberField(t*t-2,'a')
        cls.ctx=NumberFieldContext(cls.K)

    def test_integral_basis_is_an_additive_basis(self):
        with self.assertRaises(ValueError):
            NumberFieldContext(self.K,[1,2*self.K.gen()])
        L=NumberField(t*t-5,'a')
        with self.assertRaises(ValueError):
            NumberFieldContext(L,[1,L.gen()])
        NumberFieldContext(L,L.integral_basis())

    def test_fractional_crt_uses_local_congruence(self):
        ctx=self.ctx
        ray=RayConditions(ctx,ctx.K.ideal(2),tau=1)
        ideal=ctx.K.ideal(QQ(3)/5)
        gamma=ctx.K(QQ(1)/3)
        shift=ray.affine_shift(ideal,gamma)
        self.assertIn(shift-gamma,ideal)
        self.assertTrue(ray.contains(shift))
        # A fractional shift can be locally congruent without lying in the integral modulus.
        self.assertNotIn(shift-1,ray.modulus)
        with self.assertRaises(EmptySupport): ray.affine_shift(ctx.K.ideal(1),QQ(1)/2)

    def test_algorithm1_nontrivial_ray_negative_scaling_and_fractional_shift(self):
        ctx=self.ctx
        ray=RayConditions(ctx,ctx.K.ideal(2),real_places=[0,1],tau=-1)
        prepared=Algorithm1(ctx,ctx.K.ideal(QQ(3)/5),gamma=QQ(1)/3,ray=ray,x=[-2,QQ(1)/2])
        for seed in range(5):
            sample=prepared.sample(RandomBits(seed))
            self.assertTrue(prepared.verify(sample.element))
            self.assertTrue(all(ctx.places[i](sample.element).real()<0 for i in (0,1)))
            self.assertTrue(sample.report['uniform_conditional_on_completion'])

    def test_algorithm1_mixed_and_complex_signatures(self):
        for f in (t*t+1,t*t*t-2):
            K=NumberField(f,'a'); ctx=NumberFieldContext(K)
            scales=[1]*ctx.r1+[(QQ(1)/2,QQ(1)/3)]*ctx.r2
            prepared=Algorithm1(ctx,K.ideal(3),x=scales)
            sample=prepared.sample(RandomBits(6))
            self.assertTrue(prepared.verify(sample.element))
            self.assertEqual(len(sample.point(ctx)),ctx.n)

    def test_algorithm2_replay_and_prime_ideal_identity(self):
        ctx=self.ctx
        ray=RayConditions(ctx,ctx.K.ideal(2),real_places=[0],tau=-1)
        args=dict(context=ctx,ideal=ctx.K.ideal(3),walk=WalkParameters(19,3),ray=ray,
                  require_mixing_bound=False)
        a=algorithm2(**args,rng=RandomBits(123))
        b=algorithm2(**args,rng=RandomBits(123))
        self.assertEqual(a.element,b.element)
        walk_ideal=ctx.K.ideal(3)
        for p in a.prime_factors:
            self.assertTrue(p.is_prime())
            self.assertLessEqual(p.norm(),19)
            self.assertTrue(ray.coprime(p))
            walk_ideal*=p
        self.assertEqual(walk_ideal,a.walk_ideal)
        self.assertIn(a.element,walk_ideal)
        self.assertTrue(ray.contains(a.element))
        self.assertEqual(ctx.K.ideal(a.element),a.output_ideal()*a.input_ideal)
        self.assertEqual(a.report['mixing_guarantee'],'not_established')

    def test_algorithm2_imaginary_quadratic_has_zero_dimensional_gaussian(self):
        K=NumberField(t*t+1,'a'); ctx=NumberFieldContext(K)
        a=algorithm2(ctx,K.ideal(1),WalkParameters(13,2),rng=RandomBits(9),
                     require_mixing_bound=False)
        self.assertEqual(a.log_distortion,(0,))
        self.assertEqual(a.distortion,(1,))

    def test_algorithm2_refuses_unjustified_mixing_claim(self):
        with self.assertRaises(UnresolvedBound):
            algorithm2(self.ctx,self.K.ideal(1),WalkParameters(7,2))
        with self.assertRaises(ValueError):
            WalkParameters(7,2,mixing_l1_bound=QQ(1)/1000)
        with self.assertRaises(ValueError):
            algorithm2(self.ctx,self.K.ideal(1),WalkParameters(7,2,subgroup_index=2),
                       require_mixing_bound=False)
        generated=WalkParameters.from_prime_bound(self.ctx,RayConditions(self.ctx),19)
        self.assertGreater(generated.steps,20)
        self.assertIsNone(generated.mixing_l1_bound)

    def test_algorithm2_proper_subgroup_oracle(self):
        # Norm modulo 5 modulo squares gives a ray-class character for modulus 5.
        ctx=self.ctx
        ray=RayConditions(ctx,ctx.K.ideal(5))
        def oracle(ideal):
            norm=QQ(ideal.norm())
            residue=(norm.numerator()*norm.denominator().inverse_mod(5))%5
            return residue in (1,4)
        a=algorithm2(ctx,ctx.K.ideal(3),WalkParameters(19,3,subgroup_index=2),
                     ray=ray,subgroup_oracle=oracle,rng=RandomBits(3),require_mixing_bound=False)
        self.assertTrue(all(oracle(p) and ray.coprime(p) for p in a.prime_factors))
        self.assertTrue(ray.contains(a.element))

    def test_core_does_not_call_pari_class_or_unit_groups(self):
        from sage.rings.number_field.number_field import NumberField_generic
        with patch.object(NumberField_generic,'pari_bnf',side_effect=AssertionError('unexpected class-group computation')):
            a=algorithm2(self.ctx,self.K.ideal(1),WalkParameters(11,2),
                         rng=RandomBits(17),require_mixing_bound=False)
            self.assertIn(a.element,self.K.ideal(1))

    def test_resource_failure_returns_no_fake_sample(self):
        class AlwaysOne:
            def integer(self,lo,hi): return 1
        sampler=PrimeIdealSampler(self.ctx,RayConditions(self.ctx),7,limits=Limits(attempts=2))
        with self.assertRaises(ResourceLimit): sampler.sample(AlwaysOne())

    def test_zero_is_excluded_even_with_trivial_ray(self):
        prepared=Algorithm1(self.ctx,self.K.ideal(1))
        self.assertFalse(prepared.verify(0))

    def test_prime_sampler_exact_mass_with_splitting_and_ramification(self):
        K=NumberField(t*t+1,'i'); ctx=NumberFieldContext(K)
        bound=13
        sampler=PrimeIdealSampler(ctx,RayConditions(ctx),bound,limits=Limits(attempts=1))
        masses=defaultdict(QQ)
        class Choices:
            def __init__(self,p,j,r): self.p,self.j,self.r=p,j,r
            def integer(self,lo,hi): return self.p
            def below(self,n): return self.j
            def rational_event(self,numerator,denominator): return self.r<numerator
        for p in range(2,bound+1):
            if not ZZ(p).is_prime(): continue
            eligible=sampler.eligible_above(p)
            k=len(eligible)
            for j in range(k):
                for r in range(ctx.n):
                    try:
                        q,_=sampler.sample(Choices(p,j,r))
                        masses[q]+=QQ(1)/(bound*k*ctx.n)
                    except ResourceLimit: pass
        self.assertGreater(len(masses),3)
        self.assertEqual(set(masses.values()),{QQ(1)/(ctx.n*bound)})


if __name__=='__main__': unittest.main()
