"""Algorithms 3–6 and independent completeness/compact-representation checks."""

import unittest
from unittest.mock import patch
from sage.all import (QQ,ZZ,RealBallField,PolynomialRing,NumberField,matrix,vector,identity_matrix)

from rigorous_cnt import NumberFieldContext,RandomBits,Limits,ResourceLimit,UnresolvedBound
from rigorous_cnt.factor_base import FactorBase,prime_ideals_up_to
from rigorous_cnt.compact import CompactElement
from rigorous_cnt.residue import dedekind_residue
from rigorous_cnt.postprocess import bkp_basis
from rigorous_cnt.relations import (SUnitSettings,RelationEngine,UniformResidueUnits,
    IdealRelation,SUnitObservation,algorithm3,algorithm4,algorithm5)
from rigorous_cnt.sunits import certify_class_generation,SUnitPostprocessor,algorithm6,_certify_ray_generation
from rigorous_cnt.numerics import endpoints


R=PolynomialRing(QQ,'t');t=R.gen()


class ResidueTests(unittest.TestCase):
    def test_known_residues_without_class_group_calls(self):
        from sage.rings.number_field.number_field import NumberField_generic
        for f in (t*t-2,t*t+1,t*t+5):
            K=NumberField(f,'a');ctx=NumberFieldContext(K)
            with patch.object(NumberField_generic,'pari_bnf',side_effect=AssertionError('bnf forbidden')):
                c=dedekind_residue(ctx)
            RB=RealBallField(512)
            if f==t*t-2:
                known=(1+RB(2).sqrt()).log()/RB(2).sqrt()
            elif f==t*t+1:
                known=RB.pi()/4
            else:
                known=2*RB.pi()/RB(20).sqrt()
            lo,hi=endpoints(known)
            self.assertLessEqual(c.rho_lower,lo)
            self.assertGreaterEqual(c.rho_upper,hi)
            self.assertLess(c.rho_upper/c.rho_lower,QQ(5)/4)


class CompactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.K=NumberField(t*t-2,'a');cls.ctx=NumberFieldContext(cls.K)

    def test_large_cancellation_and_modular_recovery(self):
        u=self.K.gen()+1;N=ZZ(10)**50
        c=CompactElement(self.ctx,[(u*u,N),(u,1-2*N)])
        with self.assertRaises(ResourceLimit):c.expand()
        self.assertEqual(c.materialize(),u)
        for x in c.logarithms(128):
            lo,hi=endpoints(x)
            self.assertLessEqual(hi-lo,QQ(2)**(-128))

    def test_fractional_reconstruction_and_torsion_phase(self):
        N=ZZ(10)**30
        c=CompactElement(self.ctx,[(self.K(4),N),(self.K(2),-2*N-3)])
        self.assertEqual(c.materialize(),QQ(1)/8)
        K=NumberField(t*t+1,'i');ctx=NumberFieldContext(K)
        c=CompactElement(ctx,[(K.gen(),N+1)])
        self.assertEqual(c.materialize(),K.gen())

    def test_reconstruction_with_non_power_integral_basis(self):
        K=NumberField(t*t-5,'a');ctx=NumberFieldContext(K);u=(1+K.gen())/2
        c=CompactElement(ctx,[(u*u,100000),(u,-199999)])
        self.assertEqual(c.materialize(),u)


class PostprocessingTests(unittest.TestCase):
    def test_bkp_unknown_rank_and_integer_span(self):
        def oracle(p):
            RB=RealBallField(p)
            return matrix(RB,[[2,2*RB(2).sqrt()],[3,3*RB(2).sqrt()],[0,0]])
        c=bkp_basis(oracle,QQ(1)/10,2,integer_columns=1)
        self.assertEqual(c.rank,1)
        self.assertEqual(abs((c.transform*vector(ZZ,[2,3,0]))[0]),1)
        self.assertTrue(c.report['rank_and_integer_span_certified'])

    def test_index_two_is_not_mistaken_for_completeness(self):
        K=NumberField(t*t-2,'a');a=K.gen();ctx=NumberFieldContext(K)
        fb=FactorBase(ctx,K.primes_above(2))
        post=SUnitPostprocessor(fb,dedekind_residue(ctx),certify_class_generation(fb))
        def observations(elements):
            return [SUnitObservation(CompactElement.of(ctx,e),tuple(e.valuation(P) for P in fb.primes),fb)
                    for e in elements]
        bad=post.check(observations([(a+1)**2,a]))
        good=post.check(observations([a+1,a]))
        self.assertFalse(bad.complete);self.assertEqual(bad.index_interval,(2,2))
        self.assertTrue(good.complete);self.assertEqual(good.index_interval,(1,1))

    def test_euclidean_slope_factor(self):
        K=NumberField(t*t+1,'i');ctx=NumberFieldContext(K);fb=FactorBase(ctx,K.primes_above(3))
        post=SUnitPostprocessor(fb,dedekind_residue(ctx),certify_class_generation(fb))
        c=post.check([SUnitObservation(CompactElement.of(ctx,3),(1,),fb)])
        self.assertTrue(c.complete)
        RB=RealBallField(512);expected=(1+RB(9).log()**2).sqrt()
        lo,hi=endpoints(expected)
        self.assertLessEqual(c.euclidean_volume[0],lo)
        self.assertGreaterEqual(c.euclidean_volume[1],hi)
        self.assertGreater(c.euclidean_volume[0],2)

    def test_class_generation_requires_actual_witnesses(self):
        K=NumberField(t*t+5,'a');a=K.gen();ctx=NumberFieldContext(K)
        P=K.ideal(2,a+1);Q=K.ideal(3,a+1);fb=FactorBase(ctx,[Q])
        with self.assertRaises(UnresolvedBound):certify_class_generation(fb)
        witness=IdealRelation(a+1,P,(1,),fb,0,{})
        self.assertTrue(witness.verify())
        cert=certify_class_generation(fb,[witness])
        post=SUnitPostprocessor(fb,dedekind_residue(ctx),cert)
        result=post.check([SUnitObservation(CompactElement.of(ctx,a-2),(2,),fb)])
        self.assertTrue(result.complete)

    def test_incorrect_observation_is_rejected(self):
        K=NumberField(t*t-2,'a');ctx=NumberFieldContext(K);fb=FactorBase(ctx,K.primes_above(2))
        post=SUnitPostprocessor(fb,dedekind_residue(ctx),certify_class_generation(fb))
        with self.assertRaises(ValueError):
            post.check([SUnitObservation(CompactElement.of(ctx,3),(0,),fb)])

    def test_ray_generation_is_distinct_from_class_generation(self):
        K=NumberField(t*t+1,'i');ctx=NumberFieldContext(K)
        engine=RelationEngine(ctx,SUnitSettings(modulus_override=K.ideal(3)))
        # O/(3) has eight units, while roots of unity supply only four.
        empty=FactorBase(ctx,[])
        empty_post=SUnitPostprocessor(empty,engine.residue,certify_class_generation(empty))
        check=empty_post.check([])
        self.assertTrue(check.complete)
        self.assertFalse(_certify_ray_generation(engine,check))
        # 1+i has residue order eight and supplies the missing ray classes.
        fb=FactorBase(ctx,[K.ideal(1+K.gen())])
        post=SUnitPostprocessor(fb,engine.residue,certify_class_generation(fb))
        check=post.check([SUnitObservation(CompactElement.of(ctx,1+K.gen()),(1,),fb)])
        self.assertTrue(check.complete)
        self.assertTrue(_certify_ray_generation(engine,check))


class RelationTests(unittest.TestCase):
    def test_uniform_residue_units_prime_power_and_crt(self):
        K=NumberField(t*t-2,'a');ctx=NumberFieldContext(K)
        sampler=UniformResidueUnits(ctx,K.ideal(12))
        for seed in range(20):
            a=sampler.sample(RandomBits(seed))
            self.assertTrue(all(a.valuation(P)==0 for P,e in K.ideal(12).factor()))
        # Enumerate the unique accepted representatives for the ramified quotient O/(2).
        sampler=UniformResidueUnits(ctx,K.ideal(2))
        class Bits:
            def __init__(self,a,b):self.values=iter([a,b])
            def integer(self,lo,hi):return next(self.values)
        outputs=[]
        for a,b in ((1,0),(1,1)):
            outputs.append(sampler.sample(Bits(a,b)))
        self.assertNotIn(outputs[0]-outputs[1],K.ideal(2))

    def test_algorithm3_and_exceptional_algorithm5(self):
        K=NumberField(t*t+1,'i');ctx=NumberFieldContext(K);P=K.prime_above(2)
        settings=SUnitSettings(smooth_bound=59,generation_radius=2,modulus_override=P,relation_attempts=3000)
        engine=RelationEngine(ctx,settings)
        fb=FactorBase(ctx,[Q for Q in prime_ideals_up_to(ctx,59) if Q!=P])
        relation=algorithm3(engine,K.ideal(3),fb,rng=RandomBits(11))
        self.assertTrue(relation.verify())
        exceptional=algorithm5(engine,P,fb,rng=RandomBits(4))
        self.assertTrue(exceptional.verify())
        self.assertEqual(exceptional.element.valuation(P),1)


class EndToEndTests(unittest.TestCase):
    def test_full_solver_without_existing_group_solver(self):
        from sage.rings.number_field.number_field import NumberField_generic
        K=NumberField(t*t+1,'i');ctx=NumberFieldContext(K);S=K.primes_above(3)
        settings=SUnitSettings(smooth_bound=59,generation_radius=2,max_relations=60)
        with patch.object(NumberField_generic,'pari_bnf',side_effect=AssertionError('bnf forbidden')):
            result=algorithm6(ctx,S,settings=settings,rng=RandomBits(7))
        self.assertTrue(result.verify());self.assertEqual(result.rank,1)
        self.assertEqual(result.torsion_order,4)
        materialized=result.materialize_generators()
        reference=K.S_unit_group(proof=True,S=S)
        coordinates=[reference.log(g) for g in materialized]
        self.assertEqual(abs(ZZ(coordinates[0][1])),1)
        self.assertFalse(result.report['paper_bit_complexity_certified'])


if __name__=='__main__':unittest.main()
