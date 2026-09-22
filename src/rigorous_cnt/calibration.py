"""Optional small-field spectral calibration of Algorithm 2's random walk.

This is explicit, potentially expensive preprocessing using certified PARI class
and unit data. It is not the paper's asymptotically efficient parameter selector.
The sampler itself continues to execute Algorithms 1–2.
"""

from dataclasses import dataclass,field
from itertools import product
from math import prod
from sage.all import (QQ,ZZ,AA,RealBallField,ComplexBallField,matrix,vector,
                      identity_matrix,diagonal_matrix,prime_range,pari)
from .runtime import Limits,ResourceLimit
from .numerics import endpoints,rational,precisions
from .parameters import WalkParameters


def _ray_key(ray):
    return (tuple(tuple(QQ(x) for x in a.vector()) for a in ray.modulus.basis()),ray.real_places)


def _log_vector(context,a,R):
    if not a:
        raise ValueError("logarithmic embedding requires a nonzero element")
    values=[]
    for emb,w in zip(context.places,context.weights):
        z=emb(a)
        norm2=AA(z.real())**2+AA(z.imag())**2
        v=R(norm2)
        if not v > 0:
            raise ArithmeticError("embedding precision does not separate zero")
        values.append(v.log()*(QQ(w)/2))
    return vector(R,values)


def _centered_log(context,a,R):
    v=_log_vector(context,a,R)
    return v-vector(R,[sum(v)/len(v)]*len(v))


def _pari_element(context,nf,value):
    if value.type()=='t_COL':
        value=nf.nfbasistoalg(value)
    return context.K(value)


def _pari_ideal(context,nf,value):
    H=nf.idealhnf(value)
    return context.K.ideal([_pari_element(context,nf,c) for c in H])


def _pari_hnf(nf,ideal):
    """Convert through algebraic elements, not an assumed common integral basis."""
    columns=[list(nf.nfalgtobasis(a)) for a in ideal.basis()]
    return nf.idealhnf(pari(matrix(QQ,columns).transpose()))


def _arch_flags(context,ray,nf):
    if not ray.real_places:
        return [0]*context.r1
    if ray.real_places==tuple(range(context.r1)):
        return [1]*context.r1
    exact=[AA(e(context.K.gen()).real()) for e in context.places[:context.r1]]
    matched=[]
    for z in pari('f -> f.roots')(nf):
        if z.imag()!=0:
            continue
        numerical=z.real().sage()
        q=numerical.exact_rational() if hasattr(numerical,'exact_rational') else QQ(numerical)
        i=min(range(context.r1),key=lambda j:abs(exact[j]-q))
        separation=min((abs(exact[i]-v) for j,v in enumerate(exact) if j!=i),default=AA(1))
        if abs(exact[i]-q)>=separation/4 or i in matched:
            raise ResourceLimit("PARI root precision does not distinguish the requested real places")
        matched.append(i)
    if len(matched)!=context.r1:
        raise ArithmeticError("PARI real-place count differs from the number-field signature")
    return [int(i in ray.real_places) for i in matched]


def _prepare_group(context,ray,max_group_order):
    K=context.K
    # Keeping the primitive element unchanged makes the exact PARI/Sage maps explicit.
    if list(K.pari_polynomial().Vec()) != list(reversed(K.defining_polynomial().list())):
        raise ValueError("spectral calibration requires an unchanged integral PARI defining polynomial")
    bnf=K.pari_bnf(proof=True)
    nf=pari('g -> g.nf')(bnf)
    modulus=[_pari_hnf(nf,ray.modulus),_arch_flags(context,ray,nf)]
    bnr=bnf.bnrinit(modulus,1)
    cyc=tuple(int(x) for x in pari('g -> g.cyc')(bnr))
    if prod(cyc)>max_group_order:
        raise ResourceLimit("ray class group exceeds the calibration state limit")
    generators=tuple(_pari_ideal(context,nf,g) for g in pari('g -> g.gen')(bnr))
    if len(generators)!=len(cyc):
        raise ArithmeticError("inconsistent ray class group data")
    unit_group=K.unit_group(proof=True)
    units=tuple(unit_group.gens_values())
    orders=tuple(unit_group.gens_orders())
    rank=len(context.places)-1
    if len(units)!=rank+1 or not orders[0] or any(orders[1:]):
        raise ArithmeticError("unexpected fundamental unit representation")
    finite_orders=[]
    columns=[[] for _ in units]
    if ray.modulus!=K.ideal(1):
        bid=nf.idealstar(_pari_hnf(nf,ray.modulus),2)
        finite_orders=[int(x) for x in pari('g -> g.cyc')(bid)]
        for j,u in enumerate(units):
            columns[j]=[ZZ(x) for x in nf.ideallog(u,bid)]
    finite_orders.extend([2]*len(ray.real_places))
    for j,u in enumerate(units):
        columns[j].extend(ZZ(context.places[i](u).real()<0) for i in ray.real_places)
    if finite_orders:
        E=matrix(ZZ,len(finite_orders),len(units),lambda i,j:columns[j][i])
        augmented=E.augment(-diagonal_matrix(ZZ,finite_orders))
        kernel=augmented.right_kernel_matrix().matrix_from_columns(range(len(units)))
        free=kernel.matrix_from_columns(range(1,len(units)))
        if rank:
            H,U=free.hermite_form(transformation=True,include_zero_rows=False)
            W=U*kernel
            if H.nrows()!=rank or H.rank()!=rank:
                raise ArithmeticError("ray-unit kernel did not have full free rank")
        else:
            W=matrix(ZZ,0,len(units))
        for row in W.rows():
            image=E*row
            if any(image[i]%order for i,order in enumerate(finite_orders)):
                raise ArithmeticError("ray-unit exponent failed exact residue validation")
    else:
        W=matrix(ZZ,rank,len(units),lambda i,j:int(j==i+1))
    principal_generators=[]
    for g,d in zip(generators,cyc):
        exponents,alpha=bnr.bnrisprincipal(_pari_hnf(nf,g**d),1)
        alpha=_pari_element(context,nf,alpha)
        if any(exponents) or K.ideal(alpha)!=g**d:
            raise ArithmeticError("ray class group relation did not reconstruct exactly")
        difference=alpha-1
        if difference and any(difference.valuation(p)<e for p,e in ray.factorization):
            raise ArithmeticError("ray class group relation failed local congruence")
        if any(context.places[i](alpha).real()<=0 for i in ray.real_places):
            raise ArithmeticError("ray class group relation failed positivity")
        principal_generators.append(alpha)
    return nf,bnr,cyc,generators,units,W,tuple(principal_generators)


def _prime_data(context,ray,bound,nf,bnr,cyc,generators,cache):
    out=[]
    for p in prime_range(bound+1):
        for P in context.K.primes_above(p):
            if P.norm()>bound or not ray.coprime(P):
                continue
            if P not in cache:
                exponents,alpha=bnr.bnrisprincipal(_pari_hnf(nf,P),1)
                j=tuple(ZZ(x) for x in exponents)
                alpha=_pari_element(context,nf,alpha)
                reconstructed=context.K.ideal(alpha)
                for g,e in zip(generators,j): reconstructed*=g**e
                if len(j)!=len(cyc) or reconstructed!=P:
                    raise ArithmeticError("prime ray relation failed exact reconstruction")
                # bnrisprincipal promises a principal ray factor; check it too.
                delta=alpha-1
                if delta and any(delta.valuation(q)<e for q,e in ray.factorization):
                    raise ArithmeticError("prime ray relation failed local congruence")
                if any(context.places[i](alpha).real()<=0 for i in ray.real_places):
                    raise ArithmeticError("prime ray relation failed positivity")
                cache[P]=(P,j,alpha)
            out.append(cache[P])
    return out


@dataclass(frozen=True)
class SpectralMixingCertificate:
    prime_bound: int
    steps: int
    l1_upper: object
    spectral_upper: object
    cutoff_radius: object
    unit_rank: int
    volume_upper: object
    low_characters: int
    prime_count: int
    cyclic_orders: tuple
    precision: int
    _context: object=field(repr=False,compare=False)
    _ray: tuple=field(repr=False)

    def validate_for(self,context,ray,walk,epsilon):
        if (context is not self._context or _ray_key(ray)!=self._ray
                or walk.prime_bound!=self.prime_bound or walk.steps!=self.steps
                or walk.subgroup_index!=1):
            raise ValueError("mixing certificate does not match this field, ray, or random walk")
        if rational(walk.mixing_l1_bound)!=self.l1_upper or self.l1_upper>epsilon/2:
            raise ValueError("mixing certificate does not meet the requested error budget")


def calibrate_walk(context,ray=None,epsilon=QQ(1)/100,initial_prime_bound=19,
                   max_prime_bound=4096,max_group_order=64,max_frequency_vectors=50000,
                   max_characters=100000,limits=Limits()):
    """Verify a finite set of Hecke eigenvalues and bound the Gaussian tail.

    Uses Appendix A.1's Fourier bound, with actual ball-enclosed eigenvalues in
    place of its asymptotic analytic estimate. Class/unit computation is explicit
    preprocessing; this routine makes no claim to the paper's complexity bound.
    Only G=Pic^0_m is supported by this calibrator; the sampler accepts general G.
    """
    from .arithmetic import RayConditions
    ray=ray or RayConditions(context)
    eps=rational(epsilon)
    if ray.context is not context or not 0<eps<min(QQ(1),QQ(context.n)/20):
        raise ValueError("invalid context or Algorithm 2 error parameter")
    if not 2<=initial_prime_bound<=max_prime_bound:
        raise ValueError("invalid prime-bound range")
    if min(max_group_order,max_frequency_vectors,max_characters)<1:
        raise ValueError("calibration budgets must be positive")
    nf,bnr,cyc,generators,units,W,relations=_prepare_group(context,ray,max_group_order)
    n,m=context.n,len(context.places)
    rank=m-1
    volume=ZZ(context.discriminant*ray.modulus.norm()*ZZ(2)**len(ray.real_places))
    s=QQ(1)/(n*n)
    sprime=QQ(1)/(2800*n*n)
    walk_error=eps/2
    # Tail: 2 V s'^(-r) exp(-2 R^2 s^2) <= walk_error^2/2.
    for precision in precisions(limits):
        RB=RealBallField(precision)
        try:
            log_units=matrix(RB,[_log_vector(context,u,RB) for u in units])
            L=W*log_units
            if rank:
                dual=L.transpose()*(L*L.transpose()).inverse()
                target=(RB(4*volume/sprime**rank/(walk_error**2)).log()/(2*s*s)).sqrt()
                radius=QQ(max(1,endpoints(target)[1].ceil(),
                              endpoints((RB(rank)/(2*s*s)).sqrt())[1].ceil()))
                coefficient_bounds=[int(endpoints(sum(x*x for x in row).sqrt()*radius)[1].ceil())
                                    for row in L.rows()]
            else:
                dual=matrix(RB,m,0)
                radius=QQ(0)
                coefficient_bounds=[]
            break
        except (ZeroDivisionError,ArithmeticError):
            continue
    number_vectors=prod(2*k+1 for k in coefficient_bounds)
    if number_vectors>max_frequency_vectors:
        raise ResourceLimit(f"spectral calibration needs {number_vectors} frequency candidates")
    h=prod(cyc)
    if number_vectors*h>max_characters:
        raise ResourceLimit("spectral calibration exceeds the character budget")
    RB=RealBallField(precision)
    CB=ComplexBallField(precision)
    relation_logs=[_centered_log(context,a,RB) for a in relations]
    frequency_data=[]
    for k in product(*(range(-b,b+1) for b in coefficient_bounds)):
        u=dual*vector(ZZ,k)
        if rank and endpoints(sum(x*x for x in u))[0]>radius**2:
            continue
        base=tuple(u.dot_product(v)/d for v,d in zip(relation_logs,cyc))
        for finite in product(*(range(d) for d in cyc)):
            if not any(k) and not any(finite):
                continue
            phase=tuple(x+RB(QQ(j)/d) for x,j,d in zip(base,finite,cyc))
            frequency_data.append((u,phase))
    cache={}
    bound=int(initial_prime_bound)
    while bound<=max_prime_bound:
        primes=_prime_data(context,ray,bound,nf,bnr,cyc,generators,cache)
        if primes:
            prime_logs=[_centered_log(context,a,RB) for P,j,a in primes]
            spectral=QQ(0)
            for u,phase in frequency_data:
                re,im=RB(0),RB(0)
                for (P,j,a),v in zip(primes,prime_logs):
                    argument=u.dot_product(v)+sum(e*z for e,z in zip(j,phase))
                    angle=2*RB.pi()*argument
                    re+=angle.cos(); im+=angle.sin()
                upper=endpoints(CB(re/len(primes),im/len(primes)).abs())[1]
                spectral=max(spectral,upper)
                if spectral>QQ(3)/4:
                    break
            if spectral<=QQ(3)/4:
                if not frequency_data:
                    steps,l1=0,QQ(0)
                else:
                    factor=QQ(2*volume)/sprime**rank
                    if spectral==0:
                        steps=1
                    else:
                        quotient=RB(2*factor/(walk_error**2)).log()/(-2*RB(spectral).log())
                        steps=max(1,int(endpoints(quotient)[1].ceil()))
                    # Directly check the combined low-frequency and tail bound.
                    tail=RB(0) if not rank else (-2*RB(radius)**2*s*s).exp()
                    for _ in range(limits.lattice_steps):
                        total=(RB(factor)*(RB(spectral)**(2*steps)+tail)).sqrt()
                        l1=endpoints(total)[1]
                        if l1<=walk_error:
                            break
                        steps+=1
                    else:
                        raise ResourceLimit("spectral mixing step search exceeded its budget")
                certificate=SpectralMixingCertificate(bound,steps,l1,spectral,radius,rank,volume,
                    len(frequency_data),len(primes),cyc,precision,context,_ray_key(ray))
                return WalkParameters(bound,steps,1,l1,
                    "Ball-verified finite Hecke spectrum plus Appendix A.1 Gaussian tail; "
                    "certified PARI class/unit preprocessing; no paper complexity claim.",certificate)
        if bound==max_prime_bound:
            break
        bound=min(max_prime_bound,2*bound)
    raise ResourceLimit("no certified spectral gap found within the prime-bound budget")
