"""Theorem 13.2 Sample and Algorithms 3–5.

Relations are always checked exactly. Distribution/complexity hypotheses are
reported separately from these algebraic certificates.
"""

from dataclasses import dataclass,field
from sage.all import QQ,ZZ,AA,RealBallField,matrix,vector,prime_range
from .arithmetic import RayConditions
from .sampling import algorithm2,paper_radius,BoxSample,IdealSample
from .primes import PrimeIdealSampler
from .parameters import WalkParameters
from .gaussian import discrete_gaussian
from .numerics import rational,endpoints,precisions,dyadic_lower
from .runtime import Limits,RandomBits,ResourceLimit,UnresolvedBound
from .residue import dedekind_residue
from .factor_base import FactorBase,prime_ideals_up_to
from .compact import CompactElement


@dataclass(frozen=True)
class SUnitSettings:
    """Computational choices; unknown theorem bounds are never marked proved.

    Defaults give a reference run with exact result certification, but without
    an asserted mixing/generation-time guarantee. Supply walk_factory and its
    error certificates to request the corresponding sampling guarantee.
    """
    smooth_bound: int=59
    walk_prime_bound: int=19
    walk_steps: int=2
    generation_radius: object=QQ(2)
    generation_radius_justification: str="to be checked against the final generating basis"
    require_mixing_bound: bool=False
    walk_factory: object=field(default=None,repr=False,compare=False)
    modulus_override: object=field(default=None,repr=False,compare=False)
    max_relations: int=200
    relation_attempts: int=20000
    max_ideal_bits: int=65536
    max_minkowski_bound: int=10000
    max_factor_base_bound: int=10000
    check_every: int=1

    def __post_init__(self):
        for name in ('smooth_bound','walk_prime_bound','walk_steps','max_relations',
                     'relation_attempts','max_ideal_bits','max_minkowski_bound','max_factor_base_bound','check_every'):
            value=getattr(self,name)
            if ZZ(value)!=value or value<0:
                raise ValueError(f'{name} must be a nonnegative integer')
        if self.smooth_bound<2 or self.walk_prime_bound<2 or min(
                self.max_relations,self.relation_attempts,self.max_ideal_bits,self.check_every)<1:
            raise ValueError("invalid relation-generation budget")
        if rational(self.generation_radius)<=0:
            raise ValueError("generation_radius must be a positive rational")
        if self.require_mixing_bound and self.walk_factory is None:
            raise UnresolvedBound("certified mixing requires a walk_factory")


class UniformResidueUnits:
    """Uniform (O_K/m)^* using prime-power quotient boxes and exact CRT."""
    def __init__(self,context,modulus,limits=Limits()):
        self.context=context;self.modulus=context.K.ideal(modulus);self.limits=limits
        self.parts=[]
        O=context.integral_basis
        for P,e in self.modulus.factor():
            Q=P**e
            rows=matrix(QQ,[context.coordinates_in_basis(a,O) for a in Q.basis()])
            if any(x.denominator()!=1 for x in rows.list()):
                raise ArithmeticError("integral quotient had fractional integral-basis coordinates")
            H=matrix(ZZ,rows).hermite_form()
            diagonals=tuple(ZZ(H[i,i]) for i in range(context.n))
            if any(d<=0 for d in diagonals):
                raise ArithmeticError("invalid quotient HNF")
            complement=self.modulus/Q
            weight=context.K(1) if complement==context.unit_ideal else complement.element_1_mod(Q)
            self.parts.append((P,Q,diagonals,weight))

    def sample(self,rng):
        if not self.parts:
            return self.context.K(1)
        out=self.context.K(0)
        for P,Q,diagonals,weight in self.parts:
            for _ in range(self.limits.attempts):
                value=sum((rng.integer(0,d-1)*a for d,a in zip(diagonals,self.context.integral_basis)),self.context.K(0))
                if value not in P:
                    out+=weight*value
                    break
            else:
                raise ResourceLimit("residue-unit rejection budget exhausted")
        if any(out.valuation(P)!=0 for P,Q,ds,w in self.parts):
            raise ArithmeticError("CRT residue is not a unit modulo the modulus")
        return out


class SpecializedSampler:
    """Theorem 13.2: random tau in the finite residue unit group, no real ray."""
    def __init__(self,context,modulus,block_size,omega,settings,limits=Limits()):
        self.context=context;self.modulus=context.K.ideal(modulus)
        self.block_size=block_size;self.omega=omega;self.settings=settings;self.limits=limits
        self.radius=paper_radius(context,context.unit_ideal,context.scales(),self.modulus,block_size,omega)
        self.epsilon=dyadic_lower(lambda R:1/(1200*context.discriminant*R(self.radius)**context.n),limits)
        ray=RayConditions(context,self.modulus)
        if settings.walk_factory is None:
            self.walk=WalkParameters(settings.walk_prime_bound,settings.walk_steps)
        else:
            self.walk=settings.walk_factory(context,ray,self.epsilon)
        self.walk.check(self.epsilon,settings.require_mixing_bound)
        self.units=UniformResidueUnits(context,self.modulus,limits)
        self.primes=PrimeIdealSampler(context,ray,self.walk.prime_bound,limits=limits)
        self._normalizations={}

    def _normalize(self,ideal):
        if ideal not in self._normalizations:
            entries=[q for a in ideal.basis() for q in self.context.coordinates_in_basis(a,self.context.integral_basis)]
            denominator=ZZ(1)
            for q in entries:denominator=denominator.lcm(q.denominator())
            content=ZZ(0)
            for q in entries:content=content.gcd(ZZ(q*denominator))
            scalar=QQ(abs(content))/denominator
            if any(self.context.K(scalar).valuation(P)!=0 for P,e in self.modulus.factor()):
                scalar=QQ(1)
            self._normalizations[ideal]=(scalar,ideal/scalar)
        return self._normalizations[ideal]

    def sample(self,ideal,y=None,rng=None):
        rng=rng or RandomBits()
        tau=self.units.sample(rng)
        scalar,normalized=self._normalize(self.context.K.ideal(ideal))
        ray=RayConditions(self.context,self.modulus,tau=tau/scalar)
        original_y=self.context.scales(y)
        common=abs(original_y[0][0]) or abs(original_y[0][1])
        normalized_y=tuple((re/common,im/common) for re,im in original_y)
        result=algorithm2(self.context,normalized,self.walk,ray=ray,y=normalized_y,epsilon=self.epsilon,
            block_size=self.block_size,omega=self.omega,rng=rng,limits=self.limits,
            require_mixing_bound=self.settings.require_mixing_bound,prime_sampler=self.primes)
        if scalar==1 and common==1:return result
        # Positive rational scaling is an exact change of variables in every
        # lattice and box. It preserves the output law, including the local ray.
        scale=tuple((common*re,common*im) for re,im in result.box_sample.scale)
        box=BoxSample(scalar*result.box_sample.element,scale,
            scalar*common*result.box_sample.radius,result.box_sample.coefficients,
            dict(result.box_sample.report,rational_normalization=str(scalar),ambient_normalization=str(common)))
        return IdealSample(scalar*result.element,self.context.K.ideal(ideal),scalar*result.walk_ideal,
            result.prime_factors,result.log_distortion,result.distortion,result.grid_spacing,box,
            dict(result.report,rational_normalization=str(scalar)))


def _x_interval(context,R):
    logD=R(context.discriminant).log()
    a=(logD/logD.log()**2)**(QQ(2)/3)
    b=(R(context.n)/R(context.n).log())**(QQ(2)/3)
    alo,ahi=endpoints(a);blo,bhi=endpoints(b)
    return R(max(alo,blo)).union(R(max(ahi,bhi)))


class RelationEngine:
    """Shared field parameters, residue certificate, and sampling subroutines."""
    def __init__(self,context,settings=None,residue=None,limits=Limits()):
        self.context=context;self.settings=settings or SUnitSettings();self.limits=limits
        self.residue=residue or dedekind_residue(context)
        self.residue.validate_for(context)
        if self.residue.rho_upper>=2*self.residue.rho_lower:
            raise ValueError("Algorithm 3 needs a factor-two residue approximation")
        self.block_size=max(2,min(context.n,int((AA(context.n)**(QQ(2)/3)).ceil())))
        self.rho_approx=self.residue.rho_lower
        chosen=None
        for p in precisions(limits):
            R=RealBallField(p);x=_x_interval(context,R)
            threshold=(x*x.log()**2).exp()
            lo,hi=endpoints(threshold)
            if self.rho_approx<=lo:
                chosen=context.unit_ideal;self.modulus_branch='small_residue';break
            if self.rho_approx>hi:
                xlo,xhi=endpoints(x)
                if xhi.ceil()>self.settings.max_factor_base_bound:
                    raise ResourceLimit("modulus selection exceeds the prime enumeration budget")
                factors=[];ambiguous=False
                for P in prime_ideals_up_to(context,xhi.ceil()):
                    if P.norm()<xlo: factors.append(P)
                    elif P.norm()<=xhi: ambiguous=True;break
                if ambiguous: continue
                chosen=context.unit_ideal
                for P in factors:chosen*=P
                self.modulus_branch='large_residue';break
        self.modulus=chosen
        if self.settings.modulus_override is not None:
            self.modulus=context.K.ideal(self.settings.modulus_override)
            if not self.modulus or not self.modulus.is_integral():
                raise ValueError("modulus_override must be a nonzero integral ideal")
            if any(e!=1 for P,e in self.modulus.factor()):
                raise ValueError("Algorithm 5 requires a squarefree modulus")
            self.modulus_branch='explicit_override'
        self.exceptional_primes=tuple(P for P,e in self.modulus.factor())
        R=RealBallField(192);x=_x_interval(context,R)
        self.x_upper=endpoints(x)[1]
        r0=paper_radius(context,context.unit_ideal,context.scales(),self.modulus,self.block_size,1)
        target=max(QQ(self.settings.smooth_bound),QQ(self.settings.walk_prime_bound),10*self.x_upper**2)
        self.omega=1
        while endpoints(R(self.omega*r0)**context.n/R(context.n).exp())[0]<target:
            self.omega+=1
            if self.omega>limits.lattice_steps:
                raise ResourceLimit("radius parameter search exceeded its budget")
        self._samplers={}
        self._density_checks={}
        self.statistics={'sampling_calls':0,'relations':0,'smoothness_rejections':0}

    def sampler(self,modulus=None):
        modulus=self.modulus if modulus is None else self.context.K.ideal(modulus)
        if modulus not in self._samplers:
            self._samplers[modulus]=SpecializedSampler(self.context,modulus,self.block_size,
                                                     self.omega,self.settings,self.limits)
        return self._samplers[modulus]

    def smoothness_diagnostics(self,factor_base,sampler,mixing):
        key=(factor_base,sampler.modulus)
        if key not in self._density_checks:
            R=RealBallField(192)
            B=max(self.settings.smooth_bound,int(sampler.walk.prime_bound))
            usable=[P for P in factor_base.primes if P.norm()<=B and sampler.units.modulus.valuation(P)==0]
            count_condition=QQ(len(usable))>=endpoints(R(B)/(4*R(B).log()))[1]
            radius_condition=endpoints(R(sampler.radius)**self.context.n/R(self.context.n).exp())[0]>=B
            ray=RayConditions(self.context,sampler.modulus)
            if sampler.walk.prime_bound>self.settings.max_factor_base_bound:
                absorbed=False
            else:
                walk_primes=prime_ideals_up_to(self.context,sampler.walk.prime_bound)
                absorbed=all(P in factor_base.index for P in walk_primes if ray.coprime(P))
            probability=None
            if count_condition and radius_condition and absorbed:
                u=self.context.n*R(sampler.radius).log()/R(B).log()
                log_density=((1-u)*(4*R(B).log()).log()-R(self.residue.rho_upper).log()
                             -R(B).log()-u*u.log())
                modulus_factor=QQ(1)
                for P,e in sampler.modulus.factor():modulus_factor*=QQ(P.norm())/(P.norm()-1)
                probability=endpoints(R(modulus_factor)*log_density.exp()/6)[0]
            self._density_checks[key]={
                'prime_count_condition_verified':bool(count_condition),
                'radius_condition_verified':bool(radius_condition),
                'walk_factors_absorbed':bool(absorbed),
                'conditional_success_probability_lower_bound':None if probability is None else str(probability),
            }
        result=dict(self._density_checks[key])
        result['success_probability_status']=(
            'verified_under_GRH' if mixing=='verified_spectral_bound' else 'conditional_on_supplied_mixing') if (
            result['conditional_success_probability_lower_bound'] is not None and mixing!='not_established') else 'not_established'
        return result


@dataclass(frozen=True)
class IdealRelation:
    element: object
    input_ideal: object
    valuations: tuple
    factor_base: object
    attempts: int
    sample_report: dict

    def verify(self):
        K=self.factor_base.context.K
        return (self.element in self.input_ideal and all(e>=0 for e in self.valuations)
                and K.ideal(self.element)==self.input_ideal*self.factor_base.ideal(self.valuations))


def algorithm3(engine,ideal,factor_base,y=None,rng=None,modulus=None):
    """Find (alpha)=ideal*product(P^v), v>=0, by repeated specialized sampling."""
    context=engine.context
    if factor_base.context is not context:
        raise ValueError("factor-base context mismatch")
    ideal=context.K.ideal(ideal);rng=rng or RandomBits()
    sampler=engine.sampler(modulus)
    if not ideal or not RayConditions(context,sampler.modulus).coprime(ideal):
        raise ValueError("Algorithm 3 needs a nonzero ideal coprime to the modulus")
    for attempt in range(1,engine.settings.relation_attempts+1):
        sampled=sampler.sample(ideal,y,rng)
        engine.statistics['sampling_calls']+=1
        quotient=context.K.ideal(sampled.element)/ideal
        values=factor_base.factor(quotient)
        if values is None:
            engine.statistics['smoothness_rejections']+=1
            continue
        report=dict(sampled.report)
        report.update(engine.smoothness_diagnostics(factor_base,sampler,report['mixing_guarantee']))
        relation=IdealRelation(sampled.element,ideal,values,factor_base,attempt,report)
        if not relation.verify():
            raise ArithmeticError("Algorithm 3 relation failed exact reconstruction")
        engine.statistics['relations']+=1
        return relation
    raise ResourceLimit("Algorithm 3 exhausted its relation search budget")


def _exp_divisor_approx(values,weights,N,limits):
    for p in precisions(limits):
        R=RealBallField(p)
        targets=[(R(a)/w).exp() for a,w in zip(values,weights)]
        lower=min(endpoints(t)[0] for t in targets)
        if lower<=0:continue
        mid=tuple(t.mid().exact_rational() for t in targets)
        if all(endpoints(abs(R(y)-z))[1]<=lower/N for y,z in zip(mid,targets)):
            return mid


@dataclass(frozen=True)
class SUnitObservation:
    element: object
    valuations: tuple
    factor_base: object
    report: dict=field(default_factory=dict)

    def logarithmic_row(self,precision):
        R=RealBallField(precision)
        return vector(R,[-v for v in self.valuations]+list(self.element.logarithms(precision)))


def algorithm4(engine,factor_base,rng=None):
    """Gaussian S-divisor input, followed by Algorithm 3, as in §18.5."""
    context=engine.context;rng=rng or RandomBits()
    if factor_base.context is not context:
        raise ValueError("factor-base context mismatch")
    if any(P in factor_base.index for P in engine.exceptional_primes):
        raise ValueError("Algorithm 4's factor base must omit the modulus primes")
    m=len(context.places);d=len(factor_base)+m
    R=RealBallField(192)
    width=3*max(endpoints(R(d).log().sqrt())[1],rational(engine.settings.generation_radius))
    sigma=QQ(width.ceil())
    grid=R(engine.omega)**context.n*(11+16*R(context.discriminant).log()**2+9*context.n**2).exp()
    N=ZZ(endpoints(grid)[1].ceil())
    B=matrix.diagonal(QQ,[1]*len(factor_base)+[QQ(1)/N]*m)
    gaussian=discrete_gaussian(B,sigma,QQ(1)/100,rng=rng,limits=engine.limits)
    finite=tuple(ZZ(x) for x in gaussian.value[:len(factor_base)])
    infinite=tuple(QQ(x) for x in gaussian.value[len(factor_base):])
    size_bound=sum(abs(e)*ZZ(P.norm()).nbits() for P,e in zip(factor_base.primes,finite))
    if size_bound>engine.settings.max_ideal_bits:
        raise ResourceLimit("Algorithm 4 input ideal exceeds the coefficient budget")
    y=_exp_divisor_approx(infinite,context.weights,N,engine.limits)
    ideal=factor_base.ideal(finite)
    relation=algorithm3(engine,ideal,factor_base,y=y,rng=rng)
    valuations=tuple(a+v for a,v in zip(finite,relation.valuations))
    if context.K.ideal(relation.element)!=factor_base.ideal(valuations):
        raise ArithmeticError("Algorithm 4 output is not the asserted S-unit")
    return SUnitObservation(CompactElement.of(context,relation.element),valuations,factor_base,{
        'algorithm':4,'exact_principal_ideal_verified':True,
        'gaussian_sd_bound':'1/100','sigma':str(sigma),'divisor_grid_bits':int(N.nbits()),
        'input_valuations':tuple(int(x) for x in finite),'relation_attempts':relation.attempts,
        'mixing_guarantee':relation.sample_report['mixing_guarantee'],
        'success_probability_status':relation.sample_report['success_probability_status'],
        'generation_radius_bound':'not_yet_verified',
    })


def algorithm5(engine,prime,factor_base,rng=None):
    """Find (alpha)=q*product(P^v) with v_q(alpha)=1 for q|m0."""
    prime=engine.context.K.ideal(prime)
    if prime not in engine.exceptional_primes or prime in factor_base.index:
        raise ValueError("Algorithm 5 requires a modulus prime outside the factor base")
    relation=algorithm3(engine,prime,factor_base,rng=rng,modulus=engine.modulus/prime)
    if relation.element.valuation(prime)!=1:
        raise ArithmeticError("exceptional S-unit did not have valuation one")
    return relation
