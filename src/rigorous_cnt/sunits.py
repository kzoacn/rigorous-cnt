"""Algorithm 6, complete-generation certificates, and restriction to arbitrary S."""

from dataclasses import dataclass,field
from math import factorial
from sage.all import QQ,ZZ,RealBallField,matrix,vector,diagonal_matrix,pari
from .factor_base import FactorBase,prime_ideals_up_to
from .compact import CompactElement,combine_compact
from .relations import RelationEngine,SUnitSettings,SUnitObservation,algorithm3,algorithm4,algorithm5
from .postprocess import bkp_basis,kessler_lower_bound,gram_volume
from .numerics import endpoints
from .runtime import Limits,RandomBits,ResourceLimit,UnresolvedBound


def minkowski_prime_bound(context):
    R=RealBallField(192)
    value=(4/R.pi())**context.r2*QQ(factorial(context.n))/ZZ(context.n)**context.n*R(context.discriminant).sqrt()
    return max(1,int(endpoints(value)[1].floor()))


@dataclass(frozen=True)
class ClassGenerationCertificate:
    factor_base: object
    minkowski_bound: int
    witnesses: tuple

    def verify(self):
        base=self.factor_base
        for P in prime_ideals_up_to(base.context,self.minkowski_bound):
            if P in base.index:continue
            if not any(r.input_ideal==P and r.factor_base is base and r.verify() for r in self.witnesses):
                return False
        return True


def certify_class_generation(factor_base,witnesses=(),max_bound=10000):
    bound=minkowski_prime_bound(factor_base.context)
    if bound>max_bound:
        raise ResourceLimit("Minkowski class-generation certificate exceeds its enumeration budget")
    cert=ClassGenerationCertificate(factor_base,bound,tuple(witnesses))
    if not cert.verify():
        raise UnresolvedBound("factor base does not yet account for all Minkowski-bound prime classes")
    return cert


def _certify_ray_generation(engine,ordinary_check):
    """T spans the ray group iff T-units surject onto the finite residue units.

    Ordinary class generation is already certified. The kernel of ray class ->
    class is therefore accounted for by the residue images of the T-units.
    """
    if engine.modulus==engine.context.unit_ideal:return True
    nf=engine.context.K.pari_nf()
    bid=nf.idealstar(engine.modulus.pari_hnf(),2)
    orders=tuple(ZZ(x) for x in pari('g -> g.cyc')(bid))
    if not orders:return True
    images=[]
    for element in ordinary_check.basis_elements:
        value=vector(ZZ,[0]*len(orders))
        for a,e in element.factors:
            value+=e*vector(ZZ,list(nf.ideallog(a,bid)))
        images.append([x%order for x,order in zip(value,orders)])
    images.append(list(nf.ideallog(engine.residue.torsion_generator,bid)))
    rows=matrix(ZZ,matrix(ZZ,images).stack(diagonal_matrix(ZZ,orders)),sparse=False)
    H=rows.hermite_form(include_zero_rows=False)
    return abs(H.det())==1


@dataclass(frozen=True)
class GenerationCheck:
    complete: bool
    rank: int
    expected_rank: int
    transform: object
    basis_elements: tuple
    basis_valuations: object
    euclidean_volume: object
    target_euclidean_volume: object
    index_interval: object
    precision: int
    report: dict


class SUnitPostprocessor:
    """BKP basis extraction followed by an independently bounded index test."""
    def __init__(self,factor_base,residue,class_certificate,limits=Limits()):
        self.factor_base=factor_base;self.context=factor_base.context
        residue.validate_for(self.context)
        if class_certificate.factor_base is not factor_base or not class_certificate.verify():
            raise ValueError("a matching verified class-generation certificate is required")
        self.residue=residue;self.class_certificate=class_certificate;self.limits=limits
        self.minimum=kessler_lower_bound(self.context,limits)
        self._verified={};self._logs={}

    def _validate(self,observation):
        if observation.factor_base is not self.factor_base or observation.element.context is not self.context:
            raise ValueError("observation belongs to a different factor base or context")
        if len(observation.valuations)!=len(self.factor_base):
            raise ValueError("observation valuation dimension mismatch")
        key=(observation.element.factors,tuple(observation.valuations))
        if key in self._verified:return
        values=[ZZ(0)]*len(self.factor_base)
        for a,e in observation.element.factors:
            factored=self.factor_base.factor(self.context.K.ideal(a),nonnegative=False)
            if factored is None:
                raise ValueError("a power-product base is not supported on this factor base")
            values=[v+e*w for v,w in zip(values,factored)]
        if tuple(values)!=tuple(observation.valuations):
            raise ValueError("observation's exact ideal factorization is incorrect")
        self._verified[key]=True

    def _row(self,observation,p):
        key=(observation.element.factors,p)
        if key not in self._logs:
            self._logs[key]=tuple(observation.element.logarithms(p,self.limits))
        R=RealBallField(p)
        return vector(R,[-e for e in observation.valuations]+list(self._logs[key]))

    def check(self,observations):
        observations=tuple(observations)
        for observation in observations:self._validate(observation)
        d=len(self.factor_base)+len(self.context.places)
        expected=d-1
        def oracle(p):
            return matrix(RealBallField(p),len(observations),d,
                          [x for o in observations for x in self._row(o,p)])
        reduced=bkp_basis(oracle,self.minimum,max(1,expected),len(self.factor_base),self.limits)
        M=reduced.transform
        elements=tuple(combine_compact(self.context,[o.element for o in observations],row) for row in M.rows())
        original_values=matrix(ZZ,len(observations),len(self.factor_base),
                               [x for o in observations for x in o.valuations])
        values=M*original_values
        base_report={'BKP':reduced.report,'result_assumption':self.residue.assumption,
                     'class_generation':'Minkowski prime witnesses',
                     'volume_convention':'Euclidean Gram volume; slope factor included'}
        if reduced.rank!=expected:
            return GenerationCheck(False,reduced.rank,expected,M,elements,values,None,None,None,
                                   reduced.precision,dict(base_report,state='rank_deficient'))
        low,high,p=gram_volume(lambda bits:M*oracle(bits),self.limits,reduced.precision)
        R=RealBallField(max(p,self.residue.precision))
        # The normal to the degree-zero hyperplane is (log N(P),...,1,...,1).
        # Its length is needed for the ambient Euclidean Gram determinant.
        slope=(R(len(self.context.places))+sum((R(P.norm()).log()**2 for P in self.factor_base.primes),R(0))).sqrt()
        target_lower=endpoints(R(self.residue.hR_lower)*slope)[0]
        target_upper=endpoints(R(self.residue.hR_upper)*slope)[1]
        index_low=QQ(low)/target_upper
        index_high=QQ(high)/target_lower
        integer_low=max(1,int(index_low.ceil()))
        integer_high=int(index_high.floor())
        if integer_high<integer_low:
            raise ArithmeticError("covolume intervals contain no positive integer index")
        complete=integer_low==integer_high==1
        state='complete' if complete else 'proper_sublattice' if integer_low>1 else 'undetermined'
        return GenerationCheck(complete,reduced.rank,expected,M,elements,values,(low,high),
            (target_lower,target_upper),(integer_low,integer_high),p,dict(base_report,state=state))


@dataclass(frozen=True)
class SUnitSystem:
    context: object
    factor_base: object
    generators: tuple
    valuations: object
    torsion_order: int
    torsion_generator: object
    exponent_matrix: object
    source_elements: tuple
    working_factor_base: object
    source_valuations: object
    report: dict
    generation_check: object=field(repr=False)

    @property
    def rank(self):return len(self.generators)

    def materialize_generators(self,max_bits=100000):
        return tuple(g.materialize(self.working_factor_base,max_bits=max_bits) for g in self.generators)

    def to_dict(self):
        """Exact, reviewable data; loading this is not a substitute for certification."""
        def ideal(P):return [[str(q) for q in a.vector()] for a in P.basis()]
        check=self.generation_check
        return {
            'format':'rigorous-cnt-sunit-v1',
            'defining_polynomial':[str(q) for q in self.context.K.defining_polynomial()],
            'requested_prime_ideals':[ideal(P) for P in self.factor_base.primes],
            'working_prime_ideals':[ideal(P) for P in self.working_factor_base.primes],
            'source_elements':[[{'coefficients':[str(q) for q in a.vector()],'exponent':str(e)}
                                for a,e in g.factors] for g in self.source_elements],
            'exponent_matrix':[[str(x) for x in row] for row in self.exponent_matrix.rows()],
            'source_valuations':[[str(x) for x in row] for row in self.source_valuations.rows()],
            'valuations':[[str(x) for x in row] for row in self.valuations.rows()],
            'torsion_order':self.torsion_order,
            'torsion_generator':[str(q) for q in self.torsion_generator.vector()],
            'completion_index_interval':list(check.index_interval),
            'sampled_lattice_volume':[str(x) for x in check.euclidean_volume],
            'target_lattice_volume':[str(x) for x in check.target_euclidean_volume],
            'report':self.report,
        }

    def verify(self):
        if self.rank!=len(self.factor_base)+len(self.context.places)-1:
            return False
        if self.torsion_generator**self.torsion_order!=1:return False
        for element,row in zip(self.source_elements,self.source_valuations.rows()):
            values=[ZZ(0)]*len(self.working_factor_base)
            for a,e in element.factors:
                v=self.working_factor_base.factor(self.context.K.ideal(a),nonnegative=False)
                if v is None:return False
                values=[x+e*y for x,y in zip(values,v)]
            if tuple(values)!=tuple(row):return False
        full_values=self.exponent_matrix*self.source_valuations
        for j,P in enumerate(self.working_factor_base.primes):
            column=full_values.column(j)
            if P in self.factor_base.index:
                if column!=self.valuations.column(self.factor_base.index[P]):return False
            elif any(column):return False
        expected=tuple(combine_compact(self.context,self.source_elements,row)
                       for row in self.exponent_matrix.rows())
        return expected==self.generators and self.generation_check.complete


def algorithm6(context,S=(),settings=None,residue=None,rng=None,limits=Limits(),progress=None):
    """Compute a certified fundamental S-unit system in compact representation.

    Uses Algorithms 3–5 for every collected generator. The only external group
    invariant used for completion is an independently computed GRH residue bound.
    The reference sampling parameters do not imply the paper's running-time bound.
    """
    settings=settings or SUnitSettings();rng=rng or RandomBits()
    engine=RelationEngine(context,settings,residue,limits)
    requested=FactorBase(context,S)
    minkowski=minkowski_prime_bound(context)
    if minkowski>settings.max_minkowski_bound:
        raise ResourceLimit("the class-generation proof requires too many small primes")
    actual_walk_bounds=[int(engine.sampler().walk.prime_bound)]
    actual_walk_bounds.extend(int(engine.sampler(engine.modulus/P).walk.prime_bound)
                              for P in engine.exceptional_primes)
    bound=max(settings.smooth_bound,max(actual_walk_bounds),minkowski)
    if bound>settings.max_factor_base_bound:
        raise ResourceLimit("working factor base exceeds the prime enumeration budget")
    all_primes=set(prime_ideals_up_to(context,bound))|set(requested.primes)|set(engine.exceptional_primes)
    working=FactorBase(context,all_primes)
    ordinary=FactorBase(context,[P for P in working.primes if P not in engine.exceptional_primes])
    exceptionals={}
    witnesses=[]
    for P in prime_ideals_up_to(context,minkowski):
        if P not in ordinary.index:
            relation=algorithm5(engine,P,ordinary,rng)
            exceptionals[P]=relation;witnesses.append(relation)
    class_cert=certify_class_generation(ordinary,witnesses,settings.max_minkowski_bound)
    post=SUnitPostprocessor(ordinary,engine.residue,class_cert,limits)
    observations=[]
    rank=len(ordinary)+len(context.places)-1
    check=post.check(()) if rank==0 else None
    for iteration in range(1,settings.max_relations+1):
        if check is not None and check.complete:break
        observation=algorithm4(engine,ordinary,rng)
        observations.append(observation)
        if len(observations)>=rank and (len(observations)-rank)%settings.check_every==0:
            check=post.check(observations)
        if progress is not None:
            progress({'relations':len(observations),'target_rank':rank,
                      'sample_attempts':engine.statistics['sampling_calls'],
                      'rank':None if check is None else check.rank,
                      'index':None if check is None else check.index_interval,
                      'complete':bool(check and check.complete)})
    else:
        if not (check and check.complete):
            raise ResourceLimit("Algorithm 6 exhausted its relation budget without a complete-generation certificate")
    if check is None or not check.complete:
        raise ResourceLimit("Algorithm 6 did not certify a full generating set")
    # Compute any remaining exceptional units after the ordinary T-unit lattice.
    for P in engine.exceptional_primes:
        if P not in exceptionals:
            exceptionals[P]=algorithm5(engine,P,ordinary,rng)
    source=[o.element for o in observations]
    source_values=[]
    for o in observations:
        source_values.append([o.valuations[ordinary.index[P]] if P in ordinary.index else 0
                              for P in working.primes])
    num_ordinary=len(source)
    for P in engine.exceptional_primes:
        relation=exceptionals[P]
        source.append(CompactElement.of(context,relation.element))
        source_values.append([1 if Q==P else relation.valuations[ordinary.index[Q]] if Q in ordinary.index else 0
                              for Q in working.primes])
    extended=matrix(ZZ,check.rank+len(engine.exceptional_primes),len(source))
    for i,row in enumerate(check.transform.rows()):
        for j,value in enumerate(row):extended[i,j]=value
    for i,P in enumerate(engine.exceptional_primes):extended[check.rank+i,num_ordinary+i]=1
    valuations=matrix(ZZ,len(source),len(working),[x for row in source_values for x in row])
    full_values=extended*valuations
    extras=[j for j,P in enumerate(working.primes) if P not in requested.index]
    if extras:
        extra_values=full_values.matrix_from_columns(extras)
        restriction=extra_values.transpose().right_kernel_matrix()
        final_matrix=restriction*extended
    else:
        final_matrix=extended
    final_full_values=final_matrix*valuations
    wanted_columns=[working.index[P] for P in requested.primes]
    final_values=final_full_values.matrix_from_columns(wanted_columns)
    generators=tuple(combine_compact(context,source,row) for row in final_matrix.rows())
    # The returned basis also gives an a posteriori check of the chosen radius
    # for the ordinary T-unit Gaussian (it does not certify unrelated time bounds).
    radius_verified=False
    if check.basis_elements:
        R=RealBallField(max(192,check.precision))
        maximum=QQ(0)
        for element,vals in zip(check.basis_elements,check.basis_valuations.rows()):
            row=vector(R,list(vals)+list(element.logarithms(R.precision())))
            maximum=max(maximum,endpoints(sum(x*x for x in row).sqrt())[1])
        radius_verified=maximum<=QQ(settings.generation_radius)
    else:radius_verified=True
    ray_generation=_certify_ray_generation(engine,check)
    mixing_statuses={o.report['mixing_guarantee'] for o in observations}
    mixing_statuses.update(r.sample_report['mixing_guarantee'] for r in exceptionals.values())
    mixing_status=('not_used' if not mixing_statuses else 'not_established' if 'not_established' in mixing_statuses
                   else 'verified_spectral_bound' if mixing_statuses=={'verified_spectral_bound'}
                   else 'conditional_on_supplied_bound')
    result=SUnitSystem(context,requested,generators,final_values,engine.residue.torsion_order,
        engine.residue.torsion_generator,final_matrix,tuple(source),working,valuations,{
            'algorithm':6,'complete_generation':'certified_under_GRH',
            'residue_method':'Belabas–Friedman prime-power sum','residue_cutoff':engine.residue.cutoff,
            'relations':len(observations),'exceptional_units':len(exceptionals),
            'sampling_calls':engine.statistics['sampling_calls'],'working_primes':len(working),
            'requested_primes':len(requested),'modulus_branch':engine.modulus_branch,
            'generation_radius_verified_after_completion':radius_verified,
            'ray_class_generation_verified_after_completion':ray_generation,
            'sampling_mixing_guarantee':mixing_status,
            'restriction':'saturated integer kernel of extra valuations',
            'paper_bit_complexity_certified':False,
        },check)
    if not result.verify():
        raise ArithmeticError("final compact S-unit system failed exact validation")
    return result
