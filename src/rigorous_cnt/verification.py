"""Recheck exported S-unit data without trusting its stored success flags."""

from sage.all import QQ,ZZ,RealBallField,PolynomialRing,NumberField,matrix,vector
from .arithmetic import NumberFieldContext
from .factor_base import FactorBase
from .compact import CompactElement,combine_compact
from .relations import SUnitObservation
from .residue import dedekind_residue
from .sunits import SUnitPostprocessor,certify_class_generation
from .postprocess import bkp_basis,kessler_lower_bound,gram_volume
from .runtime import Limits


def verify_sunit_record(data,limits=Limits()):
    """Recompute support, full source lattice, restriction, and final index.

    Requires the working factor base exported by algorithm6, including its
    Minkowski prime cover. Stored report/volume/index claims are ignored.
    """
    if data.get('format')!='rigorous-cnt-sunit-v1':
        raise ValueError("unsupported S-unit record format")
    polynomial_ring=PolynomialRing(QQ,'t')
    polynomial=polynomial_ring([QQ(x) for x in data['defining_polynomial']])
    K=NumberField(polynomial,'a');context=NumberFieldContext(K)
    def ideal(rows):
        return K.ideal([K([QQ(x) for x in row]) for row in rows])
    working=FactorBase(context,[ideal(rows) for rows in data['working_prime_ideals']])
    requested=FactorBase(context,[ideal(rows) for rows in data['requested_prime_ideals']])
    if any(P not in working.index for P in requested.primes):
        raise ValueError("requested primes are missing from the working base")
    source=tuple(CompactElement(context,[(K([QQ(x) for x in term['coefficients']]),ZZ(term['exponent']))
                                         for term in factors]) for factors in data['source_elements'])
    stored_values=data['source_valuations']
    if len(stored_values)!=len(source) or any(len(row)!=len(working) for row in stored_values):
        raise ValueError("source valuation dimensions are inconsistent")
    values=matrix(ZZ,len(source),len(working),[ZZ(x) for row in stored_values for x in row])
    observations=[SUnitObservation(g,tuple(row),working) for g,row in zip(source,values.rows())]
    residue=dedekind_residue(context)
    post=SUnitPostprocessor(working,residue,certify_class_generation(working),limits)
    source_check=post.check(observations)
    if not source_check.complete:
        return {'verified':False,'reason':'source generators do not span the full working S-unit group',
                'source_index_interval':source_check.index_interval}
    expected_rank=len(requested)+len(context.places)-1
    stored_matrix=data['exponent_matrix']
    if len(stored_matrix)!=expected_rank or any(len(row)!=len(source) for row in stored_matrix):
        raise ValueError("final exponent matrix dimensions are inconsistent")
    M=matrix(ZZ,expected_rank,len(source),[ZZ(x) for row in stored_matrix for x in row])
    final_values=M*values
    extras=[j for j,P in enumerate(working.primes) if P not in requested.index]
    if any(final_values[i,j] for i in range(expected_rank) for j in extras):
        return {'verified':False,'reason':'a final generator has support outside requested S'}
    wanted=[working.index[P] for P in requested.primes]
    claimed_values=data['valuations']
    if (len(claimed_values)!=expected_rank or any(len(row)!=len(requested) for row in claimed_values)
        or final_values.matrix_from_columns(wanted)!=matrix(ZZ,expected_rank,len(requested),
            [ZZ(x) for row in claimed_values for x in row])):
        return {'verified':False,'reason':'final valuations do not match the exponent matrix'}
    N=source_check.transform
    if extras:
        extra=(N*values).matrix_from_columns(extras)
        N=extra.transpose().right_kernel_matrix()*N
    if N.nrows()!=expected_rank:
        raise ArithmeticError("independent integer restriction gave an unexpected rank")
    def log_oracle(transform):
        elements=tuple(combine_compact(context,source,row) for row in transform.rows())
        finite=transform*values
        def oracle(p):
            R=RealBallField(p)
            rows=[vector(R,[-x for x in vals]+list(g.logarithms(p,limits)))
                  for vals,g in zip(finite.rows(),elements)]
            return matrix(R,len(rows),len(working)+len(context.places),[x for row in rows for x in row])
        return oracle
    candidate=log_oracle(M)
    rank_check=bkp_basis(candidate,kessler_lower_bound(context,limits),max(1,expected_rank),len(working),limits)
    if rank_check.rank!=expected_rank:
        return {'verified':False,'reason':'final generators are rank deficient'}
    low,high,p=gram_volume(candidate,limits)
    target_low,target_high,_=gram_volume(log_oracle(N),limits)
    interval=(max(1,int((low/target_high).ceil())),int((high/target_low).floor()))
    torsion=K([QQ(x) for x in data['torsion_generator']])
    order=ZZ(data['torsion_order'])
    torsion_ok=(order==residue.torsion_order and torsion**order==1
                and all(torsion**(order//q)!=1 for q,e in order.factor()))
    return {'verified':bool(interval==(1,1) and torsion_ok),'free_rank':expected_rank,
            'target_index_interval':interval,'torsion_verified':bool(torsion_ok),
            'assumption':residue.assumption,'stored_flags_trusted':False}
