"""Independent GRH-conditional Dedekind residue intervals.

Belabas–Friedman, arXiv:1305.0035, Theorem 1. Only prime splitting and
roots of unity are used; neither class groups nor fundamental units are computed.
"""

from dataclasses import dataclass,field
from sage.all import QQ,ZZ,RealBallField,prime_range
from .numerics import rational,endpoints
from .runtime import ResourceLimit


@dataclass(frozen=True)
class ResidueCertificate:
    rho_lower: object
    rho_upper: object
    hR_lower: object
    hR_upper: object
    log_error: object
    cutoff: int
    precision: int
    torsion_order: int
    torsion_generator: object
    _context: object=field(repr=False,compare=False)
    assumption: str="GRH for zeta_K and zeta_Q (Belabas–Friedman Theorem 1)"

    def validate_for(self,context):
        if context is not self._context:
            raise ValueError("residue certificate belongs to another number-field context")
        if not 0<self.rho_lower<=self.rho_upper or not 0<self.hR_lower<=self.hR_upper:
            raise ValueError("invalid positive residue interval")


def _error_bound(R,n,D,X):
    logD=R(D).log(); rootX=R(X).sqrt()
    return (R('2.324')*logD/(rootX*R(3*X).log()) *
            ((1+R('3.88')/R(QQ(X)/9).log())*(1+2/logD.sqrt())**2
             +R('4.26')*(n-1)/(rootX*logD)))


def _weighted_prime_sum(R,cutoff,splitting):
    X=QQ(cutoff)
    sqrtX=R(X).sqrt(); logX=R(X).log()
    total=R(0)
    for p,degrees in splitting:
        if p>=X:
            continue
        logp=R(p).log()
        power,j=ZZ(p),1
        while power<X:
            coefficient=sum(f for f in degrees if j%f==0)-1
            if coefficient:
                total+=coefficient*(sqrtX*logX/(j*R(power))-logp/R(power).sqrt())
            power*=p; j+=1
    return total


def dedekind_residue(context,log_error=QQ(1)/16,max_cutoff=262144,precision=192):
    """Enclose rho_K and h_K R_K using explicit prime-power sums.

    The error bound is analytic (conditional on GRH), not an estimate inferred
    from the stability of successive numerical values.
    """
    tolerance=rational(log_error)
    if not 0<tolerance<QQ(1)/3 or precision<64:
        raise ValueError("require 0<log_error<1/3 and precision>=64")
    R=RealBallField(precision)
    X=128
    while X<=max_cutoff:
        error=_error_bound(R,context.n,context.discriminant,X)
        error_upper=endpoints(error)[1]
        if error_upper<=tolerance:
            break
        X*=2
    else:
        raise ResourceLimit("Dedekind residue cutoff exceeds the requested budget")
    splitting=[]
    for p in prime_range(X):
        degrees=tuple(int(P.residue_class_degree()) for P in context.K.primes_above(p))
        splitting.append((ZZ(p),degrees))
    B=_weighted_prime_sum(R,X,splitting)
    B9=_weighted_prime_sum(R,QQ(X)/9,splitting)
    estimate=3*(B-B9)/(2*R(X).sqrt()*R(3*X).log())
    lo,hi=endpoints(estimate)
    rho_lower=endpoints(R(lo-error_upper).exp())[0]
    rho_upper=endpoints(R(hi+error_upper).exp())[1]
    nf=context.K.pari_nf()
    roots=nf.nfrootsof1()
    torsion_order=int(roots[0])
    value=roots[1]
    if value.type()=='t_COL':
        value=nf.nfbasistoalg(value)
    torsion_generator=context.K(value)
    if torsion_generator**torsion_order!=1:
        raise ArithmeticError("root-of-unity data failed exact validation")
    factor=R(torsion_order)*R(context.discriminant).sqrt()/(R(2)**context.r1*(2*R.pi())**context.r2)
    lower=endpoints(R(rho_lower)*factor)[0]
    upper=endpoints(R(rho_upper)*factor)[1]
    result=ResidueCertificate(rho_lower,rho_upper,lower,upper,error_upper,X,precision,
                             torsion_order,torsion_generator,context)
    result.validate_for(context)
    return result
