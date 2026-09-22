"""Exact power-product representations of number-field elements."""

from dataclasses import dataclass
from sage.all import QQ,ZZ,RealBallField,vector,AA,matrix,PolynomialRing,Zmod,next_prime
from .runtime import ResourceLimit,Limits
from .numerics import endpoints


@dataclass(frozen=True,init=False)
class CompactElement:
    context: object
    factors: tuple

    def __init__(self,context,factors=()):
        combined={}
        for a,e in factors:
            a=context.K(a); exponent=ZZ(e)
            if exponent!=e:
                raise ValueError("power-product exponents must be integers")
            if not a:
                raise ValueError("a power product cannot contain zero")
            if exponent:
                combined[a]=combined.get(a,ZZ(0))+exponent
        object.__setattr__(self,'context',context)
        object.__setattr__(self,'factors',tuple((a,e) for a,e in combined.items() if e and a!=1))

    @classmethod
    def of(cls,context,a):
        return cls(context,[(a,1)])

    def __mul__(self,other):
        if not isinstance(other,CompactElement) or other.context is not self.context:
            raise ValueError("power products must share a context")
        return CompactElement(self.context,self.factors+other.factors)

    def __pow__(self,e):
        e0=ZZ(e)
        if e0!=e:
            raise ValueError("integer exponent required")
        return CompactElement(self.context,[(a,e0*k) for a,k in self.factors])

    def valuation(self,P):
        return sum((e*a.valuation(P) for a,e in self.factors),ZZ(0))

    def logarithms(self,precision=128,limits=Limits()):
        """Enclose weighted logarithms with width <= 2**(-precision).

        Working precision includes the size of the exponents; a tiny resulting
        unit may be represented by power products with enormous cancellation.
        """
        precision=int(precision)
        if precision<2:raise ValueError("at least two precision bits are required")
        amplification=sum((abs(e) for a,e in self.factors),ZZ(1))
        work=precision+int(amplification.nbits())+24
        tolerance=QQ(2)**(-precision)
        while work<=limits.precision:
            R=RealBallField(work)
            out=vector(R,[0]*len(self.context.places))
            valid=True
            for a,e in self.factors:
                for i,(emb,w) in enumerate(zip(self.context.places,self.context.weights)):
                    z=emb(a)
                    squared=AA(z.real())**2+AA(z.imag())**2
                    value=R(squared)
                    if not value>0:
                        valid=False;break
                    out[i]+=e*(QQ(w)/2)*value.log()
                if not valid:break
            if valid and all(hi-lo<=tolerance for lo,hi in map(endpoints,out)):
                return out
            work*=2
        raise ResourceLimit("compact logarithm cancellation exceeds the precision budget")

    def materialize(self,factor_base=None,max_bits=100000,precision=128,limits=Limits()):
        """Recover a small exact result without expanding enormous powers.

        Exact ideal valuations give a denominator. Certified archimedean bounds
        bound the integral-basis coefficients; modular power products then
        determine these integers uniquely. This can be much cheaper than expand.
        """
        K=self.context.K
        valuations={}
        for a,e in self.factors:
            if factor_base is None:
                factors=K.ideal(a).factor()
            else:
                if factor_base.context is not self.context:
                    raise ValueError("factor-base context mismatch")
                values=factor_base.factor(K.ideal(a),nonnegative=False)
                if values is None:
                    raise ValueError("a power-product factor lies outside the supplied factor base")
                factors=zip(factor_base.primes,values)
            for P,v in factors:valuations[P]=valuations.get(P,ZZ(0))+e*v
        denominator_powers={}
        for P,v in valuations.items():
            if v<0:
                p=ZZ(P.smallest_integer());power=ZZ((-QQ(v)/P.ramification_index()).ceil())
                denominator_powers[p]=max(denominator_powers.get(p,ZZ(0)),power)
        if sum(e*p.nbits() for p,e in denominator_powers.items())>max_bits:
            raise ResourceLimit("materialized denominator exceeds the requested bit budget")
        denominator=ZZ(1)
        for p,e in denominator_powers.items():denominator*=p**e
        logs=self.logarithms(precision,limits)
        R=RealBallField(max(precision,logs.base_ring().precision()))
        inv=self.context.embedding_matrix(self.context.integral_basis).inverse()
        inverse_bound=max(endpoints(R(sum(abs(x) for x in row)))[1] for row in inv.rows())
        max_log=max(endpoints(x)[1]/w for x,w in zip(logs,self.context.weights))
        height_log=R(max_log)+R(denominator).log()+R(inverse_bound).log()
        if endpoints(height_log/R(2).log())[1]>max_bits:
            raise ResourceLimit("materialized coefficients exceed the requested bit budget")
        coefficient_bound=max(ZZ(1),ZZ(endpoints(height_log.exp())[1].ceil()))
        basis=matrix(QQ,[a.vector() for a in self.context.integral_basis])
        inverse_basis=basis.inverse()
        polynomial=K.defining_polynomial().monic()
        # Replace negative powers by exact inverse factors, so modular arithmetic
        # only needs multiplication over Z/(ell**k), not polynomial inversion.
        factors=[(a if e>=0 else 1/a,abs(e)) for a,e in self.factors]
        denominators=[q.denominator() for a,e in factors for q in a.vector()]
        denominators += [q.denominator() for q in inverse_basis.list()+list(polynomial)]
        ell=ZZ(2)
        while any(d%ell==0 for d in denominators):ell=next_prime(ell)
        per_digit=int(ell.nbits())-1
        exponent=max(1,(int(coefficient_bound.nbits())+1+per_digit-1)//per_digit)
        modulus=ell**exponent
        ring=Zmod(modulus)
        def modq(q):
            q=QQ(q)
            return ring(q.numerator())/ring(q.denominator())
        PR=PolynomialRing(ring,'z')
        quotient=PR.quotient(PR([modq(q) for q in polynomial]),'theta')
        residue=quotient(denominator)
        for a,e in factors:
            residue*=quotient(PR([modq(q) for q in a.vector()]))**e
        coefficients=list(residue.lift())
        coefficients += [ring(0)]*(self.context.n-len(coefficients))
        change=matrix(ring,self.context.n,self.context.n,[modq(q) for q in inverse_basis.list()])
        coordinates=vector(ring,coefficients)*change
        integers=[]
        for c in coordinates:
            z=ZZ(c)
            if z>modulus//2:z-=modulus
            if abs(z)>coefficient_bound:
                raise ArithmeticError("modular reconstruction contradicts the certified height bound")
            integers.append(z)
        return sum((c*a for c,a in zip(integers,self.context.integral_basis)),K(0))/denominator

    def expand(self,max_bits=1000000):
        estimate=0
        for a,e in self.factors:
            size=max((abs(q.numerator()).nbits()+q.denominator().nbits() for q in a.vector()),default=1)
            estimate+=abs(e)*max(1,size)*self.context.n
        if estimate>max_bits:
            raise ResourceLimit("power-product expansion exceeds the requested coefficient budget")
        out=self.context.K(1)
        for a,e in self.factors:
            out*=a**e
        return out

    def __repr__(self):
        if not self.factors:
            return 'CompactElement(1)'
        return 'CompactElement('+ ' * '.join(f'({a})^{e}' for a,e in self.factors)+')'


def combine_compact(context,elements,exponents):
    if len(elements)!=len(exponents):
        raise ValueError("power-product dimension mismatch")
    factors=[]
    for element,e in zip(elements,exponents):
        element=element if isinstance(element,CompactElement) else CompactElement.of(context,element)
        if element.context is not context:
            raise ValueError("power-product context mismatch")
        factors.extend((a,ZZ(e)*k) for a,k in element.factors)
    return CompactElement(context,factors)
