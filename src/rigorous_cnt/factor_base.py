"""Prime factor bases and exact principal-ideal relation checks."""

from sage.all import QQ,ZZ,prime_range


def prime_ideals_up_to(context,bound):
    bound=ZZ(bound)
    return tuple(P for p in prime_range(bound+1) for P in context.K.primes_above(p)
                 if P.norm()<=bound)


class FactorBase:
    def __init__(self,context,primes=()):
        self.context=context
        unique=set()
        for P in primes:
            P=context.K.ideal(P)
            if not P or not P.is_integral() or not P.is_prime():
                raise ValueError("a factor base consists of nonzero prime ideals")
            unique.add(P)
        self.primes=tuple(sorted(unique,key=lambda P:(P.norm(),tuple(tuple(a.vector()) for a in P.basis()))))
        self.index={P:i for i,P in enumerate(self.primes)}
        self.rational_primes=tuple(sorted({ZZ(P.smallest_integer()) for P in self.primes}))

    def __len__(self):
        return len(self.primes)

    def factor(self,ideal,nonnegative=True):
        ideal=self.context.K.ideal(ideal)
        if not ideal:
            raise ValueError("cannot factor the zero ideal over a factor base")
        if nonnegative and not ideal.is_integral():
            return None
        norm=QQ(ideal.norm())
        for part in (abs(norm.numerator()),norm.denominator()):
            remainder=ZZ(part)
            for p in self.rational_primes:
                while remainder%p==0:
                    remainder//=p
            if remainder!=1:
                return None
        values=tuple(ZZ(ideal.valuation(P)) for P in self.primes)
        if nonnegative and any(e<0 for e in values):
            return None
        reconstructed=self.ideal(values)
        return values if reconstructed==ideal else None

    def ideal(self,exponents):
        if len(exponents)!=len(self):
            raise ValueError("factor-base dimension mismatch")
        out=self.context.K.ideal(1)
        for P,e in zip(self.primes,exponents):
            if ZZ(e)!=e:
                raise ValueError("ideal exponents must be integers")
            out*=P**ZZ(e)
        return out
