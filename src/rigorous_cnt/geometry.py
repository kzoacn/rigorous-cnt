"""Proposition 8.8: rational-grid sampling transferred to an exact lattice."""

from dataclasses import dataclass
from sage.all import AA, QQ, ZZ, matrix, vector
from .numerics import nearest
from .runtime import Limits, RandomBits, ResourceLimit


class ProductRegion:
    """A centered product of real intervals and complex disks in real coordinates."""

    def __init__(self, real_radii, complex_radii):
        self.real_radii = tuple(AA(x) for x in real_radii)
        self.complex_radii = tuple(AA(x) for x in complex_radii)
        self.n = len(self.real_radii)+2*len(self.complex_radii)
        if not self.n or any(x <= 0 for x in self.real_radii+self.complex_radii):
            raise ValueError("region radii must be positive")

    def contains(self, v, dilation=1):
        d = AA(dilation)
        if len(v) != self.n:
            raise ValueError("region coordinate dimension mismatch")
        k = 0
        for r in self.real_radii:
            if abs(v[k]) > d*r:
                return False
            k += 1
        for r in self.complex_radii:
            if v[k]**2+v[k+1]**2 > (d*r)**2:
                return False
            k += 2
        return True

    def sample_grid(self, denominator, rng, dilation=1, limits=Limits()):
        N, d = ZZ(denominator), AA(dilation)
        out = []
        for r in self.real_radii:
            bound = ZZ((N*d*r).floor())
            out.append(QQ(rng.integer(-bound,bound))/N)
        for r in self.complex_radii:
            bound = ZZ((N*d*r).floor())
            for _ in range(limits.attempts):
                x,y = ZZ(rng.integer(-bound,bound)),ZZ(rng.integer(-bound,bound))
                if x*x+y*y <= (N*d*r)**2:
                    out.extend((QQ(x)/N,QQ(y)/N))
                    break
            else:
                raise ResourceLimit("complex disk rejection budget exhausted")
        return vector(QQ,out)


@dataclass(frozen=True)
class GridCertificate:
    denominator: object
    epsilon: object
    D: object
    U: object
    approximation: object
    approximate_shift: object
    inverse_norm_upper: object
    error_norm_upper: object


class UniformLatticeSampler:
    """Prepare Proposition 8.8 for a column basis B and target S+t.

    All six geometric hypotheses are checked using exact algebraic arithmetic.
    Frobenius bounds are used to dominate operator norms. The grid denominator
    is selected from those verified quantities, rather than a transcription of
    the paper's asymptotic denominator formula.
    """

    def __init__(self,B,region,shift,epsilon=None,limits=Limits()):
        B = matrix(AA,B)
        n = B.ncols()
        if B.nrows() != n or n != region.n or not B.det():
            raise ValueError("a nonsingular square lattice basis is required")
        self.B,self.region,self.shift,self.limits = B,region,vector(AA,shift),limits
        if len(self.shift) != n:
            raise ValueError("shift dimension mismatch")
        epsilon = QQ(1)/(6*n) if epsilon is None else QQ(epsilon)
        if not 0 < epsilon < QQ(1)/5:
            raise ValueError("Proposition 8.8 requires 0<epsilon<1/5")
        D = 2*sum(sum(x*x for x in B.column(i)).sqrt() for i in range(n))
        inner_radius = min(region.real_radii+region.complex_radii)
        if D/epsilon > inner_radius:
            raise ValueError("basis is too long for the region in Proposition 8.8")
        outer = sum(r*r for r in region.real_radii+region.complex_radii).sqrt()
        shift_norm = sum(x*x for x in self.shift).sqrt()
        U = 2*max(outer,shift_norm,AA(1))
        inv = B.inverse()
        inv_norm = sum(x*x for x in inv.list()).sqrt()
        # Rounding error ||C-B||_2 <= ||C-B||_F <= n/(2N).
        required = max(AA(n)/D, AA(n)*inv_norm*U/D, AA(n)**(QQ(3)/2)/D, AA(1))
        N = ZZ(1)
        while N < required:
            N *= 2
        for _ in range(limits.lattice_steps):
            C = matrix(QQ,n,n,lambda i,j: QQ(nearest(N*B[i,j]))/N)
            t = vector(QQ,[QQ(nearest(N*x))/N for x in self.shift])
            error = sum((C[i,j]-B[i,j])**2 for i in range(n) for j in range(n)).sqrt()
            lengths = sum(AA(sum(x*x for x in C.column(i))).sqrt() for i in range(n))
            if (C.det() and error*inv_norm <= D/U and lengths <= D
                    and sum(x*x for x in t) < (2*U)**2):
                break
            N *= 2
        else:
            raise ResourceLimit("rational lattice approximation budget exhausted")
        self.C,self.C_inverse = C,C.inverse()
        self.w0 = vector(ZZ,[nearest(x) for x in self.C_inverse*t])
        self.epsilon = epsilon
        self.certificate = GridCertificate(N,epsilon,D,U,C,t,inv_norm,error)

    def sample(self,rng=None,accept=None):
        rng = rng or RandomBits()
        e,N = self.epsilon,self.certificate.denominator
        for attempt in range(1,self.limits.attempts+1):
            u = self.region.sample_grid(N,rng,1+4*e,self.limits)
            v = vector(ZZ,[nearest(x) for x in self.C_inverse*u])
            if not self.region.contains(self.C*v,1+3*e):
                continue
            result = v+self.w0
            if self.region.contains(self.B*result-self.shift) and (accept is None or accept(result)):
                return result,attempt
        raise ResourceLimit("uniform lattice rejection budget exhausted")
