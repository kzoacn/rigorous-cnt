"""Exact algebraic Gram matrices and checked short-basis construction.

LLL is used first. If necessary, block HKZ tours follow Algorithm 7's schedule.
The enumerator is a reference SVP implementation, not a claim to the Kannan
bit-complexity used in the paper. All returned geometric bounds are checked.
"""

from dataclasses import dataclass
from sage.all import AA, QQ, ZZ, RealBallField, matrix, vector, identity_matrix, block_diagonal_matrix
from .numerics import nearest,endpoints
from .runtime import Limits, ResourceLimit, UnresolvedBound


def proved_integer_lll(A):
    """Use exact NTL if fpLLL cannot finish at the eta=1/2 boundary."""
    from fpylll.util import ReductionError
    try:
        reduced,U=A.LLL(delta=QQ(3)/4,eta=QQ(1)/2,algorithm='fpLLL:proved',transformation=True)
        return reduced,U,'fpLLL:proved'
    except ReductionError:
        reduced,U=A.LLL(delta=QQ(3)/4,algorithm='NTL:LLL',transformation=True)
        return reduced,U,'NTL:LLL exact fallback'


def gram_schmidt(G):
    G = matrix(AA, G)
    n = G.nrows()
    if G.ncols() != n or G != G.transpose():
        raise ValueError("a symmetric Gram matrix is required")
    mu = identity_matrix(AA, n)
    norms = []
    for i in range(n):
        for j in range(i):
            mu[i,j] = (G[i,j]-sum(mu[i,k]*mu[j,k]*norms[k] for k in range(j)))/norms[j]
        norm = G[i,i]-sum(mu[i,k]**2*norms[k] for k in range(i))
        if norm <= 0:
            raise ValueError("Gram matrix is not positive definite")
        norms.append(norm)
    return mu, norms


def lll_transform(G, limits=Limits()):
    G = matrix(AA, G)
    n = G.nrows()
    U = identity_matrix(ZZ, n)
    k, steps = 1, 0
    while k < n:
        steps += 1
        if steps > limits.lattice_steps:
            raise ResourceLimit("exact LLL step budget exhausted")
        for j in reversed(range(k)):
            mu, norms = gram_schmidt(U*G*U.transpose())
            q = nearest(mu[k,j])
            if q:
                U.set_row(k, U.row(k)-q*U.row(j))
        mu, norms = gram_schmidt(U*G*U.transpose())
        if norms[k] >= (QQ(3)/4-mu[k,k-1]**2)*norms[k-1]:
            k += 1
        else:
            U.swap_rows(k, k-1)
            k = max(1, k-1)
    return U


def _nearby_integers(center, lo, hi):
    mid = nearest(center)
    if lo <= mid <= hi:
        yield mid
    distance = 1
    while mid-distance >= lo or mid+distance <= hi:
        if mid-distance >= lo:
            yield mid-distance
        if mid+distance <= hi:
            yield mid+distance
        distance += 1


def shortest_vector_coefficients(G, limits=Limits()):
    """Enumerate all relevant GSO branches with exact algebraic comparisons."""
    G = matrix(AA, G)
    n = G.nrows()
    if n == 0:
        raise ValueError("zero-dimensional SVP")
    mu, norms = gram_schmidt(G)
    best = [AA(G[0,0]), tuple([ZZ(1)]+[ZZ(0)]*(n-1))]
    coeff = [ZZ(0)]*n
    nodes = [0]

    def visit(k, partial):
        nodes[0] += 1
        if nodes[0] > limits.enumeration_nodes:
            raise ResourceLimit("exact SVP enumeration node budget exhausted")
        if k < 0:
            if any(coeff) and partial < best[0]:
                best[:] = [partial, tuple(coeff)]
            return
        slack = best[0]-partial
        if slack < 0:
            return
        center = -sum(mu[j,k]*coeff[j] for j in range(k+1,n))
        radius = (slack/norms[k]).sqrt()
        lo, hi = ZZ((center-radius).ceil()), ZZ((center+radius).floor())
        for z in _nearby_integers(center, lo, hi):
            coeff[k] = z
            new = partial+(z-center)**2*norms[k]
            if new <= best[0]:
                visit(k-1, new)
        coeff[k] = 0

    visit(n-1, AA(0))
    return best[1], best[0]


def hkz_transform(G, limits=Limits()):
    G = matrix(AA, G)
    n = G.nrows()
    if n <= 1:
        return identity_matrix(ZZ,n)
    coeff, length = shortest_vector_coefficients(G, limits)
    if length < G[0,0]:
        D, L, R = matrix(ZZ,n,1,coeff).smith_form()
        if abs(D[0,0]) != 1:
            raise ArithmeticError("SVP returned a nonprimitive vector")
        U = L.inverse().transpose().change_ring(ZZ)
    else:
        U = identity_matrix(ZZ,n)
    H = U*G*U.transpose()
    tail = matrix(AA,n-1,n-1,lambda i,j: H[i+1,j+1]-H[i+1,0]*H[0,j+1]/H[0,0])
    V = hkz_transform(tail, limits)
    U = block_diagonal_matrix(identity_matrix(ZZ,1),V)*U
    return _size_reduce_transform(G, U)


def _size_reduce_transform(G, U):
    U = matrix(ZZ,U)
    for i in range(U.nrows()):
        for j in reversed(range(i)):
            mu, _ = gram_schmidt(U*G*U.transpose())
            q = nearest(mu[i,j])
            if q:
                U.set_row(i,U.row(i)-q*U.row(j))
    return U


def _projected_block(G, start, size):
    mu, norms = gram_schmidt(G)
    return matrix(AA,size,size,lambda i,j: sum(
        mu[start+i,k]*mu[start+j,k]*norms[k]
        for k in range(start,min(start+i,start+j)+1)))


@dataclass(frozen=True)
class ShortBasis:
    elements: tuple
    transform: object
    squared_lengths: tuple
    length_bound: object
    hkz_tours: int
    complexity_certified: bool = False
    method: str = 'exact_algebraic_LLL_and_block_HKZ'


def _rational_lll_candidate(context,basis,B,scales,bound,limits):
    """Fast candidate generation; only a checked exact lattice basis is accepted.

    LLL is applied to an integer approximation. Its transformation is applied
    to the ORIGINAL algebraic elements, and the required geometric bound is
    certified independently. No floating-point reduction property is assumed.
    """
    precision=64
    while precision<=limits.precision:
        R=RealBallField(precision)
        numeric=matrix(R,B)
        for i in range(context.r1,context.n):
            numeric.set_row(i,R(2).sqrt()*numeric.row(i))
        max_entry=max(max(abs(lo),abs(hi)) for lo,hi in map(endpoints,numeric.list()))
        exponent=int(max_entry.numerator().nbits()-max_entry.denominator().nbits())
        power=QQ(2)**(precision//2-exponent)
        rounded=matrix(ZZ,context.n,context.n,
                       lambda i,j:nearest(power*numeric[j,i].mid().exact_rational()))
        if rounded.rank()<context.n:
            precision*=2;continue
        _,U,_=proved_integer_lll(rounded)
        elements=tuple(sum((U[i,j]*basis[j] for j in range(context.n)),context.K(0))
                       for i in range(context.n))
        coordinates=[context.embed(a,scales) for a in elements]
        upper=[]
        for row in coordinates:
            squared=sum(R(x)**2 for x in row[:context.r1])+2*sum(R(x)**2 for x in row[context.r1:])
            upper.append(endpoints(squared)[1])
        lower_bound=endpoints(R(bound)**2)[0]
        if all(x<=lower_bound for x in upper):
            if abs(U.det())!=1:
                raise ArithmeticError("approximate LLL transformation was not unimodular")
            exact_lengths=tuple(context.weighted_norm2(row) for row in coordinates)
            return ShortBasis(elements,U,exact_lengths,AA(bound),0,False,
                              'integer_LLL_candidate_with_certified_geometric_bound')
        # A failed candidate never escapes. More precision or the exact
        # block-HKZ path below is required.
        precision*=2
        if precision>2048:
            break
    return None


def short_ideal_basis(context, ideal, scales, bound, block_size=2, limits=Limits()):
    n = context.n
    if int(block_size) != block_size or not 2 <= block_size <= n:
        raise ValueError("block_size must be an integer between 2 and the degree")
    block_size = int(block_size)
    basis = tuple(ideal.basis())
    B = context.embedding_matrix(basis,scales)
    candidate=_rational_lll_candidate(context,basis,B,scales,bound,limits)
    if candidate is not None:
        return candidate
    W = matrix.diagonal(AA,[1]*context.r1+[2]*(2*context.r2))
    G = B.transpose()*W*B
    U = lll_transform(G,limits)
    tours = 0
    while True:
        H = U*G*U.transpose()
        lengths = tuple(H[i,i] for i in range(n))
        if all(x <= bound**2 for x in lengths):
            break
        if tours >= limits.lattice_steps:
            raise ResourceLimit("block-HKZ tour budget exhausted")
        before = matrix(ZZ,U)
        for start in range(n-block_size+1):
            H = U*G*U.transpose()
            V = hkz_transform(_projected_block(H,start,block_size),limits)
            block = V*U.matrix_from_rows(range(start,start+block_size))
            for j in range(block_size):
                U.set_row(start+j,block.row(j))
            U = _size_reduce_transform(G,U)
        tours += 1
        if U == before:
            raise UnresolvedBound("block-HKZ stabilized without the required short-basis bound")
    if abs(U.det()) != 1:
        raise ArithmeticError("basis transformation is not unimodular")
    new_basis = tuple(sum((U[i,j]*basis[j] for j in range(n)),context.K(0)) for i in range(n))
    return ShortBasis(new_basis,U,lengths,AA(bound),tours)
