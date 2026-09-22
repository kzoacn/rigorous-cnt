"""Algorithms 1 and 2, with exact outputs and explicit guarantee metadata."""

from dataclasses import dataclass
from sage.all import AA, QQ, ZZ, matrix, vector
from .arithmetic import NumberFieldContext,RayConditions
from .numerics import rational,nearest,norm_one_exp_approx
from .runtime import Limits,RandomBits
from .lattice import short_ideal_basis
from .geometry import ProductRegion,UniformLatticeSampler
from .gaussian import discrete_gaussian
from .parameters import WalkParameters,distortion_grid
from .primes import PrimeIdealSampler


def paper_radius(context,ideal,scales,modulus,block_size,omega):
    n,b = context.n,int(block_size)
    w = rational(omega)
    if b != block_size or not 2 <= b <= n or w < 1:
        raise ValueError("2<=block_size<=n and rational omega>=1 required")
    norm = abs(QQ(ideal.norm()))*QQ(modulus.norm())*context.scale_norm(scales)
    return (48*w*AA(b)**(QQ(2*n)/b)*AA(n)**(QQ(7)/2)
            *AA(context.discriminant)**(QQ(3)/(2*n))*AA(norm)**(QQ(1)/n))


@dataclass(frozen=True)
class BoxSample:
    element: object
    scale: tuple
    radius: object
    coefficients: tuple
    report: dict

    def point(self,context):
        """The Algorithm 1 output x*alpha, retained symbolically until requested."""
        return context.embed(self.element,self.scale)


class Algorithm1:
    """Prepare once, then sample uniformly from x*((ideal+gamma) intersect ray) intersect r B_inf."""

    def __init__(self,context,ideal,gamma=0,ray=None,x=None,block_size=2,omega=1,limits=Limits()):
        self.context,self.limits = context,limits
        self.ideal = context.K.ideal(ideal)
        if not self.ideal:
            raise ValueError("the input ideal must be nonzero")
        self.gamma = context.K(gamma)
        self.ray = ray or RayConditions(context)
        if self.ray.context is not context:
            raise ValueError("ray and algorithm must share a number-field context")
        self.scales = context.scales(x)
        self.radius = paper_radius(context,self.ideal,self.scales,self.ray.modulus,block_size,omega)
        gamma_m = self.ray.affine_shift(self.ideal,self.gamma)
        lattice_ideal = self.ideal*self.ray.modulus
        self.short_basis = short_ideal_basis(context,lattice_ideal,self.scales,
            self.radius/(24*context.n**2),block_size,limits)
        basis = self.short_basis.elements
        coords = context.coordinates_in_basis(gamma_m,basis)
        self.reduced_shift = gamma_m-sum((nearest(q)*v for q,v in zip(coords,basis)),context.K(0))
        B = context.embedding_matrix(basis,self.scales)
        real_radii = [self.radius]*context.r1
        region_center = vector(AA,[0]*context.n)
        for i in self.ray.real_places:
            # Sign of x_i*alpha_i includes the real scale, which may be negative.
            sign = 1 if self.scales[i][0]*context.places[i](self.ray.tau).real() > 0 else -1
            real_radii[i] /= 2
            region_center[i] = sign*self.radius/2
        region = ProductRegion(real_radii,[self.radius]*context.r2)
        shift = region_center-context.embed(self.reduced_shift,self.scales)
        self.grid = UniformLatticeSampler(B,region,shift,limits=limits)

    def _element(self,coordinates):
        return self.reduced_shift+sum((z*a for z,a in zip(coordinates,self.short_basis.elements)),self.context.K(0))

    def verify(self,element):
        element = self.context.K(element)
        if element-self.gamma not in self.ideal or not self.ray.contains(element):
            return False
        point = self.context.embed(element,self.scales)
        return ProductRegion([self.radius]*self.context.r1,[self.radius]*self.context.r2).contains(point)

    def sample(self,rng=None):
        rng = rng or RandomBits()
        v,attempts = self.grid.sample(rng,accept=lambda v: self.ray.contains(self._element(v)))
        element = self._element(v)
        if not self.verify(element):
            raise ArithmeticError("Algorithm 1 output failed exact validation")
        return BoxSample(element,self.scales,self.radius,tuple(v),{
            "algorithm":1,"uniform_conditional_on_completion":True,
            "exact_membership_verified":True,"geometric_hypotheses_verified":True,
            "attempts":attempts,"grid_bits":int(self.grid.certificate.denominator.nbits()),
            "basis_hkz_tours":self.short_basis.hkz_tours,
            "basis_method":self.short_basis.method,
            "paper_bit_complexity_certified":False,
            "random_bits":rng.bits_used,
        })


@dataclass(frozen=True)
class IdealSample:
    element: object
    input_ideal: object
    walk_ideal: object
    prime_factors: tuple
    log_distortion: tuple
    distortion: tuple
    grid_spacing: object
    box_sample: BoxSample
    report: dict

    def output_ideal(self):
        """(beta)*input_ideal^(-1), in the inverse input ideal class."""
        return self.input_ideal.number_field().ideal(self.element)/self.input_ideal


def algorithm2(context,ideal,walk,ray=None,y=None,epsilon=QQ(1)/100,block_size=2,
               omega=1,subgroup_oracle=None,rng=None,limits=Limits(),require_mixing_bound=True,
               prime_sampler=None):
    """Run Algorithm 2, conditional on the stated random-walk parameters.

    The default requires externally justified mixing information. Setting
    require_mixing_bound=False explicitly permits uncalibrated B,N experiments;
    their output must not be interpreted as satisfying Theorem 9.5's density bound.
    """
    if not isinstance(walk,WalkParameters):
        raise TypeError("walk must be WalkParameters")
    eps = rational(epsilon)
    delta = distortion_grid(context.n,context.discriminant,eps,omega,limits)
    conditional_mixing = walk.check(eps,require_mixing_bound)
    rng = rng or RandomBits()
    ray = ray or RayConditions(context)
    if ray.context is not context:
        raise ValueError("ray and algorithm must share a number-field context")
    verified_mixing=False
    if walk.certificate is not None:
        from .calibration import SpectralMixingCertificate
        if not isinstance(walk.certificate,SpectralMixingCertificate) or subgroup_oracle is not None:
            raise ValueError("spectral calibration currently certifies the full group without a custom oracle")
        walk.certificate.validate_for(context,ray,walk,eps)
        verified_mixing=True
    b = context.K.ideal(ideal)
    if not b or not ray.coprime(b):
        raise ValueError("a nonzero ideal coprime to the modulus is required")
    if walk.subgroup_index != 1 and subgroup_oracle is None:
        raise ValueError("a proper subgroup requires a membership oracle")
    if subgroup_oracle is not None and not subgroup_oracle(context.K.ideal(ray.tau)):
        raise ValueError("tau does not satisfy the supplied subgroup condition")
    y = context.scales(y)
    if prime_sampler is None:
        prime_sampler=PrimeIdealSampler(context,ray,walk.prime_bound,subgroup_oracle,limits)
    elif (prime_sampler.context is not context or prime_sampler.ray.modulus!=ray.modulus
          or prime_sampler.bound!=walk.prime_bound or prime_sampler.subgroup_oracle is not subgroup_oracle):
        raise ValueError("cached prime sampler does not match this random walk")
    bbar = b
    primes,prime_attempts = [],0
    for _ in range(int(walk.steps)):
        p,attempts = prime_sampler.sample(rng)
        bbar *= p
        primes.append(p)
        prime_attempts += attempts
    m = len(context.places)
    BH = matrix(QQ,m-1,m,lambda i,j: delta/context.n*(int(j==i)-int(j==i+1)))
    gaussian = discrete_gaussian(BH,QQ(1)/(context.n**2),eps/4,rng=rng,limits=limits)
    a = tuple(gaussian.value)
    A = norm_one_exp_approx(a,context.weights,delta/(2*context.n),limits)
    scales = tuple((Ai*re,Ai*im) for Ai,(re,im) in zip(A,y))
    prepared = Algorithm1(context,bbar,ray=ray,x=scales,block_size=block_size,omega=omega,limits=limits)
    sampled = prepared.sample(rng)
    beta = sampled.element
    if beta not in b or not ray.contains(beta):
        raise ArithmeticError("Algorithm 2 output failed exact validation")
    return IdealSample(beta,b,bbar,tuple(primes),a,A,delta,sampled,{
        "algorithm":2,"exact_membership_verified":True,
        "gaussian_sd_bound":str(eps/4),"distortion_error_controlled":True,
        "mixing_guarantee":("verified_spectral_bound" if verified_mixing else
                            "conditional_on_supplied_bound" if conditional_mixing else "not_established"),
        "mixing_justification":walk.justification,
        "theorem_9_5_density_bound":"conditional" if conditional_mixing else "not_established",
        "prime_bound":int(walk.prime_bound),"walk_steps":int(walk.steps),
        "prime_attempts":prime_attempts,"gaussian_attempts":gaussian.attempts,
        "paper_bit_complexity_certified":False,"random_bits":rng.bits_used,
    })
