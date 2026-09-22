"""Small runnable examples; the Algorithm 2 parameters here are uncalibrated."""

import argparse
import json
from sage.all import QQ,PolynomialRing,NumberField
from rigorous_cnt import NumberFieldContext,RayConditions,RandomBits,Algorithm1,algorithm2,WalkParameters,calibrate_walk


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--field',choices=('real','imaginary','mixed'),default='real')
    parser.add_argument('--algorithm',choices=('1','2','both'),default='both')
    parser.add_argument('--seed',type=int,default=2026)
    parser.add_argument('--prime-bound',type=int,default=19)
    parser.add_argument('--walk-steps',type=int,default=None,
                        help='override N; default uses the certified volume upper bound')
    parser.add_argument('--calibrate',action='store_true',
                        help='perform explicit spectral preprocessing to verify the mixing bound')
    args=parser.parse_args()
    ring=PolynomialRing(QQ,'t'); t=ring.gen()
    f={'real':t*t-2,'imaginary':t*t+1,'mixed':t*t*t-2}[args.field]
    K=NumberField(f,'a'); ctx=NumberFieldContext(K)
    ray=RayConditions(ctx,K.ideal(2),real_places=range(ctx.r1),tau=-1)
    ideal=K.ideal(3)
    print(f'Field: {f}; signature: {(ctx.r1,ctx.r2)}; modulus: (2); ideal: (3)')
    if args.algorithm in ('1','both'):
        sampler=Algorithm1(ctx,ideal,ray=ray)
        sample=sampler.sample(RandomBits(args.seed))
        print('Algorithm 1 alpha:',sample.element)
        print(json.dumps(sample.report,indent=2))
    if args.algorithm in ('2','both'):
        if args.calibrate:
            if args.walk_steps is not None:
                parser.error('--walk-steps cannot override certified calibration')
            walk=calibrate_walk(ctx,ray,initial_prime_bound=args.prime_bound)
            print('Algorithm 2: explicit spectral mixing bound verified; preprocessing uses class/unit data.')
        else:
            walk=(WalkParameters.from_prime_bound(ctx,ray,args.prime_bound) if args.walk_steps is None
                  else WalkParameters(args.prime_bound,args.walk_steps))
            print('Algorithm 2: B is experimental; the mixing bound is NOT established.')
        result=algorithm2(ctx,ideal,walk,ray=ray,
                          rng=RandomBits(args.seed),require_mixing_bound=args.calibrate)
        print('Algorithm 2 beta:',result.element)
        print('Walk prime norms:',[int(p.norm()) for p in result.prime_factors])
        print(json.dumps(result.report,indent=2))


if __name__=='__main__': main()
