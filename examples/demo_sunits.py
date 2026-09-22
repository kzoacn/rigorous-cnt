"""Compute Algorithm 6's S-unit system; reference sampling, independent GRH certificate."""

import argparse
import json
from pathlib import Path
from time import monotonic
from sage.all import QQ,PolynomialRing,NumberField
from rigorous_cnt import (NumberFieldContext,RandomBits,SUnitSettings,algorithm6,
                          ResourceLimit,calibrate_walk)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--field',choices=('real','imaginary','class2','mixed'),default='real')
    parser.add_argument('--primes',default='',help='comma-separated rational primes; include all primes above each')
    parser.add_argument('--smooth-bound',type=int,default=59)
    parser.add_argument('--generation-radius',type=int,default=2)
    parser.add_argument('--max-relations',type=int,default=100)
    parser.add_argument('--seed',type=int,default=2026)
    parser.add_argument('--modulus-prime',type=int,default=None,help='explicitly exercise the exceptional-unit branch')
    parser.add_argument('--calibrate-mixing',action='store_true',
                        help='optional explicit preprocessing using existing class/unit data, not used by default')
    parser.add_argument('--output',type=Path,default=None)
    args=parser.parse_args()
    R=PolynomialRing(QQ,'t');t=R.gen()
    f={'real':t*t-2,'imaginary':t*t+1,'class2':t*t+5,'mixed':t*t*t-2}[args.field]
    K=NumberField(f,'a');ctx=NumberFieldContext(K)
    S=[]
    for p in args.primes.split(','):
        if p.strip():S.extend(K.primes_above(int(p)))
    modulus=K.prime_above(args.modulus_prime) if args.modulus_prime is not None else None
    factory=(lambda context,ray,epsilon:calibrate_walk(context,ray,epsilon)) if args.calibrate_mixing else None
    settings=SUnitSettings(smooth_bound=args.smooth_bound,generation_radius=args.generation_radius,
        max_relations=args.max_relations,modulus_override=modulus,
        require_mixing_bound=args.calibrate_mixing,walk_factory=factory)
    started=monotonic()
    print(f'Field {f}; {len(S)} requested prime ideals; result certification assumes GRH.',flush=True)
    if not args.calibrate_mixing:
        print('Reference B,N: no mixing/time guarantee; complete generation is checked independently.',flush=True)
    def progress(info):
        print(f'{monotonic()-started:.1f}s: '+json.dumps(info),flush=True)
    result=algorithm6(ctx,S,settings=settings,rng=RandomBits(args.seed),progress=progress)
    print(json.dumps(result.report,indent=2),flush=True)
    print('Torsion:',result.torsion_order,result.torsion_generator,flush=True)
    print('Free rank:',result.rank,flush=True)
    try:
        print('Generators:',result.materialize_generators(max_bits=10000),flush=True)
    except ResourceLimit as error:
        print('Generators retained as compact power products:',error,flush=True)
    print('Exact representation checks:',result.verify(),flush=True)
    if args.output:
        args.output.write_text(json.dumps(result.to_dict(),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print('Saved exact data:',args.output,flush=True)


if __name__=='__main__':main()
