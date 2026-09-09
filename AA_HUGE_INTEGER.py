#!/usr/bin/env python3
"""IKERUSIKI Adaptive Address: implicit M_p=2**p-1. Python 3.9+, stdlib only.
Snapshot: https://www.mersenne.org/primes/ checked 2026-09-09.
No huge integer is expanded by default. --benchmark explicitly expands <=1M bits.
"""
import argparse, hashlib, json, math, platform, statistics, sys, time
from collections import defaultdict
from pathlib import Path

EXPONENTS = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]

def phi(groups):
    return sum(len(g)*(len(g)-1)//2 for g in groups)

def entropy(groups, n):
    out = 0.0
    for g in groups:
        p = len(g)/n
        out -= p*math.log2(p)
    return out

class Oracle:
    def __init__(self, exponents):
        self.exponents = exponents
        self.calls = 0
    def residue(self, i, modulus):
        self.calls += 1
        return (pow(2, self.exponents[i], modulus)-1) % modulus

def refine(oracle, groups, L, q):
    # Exact quotient residue without constructing n or delta:
    # ((n-ref) mod (L*q)) // L == ((n-ref)//L) mod q.
    if q == 1:
        return [list(g) for g in groups]
    out = []
    for g in groups:
        if len(g) <= 1:
            out.append(list(g)); continue
        residues = [oracle.residue(i, L*q) for i in g]
        ref = residues[0]
        buckets = defaultdict(list)
        for i, r in zip(g, residues):
            delta_mod = (r-ref) % (L*q)
            assert delta_mod % L == 0
            buckets[delta_mod//L].append(i)
        out.extend(buckets.values())
    return out

def signatures(exponents, axes):
    return [tuple((pow(2,p,b)-1)%b for b in axes) for p in exponents]

def run(exponents, initial_axes=(), axis_limit=512, max_axis=8192, mode='adaptive'):
    if not exponents or any(p<2 for p in exponents) or len(set(exponents)) != len(exponents):
        raise ValueError('Use distinct integer exponents >=2')
    t0=time.perf_counter()
    oracle=Oracle(exponents)
    selected=list(initial_axes)
    bins=defaultdict(list)
    for i in range(len(exponents)):
        bins[tuple(oracle.residue(i,b) for b in selected)].append(i)
    groups=list(bins.values())
    L=math.lcm(1,*selected)
    remaining=set(range(2,axis_limit+1))-set(selected)
    trajectory=[phi(groups)]; steps=[]; expansions=[]; evaluations=0
    while phi(groups)>0:
        before=phi(groups); cache={}; best=None
        for b in sorted(remaining):
            q=b//math.gcd(L,b)
            if q not in cache:
                cg=refine(oracle,groups,L,q)
                cache[q]=(phi(cg),entropy(cg,len(exponents)),cg)
            after,ent,cg=cache[q]; evaluations+=1
            row=((before-after,ent,-b),b,after,cg)
            if best is None or row[0]>best[0]: best=row
            if mode=='ascending': break
        if best is None or (mode=='adaptive' and best[2]==before):
            if axis_limit>=max_axis: break
            old=axis_limit; axis_limit=min(max_axis,axis_limit*2)
            remaining.update(range(old+1,axis_limit+1))
            expansions.append({'from':old,'to':axis_limit,'phi':before})
            continue
        score,b,after,newgroups=best
        # Refinement cannot merge any previously distinguished rows.
        parent={i:k for k,g in enumerate(groups) for i in g}
        assert all(len({parent[i] for i in g})==1 for g in newgroups)
        groups=newgroups; L=math.lcm(L,b); selected.append(b); remaining.remove(b)
        trajectory.append(after)
        steps.append({'axis':b,'before':before,'after':after,'entropy':score[1], 'oracle_calls_cumulative':oracle.calls})
    elapsed=time.perf_counter()-t0
    # Independent raw-residue signature check (outside selection timing).
    sig=signatures(exponents,selected)
    independent=defaultdict(list)
    for i,s in enumerate(sig): independent[s].append(i)
    assert sorted(map(sorted,groups))==sorted(map(sorted,independent.values()))
    return {'n':len(exponents),'max_exponent':str(max(exponents)), 'initial_axes':list(initial_axes),
            'selected_axes':selected,'added_axes':selected[len(initial_axes):], 'trajectory':trajectory,
            'depth':len(selected),'collision_free':phi(groups)==0,'elapsed_s':elapsed,
            'oracle_calls':oracle.calls,'candidate_evaluations':evaluations,'expansions':expansions,
            'steps':steps,'remaining_collision_exponents':[[str(exponents[i]) for i in g] for g in groups if len(g)>1],
            'addresses':[{'p':str(p),'address':list(s)} for p,s in zip(exponents,sig)]}

def measured(exponents, repeats, **kw):
    runs=[run(exponents,**kw) for _ in range(repeats)]
    r=runs[0]
    for x in runs: assert (x['selected_axes'],x['trajectory'])==(r['selected_axes'],r['trajectory'])
    r['timing']={'repeats':repeats,'median_s':statistics.median(x['elapsed_s'] for x in runs),
                 'min_s':min(x['elapsed_s'] for x in runs),'max_s':max(x['elapsed_s'] for x in runs)}
    # Cost of applying selected axes after they are known; selection is NOT included.
    ts=[]
    for _ in range(repeats):
        t=time.perf_counter(); signatures(exponents,r['selected_axes']); ts.append(time.perf_counter()-t)
    r['address_replay']={'pow_calls':len(exponents)*r['depth'],'median_s':statistics.median(ts)}
    return r

def benchmark(repeats):
    # Same residues, independent of adaptive selection; bounded explicit expansion.
    rows=[]
    for p in [32768,131072,1048576]:
        times={'modular':[],'expanded':[]}; object_bytes=0
        for _ in range(repeats):
            t=time.perf_counter(); a=[(pow(2,p,b)-1)%b for b in range(2,513)]; times['modular'].append(time.perf_counter()-t)
            t=time.perf_counter(); n=(1<<p)-1; c=[n%b for b in range(2,513)]; times['expanded'].append(time.perf_counter()-t)
            assert a==c; object_bytes=sys.getsizeof(n); del n
        rows.append({'p':p,'axes':511,'modular_median_s':statistics.median(times['modular']),
                     'expanded_build_plus_residues_median_s':statistics.median(times['expanded']),
                     'expanded_integer_object_bytes':object_bytes,'exponent_object_bytes':sys.getsizeof(p),'equal':True})
    return rows

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repeats',type=int,default=5)
    parser.add_argument('--benchmark',action='store_true')
    parser.add_argument('--output',default='results.json')
    args=parser.parse_args()
    if args.repeats<1: parser.error('repeats must be positive')
    cases={}
    cases['known52']=measured(EXPONENTS,args.repeats)
    # Consecutive exponents, explicitly generated; no primality claim or preselected axes.
    huge=[10**100+i for i in range(52)]
    cases['beyond52']=measured(huge,args.repeats)
    cases['known52_ascending']=measured(EXPONENTS,args.repeats,mode='ascending')
    cases['beyond52_ascending']=measured(huge,args.repeats,mode='ascending')
    axes=cases['known52']['selected_axes']
    cases['retain_same52']=measured(EXPONENTS,args.repeats,initial_axes=axes)
    cases['extend104']=measured(EXPONENTS+huge,args.repeats,initial_axes=axes)
    # Deliberate Blind Spot: all bases <=16 collide for these two exponents.
    period=math.lcm(*range(1,17))
    cases['expand_control']=measured([1000,1000+period],args.repeats,axis_limit=16,max_axis=512)
    result={'source_date':'2026-09-09','gimps_source':'https://www.mersenne.org/primes/',
            'upstream_blob_sha':'f00a6160787750aae76e18791f8f31600ff32834',
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'environment':{'python':sys.version,'platform':platform.platform(),'machine':platform.machine()},
            'protocol':{'initial_candidate_axes':'2..512','select':'max reduction, max float entropy, min b',
                        'residue':'M_p mod b = (pow(2,p,b)-1)%b','huge_exponents':'10**100+i, i=0..51',
                        'timing_excludes':'independent verification, JSON I/O, corpus definition; replay timed separately',
                        'expansion':'double candidate upper bound only on no progress; cap 8192; report unresolved'},
            'cases':cases,'representation_benchmark':benchmark(args.repeats) if args.benchmark else []}
    Path(args.output).write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    for name,r in cases.items():
        print(name,'axes=',r['selected_axes'],'phi=',r['trajectory'],'median_s=',round(r['timing']['median_s'],6),'pow_calls=',r['oracle_calls'])
    assert all(r['collision_free'] for r in cases.values())
    print('PASS; saved',args.output)
if __name__=='__main__': main()
