"""Compatibility and invariant checks. Bounded integers only (<=131072 bits)."""
import hashlib,json,math
from pathlib import Path
import upstream_reproducer as ref
import AA_HUGE_INTEGER as new

source=Path('upstream_reproducer.py').read_bytes()
blob=hashlib.sha1(b'blob '+str(len(source)).encode()+b'\0'+source).hexdigest()
assert blob=='f00a6160787750aae76e18791f8f31600ff32834'
class IntegerOracle:
    def __init__(self,values): self.values=values
    def residue(self,i,m): return self.values[i]%m

def compare(values,oracle):
    groups=[list(range(len(values)))];L=1;remaining=set(range(2,513));axes=[];trajectory=[new.phi(groups)];checks=0
    while new.phi(groups):
        before=new.phi(groups);best=None
        for b in sorted(remaining):
            q=b//math.gcd(L,b)
            a=ref.refine_by_q(values,groups,L,q)
            c=new.refine(oracle,groups,L,q)
            assert a==c
            assert ref.entropy_of_partition(a,len(values))==new.entropy(c,len(values))
            checks+=1
            row=((before-new.phi(c),new.entropy(c,len(values)),-b),b,c)
            if best is None or row[0]>best[0]:best=row
        score,b,groups=best
        assert score[0]>0
        axes.append(b);remaining.remove(b);L=math.lcm(L,b);trajectory.append(new.phi(groups))
    return {'axes':axes,'trajectory':trajectory,'candidate_partition_checks':checks}

out={'upstream_blob_sha':blob,'checks':[]}
for bits in [32768,65536,131072]:
    values=ref.generate_corpus(bits,100,20260812)
    r=compare(values,IntegerOracle(values))
    expected=ref.EXPECTED[(bits,20260812)]
    assert '|'.join(map(str,r['axes']))==expected[0]
    assert '->'.join(map(str,r['trajectory']))==expected[1]
    out['checks'].append({'bits':bits,'seed':20260812,**r})
p=new.EXPONENTS[:23]
r=compare([(1<<x)-1 for x in p],new.Oracle(p))
assert r['axes']==new.run(p)['selected_axes']
out['checks'].append({'small_mersenne_max_p':max(p),**r})
# Folded distance is deliberately absent. Both raw n%b and shifted (n-1)%b
# give the same equality partitions; coordinates themselves are different.
for p in range(2,80):
    for b in range(2,513): assert (pow(2,p,b)-1)%b==((1<<p)-1)%b
out['small_residue_checks']=78*511
r=new.run([1000,1000+math.lcm(*range(1,17))],axis_limit=16,max_axis=16)
assert not r['collision_free'] and not r['added_axes']
out['bounded_blind_spot_reported']=True
out['status']='PASS'
Path('validation.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
