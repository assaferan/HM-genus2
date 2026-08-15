import subprocess, re, json, os, sys
os.chdir('/scratch/home/assaferan/GitHub/HM-genus2')
rows=[json.loads(l) for l in open('sorted_output.jsonl')]
mod5={i for i,r in enumerate(rows) if r['congruence_prime']==5}
def midx(f,p):
    s=set()
    try:
        for l in open(f):
            m=re.match(r'ENTRY (\d+).*'+p,l); 
            if m: s.add(int(m.group(1)))
    except FileNotFoundError: pass
    return s
matched=(midx('kernel_sweep_out.txt',r'-> MATCH')|midx('goal1_w24_results.txt',r'-> MATCH weight'))&mod5
remaining=sorted(mod5-matched, key=lambda i: rows[i]['serre_conductor_2d'])
RES='kernel_wt_results.txt'
done=set()
try:
    for l in open(RES):
        m=re.match(r'ENTRY (\d+)',l)
        if m: done.add(int(m.group(1)))
except FileNotFoundError: pass
WEIGHTS=[(2,4),(4,2),(2,6),(6,2),(4,4),(4,6),(6,4),(6,6)]
CAP=500
print(f"remaining {len(remaining)}",flush=True)
for idx in remaining:
    if idx in done: continue
    serre=rows[idx]['serre_conductor_2d']; got=None; skipbig=False
    for (a,b) in WEIGHTS:
        try:
            out=subprocess.run(['magma','-b',f'idx:={idx}',f'k1:={a}',f'k2:={b}','kernel_wt.m'],
                               capture_output=True,text=True,timeout=CAP).stdout
        except subprocess.TimeoutExpired:
            print(f"  idx{idx}[{a},{b}] TIMEOUT",flush=True); continue
        if 'SKIP-BIG' in out: skipbig=True; print(f"  idx{idx}[{a},{b}] skipbig",flush=True); continue
        m=re.search(r'kernel: (MATCH) \(true (\d+) / ctrl (\d+)\)',out)
        if m and int(m.group(2))>0 and int(m.group(3))==0:
            got=(a,b,m.group(2)); break
        print(f"  idx{idx}[{a},{b}] {'nomatch' if 'NO-MATCH' in out else '?'}",flush=True)
    if got:
        line=f"ENTRY {idx} serre{serre} MATCH weight [{got[0]},{got[1]}] kernel true={got[2]}"
    elif skipbig:
        line=f"ENTRY {idx} serre{serre} SKIP-BIG (dim>cap at raised weight)"
    else:
        line=f"ENTRY {idx} serre{serre} NO-MATCH (tried [2,4],[4,2],[2,6],[6,2])"
    open(RES,'a').write(line+'\n'); print(" "+line,flush=True)
open(RES,'a').write('# kernel_wt pass done\n'); print("done",flush=True)
