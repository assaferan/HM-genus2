import subprocess, re, json, os, sys
os.chdir('/scratch/home/assaferan/GitHub/HM-genus2')
rows=[json.loads(l) for l in open('sorted_output.jsonl')]
mod5={i for i,r in enumerate(rows) if r['congruence_prime']==5}
def midx(f,p):
    s=set()
    try:
        for l in open(f):
            m=re.match(r'ENTRY (\d+).*'+p,l)
            if m: s.add(int(m.group(1)))
    except FileNotFoundError: pass
    return s
matched=(midx('kernel_sweep_out.txt',r'-> MATCH')|midx('goal1_w24_results.txt',r'-> MATCH weight'))&mod5
nomatch24=midx('goal1_w24_results.txt',r'-> NO MATCH')&mod5   # already failed [2,4]/[4,2]
remaining=sorted(mod5-matched, key=lambda i: rows[i]['serre_conductor_2d'])
RES='smart_pass_results.txt'
done_smart=midx(RES,r'MATCH weight')|midx(RES,r'NOMATCH-ALL')
CAP=int(sys.argv[1]) if len(sys.argv)>1 else 1800
WALL_ALL=[(2,4),(4,2),(2,6),(6,2),(4,4),(4,6),(6,4),(6,6)]
WALL_HI =[(2,6),(6,2),(4,4),(4,6),(6,4),(6,6)]
print(f"remaining: {len(remaining)}  cap={CAP}s")
for idx in remaining:
    if idx in done_smart: continue
    weights = WALL_HI if idx in nomatch24 else WALL_ALL
    got=False
    for (a,b) in weights:
        try:
            out=subprocess.run(['magma','-b',f'idx:={idx}',f'k1:={a}',f'k2:={b}','match_one.m'],
                               capture_output=True,text=True,timeout=CAP).stdout
        except subprocess.TimeoutExpired:
            print(f"  idx {idx} [{a},{b}] TIMEOUT",flush=True); continue
        m=re.search(r'MATCH orbit.*',out)
        if m:
            line=f"ENTRY {idx} serre{rows[idx]['serre_conductor_2d']} MATCH weight [{a},{b}] {m.group(0)}"
            open(RES,'a').write(line+'\n'); print("  "+line,flush=True); got=True; break
        else:
            print(f"  idx {idx} [{a},{b}] nomatch",flush=True)
    if not got:
        line=f"ENTRY {idx} serre{rows[idx]['serre_conductor_2d']} NOMATCH-ALL-WEIGHTS (cap {CAP})"
        open(RES,'a').write(line+'\n'); print("  "+line,flush=True)
open(RES,'a').write('# smart pass done\n')
print("done")
