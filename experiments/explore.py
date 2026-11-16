import requests, random, collections, json, math
U="https://dna.cedar.host.dartmouth.edu"
def reads(strands,n):
    got=[]
    while len(got)<n:
        k=min(n-len(got),5000)
        r=requests.post(U+"/api/reads",json={"key":"","strands":strands,"n":k},timeout=300).json()
        got+=r["reads"]
    return got
rng=random.Random(1)
def rs(L): return "".join(rng.choice("ACGT") for _ in range(L))
# 1. full-cover vs L, read length dist
for L in (20,30,40,60,120):
    S=[rs(L) for _ in range(50)]
    R=reads(S,10000)
    ln=[len(r) for r in R]
    full=sum(1 for r in R if len(r)==L)/len(R)
    sub=sum(1 for r in R if not any(r in s for s in S))
    print(L,"mean",sum(ln)/len(ln),"full",round(full,3),"min",min(ln),"max",max(ln),"nonsubstring",sub)
# 2. read length hist for L=300
S=[rs(300) for _ in range(30)]
R=reads(S,10000)
c=collections.Counter(len(r) for r in R)
print("L300 hist",sorted(c.items())[:40])
