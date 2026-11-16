import requests, random, collections, math
U="https://dna.cedar.host.dartmouth.edu"
def reads(strands,n):
    got=[]
    while len(got)<n:
        k=min(n-len(got),5000)
        got+=requests.post(U+"/api/reads",json={"key":"","strands":strands,"n":k},timeout=300).json()["reads"]
    return got
rng=random.Random(2)
rs=lambda L:"".join(rng.choice("ACGT") for _ in range(L))
L=200
S=[rs(L) for _ in range(40)]
R=reads(S,40000)
# locate exact reads via dict of all substrings? use find
loc=[]
for r in R:
    for j,s in enumerate(S):
        p=s.find(r)
        if p>=0: loc.append((j,p,len(r))); break
print("exact",len(loc),"of",len(R))
# end context stats: cut between s[e-1] and s[e]
for k in (1,2,3):
    cuts=collections.Counter(); occ=collections.Counter()
    for j,p,n in loc:
        e=p+n
        if 2<=e<=L-2 and e-k>=0 and e+k<=L: cuts[S[j][e-k:e+k]]+=1
        # also starts count as cuts
    # occurrences baseline
    for s in S:
        for e in range(k,L-k+1): occ[s[e-k:e+k]]+=1
    tot=sum(cuts.values()); to=sum(occ.values())
    rate={c:(cuts[c]/tot)/(occ[c]/to) for c in occ if occ[c]>=20}
    top=sorted(rate.items(),key=lambda x:-x[1])[:6]; bot=sorted(rate.items(),key=lambda x:x[1])[:4]
    print("k",k,"windows",len(rate),"top",[(a,round(b,2)) for a,b in top],"bot",[(a,round(b,2)) for a,b in bot])
# position bias of cuts
pos=collections.Counter()
for j,p,n in loc:
    pos[(p+n)//20]+=1
print("end pos by 20bin",sorted(pos.items()))
st=collections.Counter(p//20 for j,p,n in loc); print("start bins",sorted(st.items()))
# error: unlocated reads: find best match via k-mer seed + edit distance
def ed(a,b):
    prev=list(range(len(b)+1))
    for i,ca in enumerate(a,1):
        cur=[i]
        for jj,cb in enumerate(b,1): cur.append(min(prev[jj]+1,cur[-1]+1,prev[jj-1]+(ca!=cb)))
        prev=cur
    return prev[-1]
bad=[r for r in R if not any(r in s for s in S)][:300]
types=collections.Counter(); tot=0; errs=0
for r in bad:
    best=None
    for j,s in enumerate(S):
        for p in range(0,L-len(r)+3):
            pass
    # seed on first 6-mer or any 6-mer
    cands=set()
    for i in range(0,len(r)-5,3):
        sd=r[i:i+6]
        for j,s in enumerate(S):
            q=s.find(sd)
            while q>=0:
                cands.add((j,q-i)); q=s.find(sd,q+1)
    for j,st0 in cands:
        s=S[j]
        for d in range(-3,4):
            a=max(0,st0+d); w=s[a:a+len(r)+4]
            # semi-global: best over window ends
            val=min(ed(r,w[:m]) for m in range(max(1,len(r)-3),min(len(w),len(r)+3)+1))
            if best is None or val<best[0]: best=(val,j)
    if best: errs+=best[0]; tot+=len(r); types[best[0]]+=1
print("bad reads analysed",sum(types.values()),"edit-dist dist",sorted(types.items()),"per-base err in bad",errs/tot)
