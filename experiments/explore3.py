import requests, random, collections, pickle
U="https://dna.cedar.host.dartmouth.edu"
def reads(strands,n):
    got=[]
    while len(got)<n:
        k=min(n-len(got),5000)
        got+=requests.post(U+"/api/reads",json={"key":"","strands":strands,"n":k},timeout=300).json()["reads"]
    return got
rng=random.Random(3)
rs=lambda L:"".join(rng.choice("ACGT") for _ in range(L))
L=200
S=[rs(L) for _ in range(100)]
R=[]
for _ in range(5): R+=reads(S,40000)
idx={}
loc=[]
for r in R:
    for j,s in enumerate(S):
        p=s.find(r)
        if p>=0: loc.append((j,p,len(r))); break
pickle.dump((S,R,loc),open("data.pkl","wb"))
print("exact",len(loc),len(R))
# cut table by (left a bases, right b bases)
def table(a,b):
    cuts=collections.Counter(); occ=collections.Counter()
    for j,p,n in loc:
        e=p+n
        if a<=e<=L-b: cuts[S[j][e-a:e+b]]+=1
    for s in S:
        for e in range(a,L-b+1): occ[s[e-a:e+b]]+=1
    tot=sum(cuts.values()); to=sum(occ.values())
    return {c:(cuts[c]/tot)/(occ[c]/to) for c in occ if occ[c]>=30}
for a,b in ((1,0),(2,0),(3,0),(0,1),(0,2),(0,3),(2,2),(3,1)):
    t=table(a,b)
    srt=sorted(t.items(),key=lambda x:-x[1])
    vals=sorted(t.values())
    print(f"left{a} right{b}: n={len(t)} range {vals[0]:.2f}-{vals[-1]:.2f} top",[(c,round(v,1)) for c,v in srt[:6]],"bot",[(c,round(v,2)) for c,v in srt[-4:]])
