#!/usr/bin/env python3
"""Stage 6, Round 2: exact verification of the universal trace–Gram identity.

This optional certificate verifies equation (2.4) and the expansion table in
阶段6第2轮.txt, for the formal variable m. It uses only Python's standard
library and exact rational coefficients. There are no sampled matrices,
floating-point calculations, optimizers, or finite-m extrapolations.

A word is traced and summed over its matrix-index labels. Cyclic permutations
and reversal preserve its real trace. Singleton labels can be removed using
sum A_i = I; an absent summed label contributes m. Distinct-index sums are
expanded by the partition-lattice inclusion–exclusion formula proved in the
report. Positive semidefiniteness of the Gram family is proved analytically
in the report and is not delegated to this program.

Run: python3 stage6_round2_identity_certificate.py
Expected: PASS: symbolic identity for indeterminate m
"""
from fractions import Fraction as F
from collections import defaultdict,Counter
from math import factorial

def padd(a,b):
 z=[F(0)]*max(len(a),len(b))
 for i,c in enumerate(a):z[i]+=c
 for i,c in enumerate(b):z[i]+=c
 while z and not z[-1]:z.pop()
 return tuple(z)
def pmul(a,b):
 if not a or not b:return ()
 z=[F(0)]*(len(a)+len(b)-1)
 for i,c in enumerate(a):
  for j,d in enumerate(b):z[i+j]+=c*d
 while z and not z[-1]:z.pop()
 return tuple(z)
def ps(a,c):return tuple(x*c for x in a)
def partitions(n):
 if n==0:yield []
 else:
  for p in partitions(n-1):
   yield p+[[n-1]]
   for i in range(len(p)):
    yield p[:i]+[p[i]+[n-1]]+p[i+1:]
def normlabels(w):
 d={};z=[]
 for i in w:
  if i=='X':z.append('X')
  else:
   if i not in d:d[i]=chr(65+len(d))
   z.append(d[i])
 return ''.join(z)
def normalize(w,r):
 counts=Counter(w);absent=sum(i not in counts for i in range(r))
 w=tuple(i for i in w if i=='X' or counts[i]>1)
 return min(normlabels(z[i:]+z[:i]) for z in [w,w[::-1]] for i in range(len(w))),absent

def qword(p):
 if p[0]=='E':return {(p[1]+1,'X',p[1]+1):1,(0,'X',p[1]+1):-1}
 if p[0]=='F':return {(p[1]+1,'X',p[1]+1):1,(p[1]+1,'X',0):-1}
 j,l=p[1]+1,p[2]+1
 return {(j,'X',l):1,(0,'X',l):-1,(j,'X',0):-1,(0,'X',0):1}
def pattern(p,q):
 labels={};out=[]
 for z in [p,q]:
  a=[z[0]]
  for j in z[1:]:
   if j not in labels:labels[j]=len(labels)
   a.append(labels[j])
  out.append(tuple(a))
 return tuple(out)
def coef(p,q):
 if p[0]>q[0]:p,q=q,p
 a,b=p[0],q[0]
 if a==b=='E':return (F(-2),F(3,2),F(1,2)) if p[1]==q[1] else (F(-2),F(1,2))
 if a==b=='F':return (F(-2),F(3,2),F(1,2)) if p[1]==q[1] else (F(-2),F(5,2))
 if (a,b)==('E','F'):return (F(0),F(-1,2),F(-1,2)) if p[1]==q[1] else (F(0),F(-1))
 if (a,b)==('E','Z'):return (F(-1,2),F(1,2)) if p[1]==q[1] else (F(-1,2),)
 if (a,b)==('F','Z'):
  if p[1]==q[1]:return (F(1),F(-1,2))
  if p[1]==q[2]:return ()
  return (F(1,2),)
 if a==b=='Z':
  if p==q:return (F(-2),F(1))
  if p[1]==q[1] or p[2]==q[2]:return (F(-1,2),)
  return ()
 raise Exception()

labs=[('E',j) for j in range(4)]+[('F',j) for j in range(4)]+[('Z',j,l) for j in range(4) for l in range(4) if j!=l]
patterns=sorted(set(pattern(p,q) for p in labs for q in labs))
groups={z:defaultdict(tuple) for z in ['EE','EF','EZ','FF','FZ','ZZ']}
for p,q in patterns:
 coeff=coef(p,q)
 if not coeff:continue
 r=max(p[1:]+q[1:])+2
 total=defaultdict(tuple)
 for part in partitions(r):
  mu=1;mp={}
  for i,block in enumerate(part):
   mu*=(-1)**(len(block)-1)*factorial(len(block)-1)
   mp.update({j:i for j in block})
  for w,a in qword(p).items():
   for v,b in qword(q).items():
    word=w[::-1]+(0,)+v
    word=tuple('X' if j=='X' else mp[j] for j in word)
    key,power=normalize(word,len(part))
    pol=ps((F(0),)*power+coeff,mu*a*b)
    total[key]=padd(total[key],pol)
 group=''.join(sorted([p[0],q[0]]))
 for key,pol in total.items():groups[group][key]=padd(groups[group][key],pol)
allkeys=sorted(set().union(*[set(g) for g in groups.values()]))
sumout=defaultdict(tuple)
for key in allkeys:
 for g in groups.values():sumout[key]=padd(sumout[key],g[key])
expected={normalize(('X','X',0,0),1)[0]:(F(1),),normalize(('X',0,0,0,'X',1,1),2)[0]:(F(0),F(-1),F(1)),normalize((0,'X',0,1,'X',1),2)[0]:(F(1),F(-2))}
assert {key:p for key,p in sumout.items() if p}==expected
print('PASS: symbolic identity for indeterminate m')
