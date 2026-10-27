#!/usr/bin/env python3
"""Rotational quotient of the all-chord arrangement (stdlib only).

Exact cyclotomic points/areas; exploratory Decimal ordering at 60 and 120
 digits. No complete orbit of a noncentral intersection is stored.
"""
import argparse
from collections import Counter
from decimal import Decimal as D, localcontext
from functools import cmp_to_key
from itertools import combinations
import hashlib
import json
import math
from pathlib import Path
import resource
import sys
import time

from cyclotomic import Field, wedge
from explore_cells import Embedding, digest
from verify_geometry import pi, sin_cos

SCHEMA = 'rotational-quotient-v1'


class Arrangement:
    def __init__(self, n, precisions=(60,120)):
        self.n, self.f = n, Field(n)
        self.precisions = tuple(precisions)
        self.reps = [self.f.one, self.f.zero]  # orbit 1 is the fixed origin
        self.lookup = {self.f.one:0}
        self.center = False
        self.cuts = {s:{0,s} for s in range(1,n//2+1)}
        self.embeddings, self.roots = {}, {}
        for precision in precisions:
            with localcontext() as ctx:
                ctx.prec = precision
                self.embeddings[precision] = Embedding(self.f,precision)
                p = pi()
                self.roots[precision] = [tuple(reversed(sin_cos(D(2*j)*p/n,p)))
                                         for j in range(n)]
        self.quadruple_orbits = 0

    def shift(self, vertex, amount):
        if vertex == self.n:
            return vertex
        orbit, phase = divmod(vertex,self.n)
        return orbit*self.n+(phase+amount)%self.n

    def point(self, vertex):
        orbit, phase = divmod(vertex,self.n)
        return self.reps[orbit].rotate(phase)

    def sector(self, z, precision):
        with localcontext() as ctx:
            ctx.prec = precision
            x,y = self.embeddings[precision].point(z)
            guess = math.floor((math.atan2(float(y),float(x))%(2*math.pi))
                              *self.n/(2*math.pi))%self.n
            tolerance = D(10)**(-(precision-20))
            def sign(value, first, second):
                if abs(value)>tolerance:
                    return 1 if value>0 else -1
                if wedge(first,second)==0:
                    return 0
                raise ArithmeticError('unresolved rotational sector')
            for k in (guess,(guess-1)%self.n,(guess+1)%self.n):
                a,b = self.roots[precision][k],self.roots[precision][(k+1)%self.n]
                left = sign(a[0]*y-a[1]*x,self.f.root(k),z)
                right = sign(x*b[1]-y*b[0],z,self.f.root(k+1))
                if left>=0 and right>0:
                    return k
            raise ArithmeticError('point not assigned to a sector')

    def intern(self,z):
        if not z:
            self.center = True
            return self.n
        phases = [self.sector(z,p) for p in self.precisions]
        if len(set(phases))!=1:
            raise ArithmeticError('precision-dependent rotational sector')
        phase = phases[0]
        canonical = z.rotate(-phase)
        if canonical not in self.lookup:
            self.lookup[canonical] = len(self.reps)
            self.reps.append(canonical)
        return self.lookup[canonical]*self.n+phase

    def build(self):
        n = self.n
        for a,b,c in combinations(range(1,n),3):
            gaps = (a,b-a,c-b,n-c)
            if gaps != min(gaps[j:]+gaps[:j] for j in range(4)):
                continue
            self.quadruple_orbits += 1
            v = self.intern(self.f.intersection(0,b,a,c))
            for x,y in ((0,b),(a,c)):
                if y-x>n//2:
                    x,y = y,x
                length = (y-x)%n
                self.cuts[length].add(self.shift(v,-x))
                if 2*length==n:
                    self.cuts[length].add(self.shift(v,-y))
        if self.center != (n%2==0):
            raise ArithmeticError('unexpected central vertex')

    def dart(self,u,v):
        if u==self.n:
            return (1,(v//self.n)*self.n)
        return (u//self.n,self.shift(v,-(u%self.n)))

    def topology(self,precision):
        n = self.n
        with localcontext() as ctx:
            ctx.prec = precision
            tolerance = D(10)**(-(precision-20))
            embedding = self.embeddings[precision]
            rep_xy = [embedding.point(p) for p in self.reps]
            xy_cache = {}
            def xy(v):
                if v not in xy_cache:
                    orbit,phase = divmod(v,n)
                    x,y = rep_xy[orbit]; c,s = self.roots[precision][phase]
                    xy_cache[v] = (x*c-y*s,x*s+y*c)
                return xy_cache[v]
            adjacent, darts = {},set()
            for length,vertices in self.cuts.items():
                axis = max((0,1),key=lambda j:abs(xy(length)[j]-xy(0)[j]))
                ordered = sorted(vertices,key=lambda v:xy(v)[axis])
                for u,v in zip(ordered,ordered[1:]):
                    if xy(v)[axis]-xy(u)[axis]<=tolerance:
                        raise ArithmeticError('unresolved chord ordering')
                    first,reverse = self.dart(u,v),self.dart(v,u)
                    if first==reverse:
                        raise ArithmeticError('edge fixed by a nontrivial rotation')
                    for orbit,target in (first,reverse):
                        darts.add((orbit,target))
                        adjacent.setdefault(orbit,set()).add(target)
            if not {(0,1),(0,n-1)}<=darts:
                raise ArithmeticError('polygon boundary cycle missing')
            if self.center:
                adjacent[1] = {v+j for v in adjacent[1] for j in range(n)}
            successor = {}
            def sign(value):
                return 0 if abs(value)<=tolerance else (1 if value>0 else -1)
            for orbit,neighbors in adjacent.items():
                u = orbit*n
                ux,uy = xy(u)
                def compare(a,b):
                    ax,ay = xy(a)[0]-ux,xy(a)[1]-uy
                    bx,by = xy(b)[0]-ux,xy(b)[1]-uy
                    ha = 0 if sign(ay)>0 or (sign(ay)==0 and sign(ax)>0) else 1
                    hb = 0 if sign(by)>0 or (sign(by)==0 and sign(bx)>0) else 1
                    if ha!=hb:
                        return -1 if ha<hb else 1
                    cross = ax*by-ay*bx
                    if abs(cross)<=tolerance:
                        raise ArithmeticError('unresolved cyclic ordering')
                    return -1 if cross>0 else 1
                neighbors = sorted(neighbors,key=cmp_to_key(compare))
                for i,v in enumerate(neighbors):
                    successor[orbit,v] = neighbors[(i-1)%len(neighbors)]
            unseen = set(darts)
            faces = []
            for start in sorted(darts):
                if start not in unseen:
                    continue
                u,v = start[0]*n,start[1]
                edge = (u,v); visited=set();face=[];keys=[]
                while True:
                    if edge in visited:
                        if edge!=(u,v):
                            raise ArithmeticError('nonclosed face orbit')
                        break
                    visited.add(edge)
                    a,b = edge
                    face.append(a)
                    key = self.dart(a,b)
                    if key not in darts:
                        raise ArithmeticError('unknown dart')
                    keys.append(key)
                    if b==n:
                        following = successor[1,a]
                    else:
                        phase = b%n
                        following = self.shift(successor[b//n,self.shift(a,-phase)],phase)
                    edge = b,following
                if len(face)<3 or len(face)!=len(set(face)):
                    raise ArithmeticError('nonsimple face')
                counts=Counter(keys)
                if len(set(counts.values()))!=1:
                    raise ArithmeticError('inconsistent face stabilizer')
                stabilizer = next(iter(counts.values()))
                if n%stabilizer:
                    raise ArithmeticError('invalid face stabilizer')
                if not set(keys)<=unseen:
                    raise ArithmeticError('overlapping face orbits')
                unseen.difference_update(keys)
                faces.append((face,n//stabilizer))
            if unseen:
                raise ArithmeticError('unvisited darts')
            # Explicit quotient connectivity, including the fixed origin when present.
            seen=set();pending=[0]
            while pending:
                orbit=pending.pop()
                if orbit not in seen:
                    seen.add(orbit)
                    pending.extend(v//n for v in adjacent[orbit])
            if seen!=set(adjacent):
                raise ArithmeticError('disconnected quotient graph')
            vertices=(len(self.reps)-1)*n+int(self.center)
            edges=len(darts)*n//2
            face_count=sum(weight for face,weight in faces)
            if vertices-edges+face_count!=2:
                raise ArithmeticError('weighted Euler characteristic failed')
            return {'vertices':vertices,'edges':edges,'cells':face_count-1,
                    'vertex_orbits':len(self.reps)-1+int(self.center),
                    'edge_orbits':len(darts)//2,'face_orbits_including_exterior':len(faces),
                    'quotient_topology_sha256':digest(faces)},faces,embedding

    def areas(self,faces,embedding):
        normalizer = (self.n*(self.f.root(1)-self.f.root(-1))).inverse()
        edge_areas={};exact_points={};total=self.f.zero
        exterior=0;rationals=Counter();records=[];bounded_orbits=0;multiset=Counter()
        def point(v):
            if v not in exact_points:
                exact_points[v]=self.point(v)
            return exact_points[v]
        for face,weight in faces:
            area=self.f.zero
            for u,v in zip(face,face[1:]+face[:1]):
                key=self.dart(u,v)
                if key not in edge_areas:
                    a,b=key[0]*self.n,key[1]
                    edge_areas[key]=wedge(point(a),point(b))*normalizer
                area+=edge_areas[key]
            if area==-1:
                if weight!=1:
                    raise ArithmeticError('exterior orbit multiplicity')
                exterior+=1
                continue
            for precision in self.precisions:
                with localcontext() as ctx:
                    ctx.prec=precision
                    real,imag=self.embeddings[precision].point(area)
                    tolerance=D(10)**(-(precision-20))
                    if real<=tolerance or abs(imag)>tolerance:
                        raise ArithmeticError('nonpositive/nonreal face area')
            total+=area*weight
            q=area.rational()
            if q is not None:
                rationals[str(q)]+=weight
            bounded_orbits+=1
            records.append((area.key(),weight))
            multiset[area.key()] += weight
        if exterior!=1 or total!=1:
            raise ArithmeticError('exact weighted area sum or exterior failed')
        return {'exact_area_sum':'1','euler_characteristic':2,
                'rational_cells':sum(rationals.values()),
                'rational_area_counts':dict(sorted(rationals.items())),
                'bounded_cell_orbits':bounded_orbits,
                'weighted_cell_areas_sha256':digest(records),
                'exact_area_multiset_sha256':digest(sorted((key,count) for key,count in multiset.items()))}


def calculate(n,precisions=(60,120)):
    with localcontext() as ctx:
        ctx.prec=max(precisions)
        started=time.perf_counter()
        a=Arrangement(n,precisions);a.build()
        built=time.perf_counter();topology=None;faces=None;embedding=None
        for precision in precisions:
            report,these_faces,these_embedding=a.topology(precision)
            if topology is not None and report!=topology:
                raise ArithmeticError('precision-dependent quotient topology')
            topology,faces,embedding=report,these_faces,these_embedding
        ordered=time.perf_counter()
        areas=a.areas(faces,embedding)
        finished=time.perf_counter()
        report={'n':n,'quadruple_orbits':a.quadruple_orbits,**topology,**areas}
        measurements={'n':n,'build_seconds':built-started,'ordering_seconds':ordered-built,
                      'area_seconds':finished-ordered,'elapsed_seconds':finished-started,
                      'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                                       *(1 if sys.platform=='darwin' else 1024)}
        return report,measurements


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,nargs='+',default=list(range(3,101)))
    parser.add_argument('--checkpoint-dir',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--benchmark-output',type=Path)
    args=parser.parse_args()
    if any(n<3 for n in args.n) or len(set(args.n))!=len(args.n):
        parser.error('--n requires distinct integers at least 3')
    results=[];measurements=[]
    import cyclotomic, explore_cells, verify_geometry, verify_myerson
    dependencies = [Path(__file__),Path(cyclotomic.__file__),Path(explore_cells.__file__),
                    Path(verify_geometry.__file__),Path(verify_myerson.__file__)]
    code_hash=digest({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies})
    if args.checkpoint_dir:
        args.checkpoint_dir.mkdir(parents=True,exist_ok=True)
    for n in args.n:
        checkpoint=args.checkpoint_dir/f'n-{n}.json' if args.checkpoint_dir else None
        if checkpoint and checkpoint.exists():
            saved=json.loads(checkpoint.read_text())
            if saved['code_sha256']!=code_hash:
                raise SystemExit('Stale checkpoint: '+str(checkpoint))
            report,measurement=saved['report'],saved['measurement']
        else:
            report,measurement=calculate(n)
            if checkpoint:
                temporary=checkpoint.with_suffix('.tmp')
                temporary.write_text(json.dumps({'code_sha256':code_hash,'report':report,
                                                 'measurement':measurement},indent=2)+'\n')
                temporary.replace(checkpoint)
        results.append(report);measurements.append(measurement)
        print(f'Orbit cells: n={n}, cells={report["cells"]}, rational={report["rational_cells"]}, '
              f'seconds={measurement["elapsed_seconds"]:.2f}',file=sys.stderr,flush=True)
    report={'method':SCHEMA,'certified_arrangement':False,'ordering_precision_digits':[60,120],
            'precision_results_agree':True,'polygons':results}
    payload=json.dumps(report,indent=2)+'\n'
    if args.output:
        args.output.write_text(payload)
    else:
        print(payload,end='')
    if args.benchmark_output:
        args.benchmark_output.write_text(json.dumps({'peak_rss_scope':'process lifetime high-water mark',
                                                     'measurements':measurements},indent=2)+'\n')


if __name__=='__main__':
    main()
