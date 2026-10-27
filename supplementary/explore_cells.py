#!/usr/bin/env python3
"""All-diagonal cells: exact intersections/areas, numerical planar ordering.

This is exploratory, not a certified arrangement algorithm: Decimal embeddings
order intersections and incident edges. Both 60- and 120-digit runs must agree.
"""
import argparse
from collections import Counter
from decimal import Decimal as D, localcontext
from functools import cmp_to_key
import hashlib
from itertools import combinations
import json
import sys
from cyclotomic import Field, wedge
from verify_geometry import pi, sin_cos


class Embedding:
    def __init__(self, field, precision):
        self.field, self.precision = field, precision
        p = pi()
        self.powers = [tuple(reversed(sin_cos(D(2*j)*p/field.order,p)))
                       for j in range(field.degree)]

    def point(self, z):
        return tuple(sum(D(c)*self.powers[j][axis] for j,c in enumerate(z.coeff))/D(z.den)
                     for axis in (0,1))


def digest(value):
    return hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()


def build_arrangement(n):
    f = Field(n)
    points = [f.root(j) for j in range(n)]
    lookup = {p:j for j,p in enumerate(points)}
    chords = {(a,b):{a,b} for a,b in combinations(range(n),2)}
    for a,b,c,d in combinations(range(n),4):
        z = f.intersection(a,c,b,d)
        if z not in lookup:
            lookup[z] = len(points)
            points.append(z)
        vertex = lookup[z]
        chords[a,c].add(vertex)
        chords[b,d].add(vertex)
    return f, points, chords


def enumerate_faces(f, points, chords, precision, area_multiset=None):
    with localcontext() as context:
        context.prec = precision
        embedding = Embedding(f,precision)
        xy = [embedding.point(z) for z in points]
        tolerance = D(10)**(-(precision-20))
        edges = set()
        for (a,b), vertices in chords.items():
            axis = max((0,1),key=lambda j:abs(xy[b][j]-xy[a][j]))
            ordered = sorted(vertices,key=lambda j:xy[j][axis])
            for u,v in zip(ordered,ordered[1:]):
                if xy[v][axis]-xy[u][axis] <= tolerance:
                    raise ArithmeticError('unresolved ordering of distinct exact intersections')
                edges.add(tuple(sorted((u,v))))
        adjacent = [[] for _ in points]
        for u,v in sorted(edges):
            adjacent[u].append(v)
            adjacent[v].append(u)
        def sign(x):
            return 0 if abs(x) <= tolerance else (1 if x > 0 else -1)
        for u,neighbors in enumerate(adjacent):
            def compare(a,b):
                ax,ay = (xy[a][j]-xy[u][j] for j in (0,1))
                bx,by = (xy[b][j]-xy[u][j] for j in (0,1))
                ha = 0 if sign(ay)>0 or (sign(ay)==0 and sign(ax)>0) else 1
                hb = 0 if sign(by)>0 or (sign(by)==0 and sign(bx)>0) else 1
                if ha != hb:
                    return -1 if ha<hb else 1
                cross = ax*by-ay*bx
                if abs(cross) <= tolerance:
                    raise ArithmeticError('unresolved cyclic order at exact vertex')
                return -1 if cross>0 else 1
            neighbors.sort(key=cmp_to_key(compare))
        successor = {}
        for v,neighbors in enumerate(adjacent):
            for i,u in enumerate(neighbors):
                successor[u,v] = (v,neighbors[(i-1)%len(neighbors)])
        # Traversing the edge preceding the incoming edge keeps the face left.
        remaining = set(successor)
        faces = []
        for start in successor:
            if start not in remaining:
                continue
            edge, face = start, []
            while True:
                if edge not in remaining:
                    raise ArithmeticError('nonclosed half-edge orbit')
                remaining.remove(edge)
                face.append(edge[0])
                edge = successor[edge]
                if edge == start:
                    break
            if len(face)<3 or len(face)!=len(set(face)):
                raise ArithmeticError('degenerate or nonsimple face')
            j = face.index(min(face))
            faces.append(tuple(face[j:]+face[:j]))
        if len(points)-len(edges)+len(faces) != 2:
            raise ArithmeticError('Euler characteristic failed')
        # Every exact point must participate in this connected planar graph.
        seen, pending = set(), [0]
        while pending:
            u = pending.pop()
            if u not in seen:
                seen.add(u)
                pending.extend(adjacent[u])
        if len(seen)!=len(points):
            raise ArithmeticError('disconnected arrangement')
        faces.sort()
        areas = []
        normalizer = (f.order*(f.root(1)-f.root(-1))).inverse()
        edge_areas = {}
        for u,v in sorted(edges):
            edge_areas[u,v] = wedge(points[u],points[v])*normalizer
        exterior, total = 0, f.zero
        rationals = Counter()
        for face in faces:
            area = f.zero
            for u,v in zip(face,face[1:]+face[:1]):
                area += edge_areas[u,v] if u<v else -edge_areas[v,u]
            if area == -1:
                exterior += 1
                continue
            real, imag = embedding.point(area)
            if real <= tolerance or abs(imag)>tolerance:
                raise ArithmeticError('nonpositive/nonreal bounded face')
            total += area
            areas.append(area.key())
            if area_multiset is not None:
                area_multiset[area.key()] += 1
            q = area.rational()
            if q is not None:
                rationals[str(q)] += 1
        if exterior != 1 or total != 1:
            raise ArithmeticError('exterior face or exact total area failed')
        return {'vertices':len(points),'edges':len(edges),'cells':len(areas),
                'euler_characteristic':2,'exact_area_sum':'1',
                'rational_cells':sum(rationals.values()),
                'rational_area_counts':dict(sorted(rationals.items())),
                'face_cycles_sha256':digest(faces),'cell_areas_sha256':digest(areas)}


def run(ns, precisions=(60,120)):
    results = []
    for n in ns:
        f, points, chords = build_arrangement(n)
        first = None
        for precision in precisions:
            report = enumerate_faces(f,points,chords,precision)
            if first is not None and report != first:
                raise ArithmeticError(f'precision-dependent result for n={n}')
            first = report
        results.append({'n':n,**first})
        print(f'Exploratory cells: n={n}, cells={first["cells"]}, rational={first["rational_cells"]}',
              file=sys.stderr,flush=True)
    return {'method':'exact cyclotomic intersections and areas; numerical planar ordering',
            'certified_arrangement':False,'ordering_precision_digits':list(precisions),
            'precision_results_agree':True,'polygons':results}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,nargs='+',default=list(range(3,27))+[30])
    args = parser.parse_args()
    if any(n<3 for n in args.n) or len(set(args.n)) != len(args.n):
        parser.error('--n requires distinct integers at least 3')
    print(json.dumps(run(args.n),indent=2))
