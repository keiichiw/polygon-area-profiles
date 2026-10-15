#!/usr/bin/env python3
"""60-digit numerical consistency checks, not a proof of rationality/irrationality.

Independent geometry: Cartesian line intersection and polygon shoelace areas.
No third-party packages are needed. All trigonometry uses Decimal Taylor series.
"""
import argparse
import json
from decimal import Decimal as D, getcontext, localcontext


def pi():
    # Quadratically convergent Gauss--Legendre iteration, with guard digits.
    with localcontext() as context:
        context.prec += 12
        a,b,t,p = D(1),1/D(2).sqrt(),D(1)/4,D(1)
        for _ in range(8):
            next_a = (a+b)/2
            b = (a*b).sqrt()
            t -= p*(a-next_a)**2
            a,p = next_a,2*p
        value = (a+b)**2/(4*t)
    return +value


def sin_cos(x, circle_pi):
    x %= 2*circle_pi
    if x > circle_pi:
        x -= 2*circle_pi
    if x < -circle_pi:
        x += 2*circle_pi
    square = x*x
    sine, cosine, st, ct = x,D(1),x,D(1)
    for j in range(1,1000):
        st *= -square/D((2*j)*(2*j+1))
        ct *= -square/D((2*j-1)*(2*j))
        next_s, next_c = sine+st, cosine+ct
        if next_s == sine and next_c == cosine:
            return sine,cosine
        sine,cosine = next_s,next_c
    raise RuntimeError('Taylor series did not converge')


def cross(a,b):
    return a[0]*b[1] - a[1]*b[0]


def subtract(a,b):
    return a[0]-b[0],a[1]-b[1]


def intersection(a,b,c,d):
    first, second = subtract(b,a),subtract(d,c)
    factor = cross(subtract(c,a),second)/cross(first,second)
    return a[0]+factor*first[0],a[1]+factor*first[1]


def expected_values(n,k,m):
    r0 = {(6,1,2):D(1)/18, (18,3,5):D(1)/18,
          (18,4,6):D(1)/9, (30,3,10):D(1)/30, (30,5,12):D(1)/10}
    r2 = {(6,1,2):D(13)/18, (9,1,3):D(2)/3}
    r1 = {(6,1,2):D(1)/9, (12,1,4):D(1)/6}
    config = n,k,m
    a0,a2,a1 = r0.get(config),r2.get(config),r1.get(config)
    if n%2 == 0 and k == n//2-3 and m == n//2-2:
        a2 = D(1)/2+D(3)/n
    adjacent = None
    if n%2 == 0 and m == n//2-1:
        adjacent = D(1)/2-D(2)/n
    elif n == 12 and m == 3:
        adjacent = D(1)/12
    return {'a0':a0, 'a2':a2, 'a1':a1, 'adjacent':adjacent}


def run(n_max, precision, tolerance):
    getcontext().prec = precision
    circle_pi = pi()
    max_errors = {key:D(0) for key in
                  ('area_formulas','sum_one','lateral_symmetry','T_parameterizations',
                   'classified_values','E2_reflection')}
    nonmember_distances = {key:D('Infinity') for key in ('R','S','T','U')}
    counts = {'configurations':0, 'a0':0, 'a2':0, 'a1':0, 'adjacent':0}
    def record(key,value):
        max_errors[key] = max(max_errors[key],abs(value))
        assert abs(value) < tolerance, (key,str(value))
    for n in range(5,n_max+1):
        sine,cosine = {},{}
        def sc(j):
            if j not in sine:
                sine[j],cosine[j] = sin_cos(D(j)*circle_pi/(2*n),circle_pi)
            return sine[j],cosine[j]
        def sn(j):
            return sc(j)[0]
        def cs(j):
            return sc(j)[1]
        vertices = [(cs(4*j),sn(4*j)) for j in range(n)]
        vertices.append(vertices[0])
        edges = [D(0)]
        for j in range(n):
            edges.append(edges[-1]+cross(vertices[j],vertices[j+1]))
        polygon_area = edges[-1]/2
        for m in range(2,(n-1)//2+1):
            for k in range(1,m):
                counts['configurations'] += 1
                l,t = n-2*m,m-k
                end = n-m+k
                x = intersection(vertices[0],vertices[m],vertices[k],vertices[end])
                # Four independently measured boundary polygons, in cyclic order.
                ratios = []
                for start,stop in ((0,k),(k,m),(m,end),(end,n)):
                    doubled = (edges[stop]-edges[start]+cross(vertices[stop],x)
                               +cross(x,vertices[start]))
                    ratios.append(abs(doubled)/(2*polygon_area))
                a0,a1,a2,a3 = ratios
                record('sum_one',sum(ratios)-1)
                record('lateral_symmetry',a1-a3)
                r = sn(2*k)*sn(l)/(sn(4)*sn(2*k+l))
                s = sn(l)*cs(2*t-l)/(sn(4)*cs(2*t))
                lateral_t = 2*sn(l)**2*sn(2*t)/(cs(2*t)*sn(4))
                big_c = n-2*t
                lateral_t_bc = 2*sn(l)**2*cs(big_c)/(sn(4)*sn(big_c))
                adjacent_u = sn(4*m)/sn(4)
                record('T_parameterizations',lateral_t_bc-lateral_t)
                formulas = (D(k)/n-2*r/n, D(t)/n-lateral_t/n,
                            D(k+l)/n+2*s/n, D(t)/n-lateral_t/n,
                            D(m)/n-adjacent_u/n, D(big_c)/n+2*lateral_t_bc/n)
                measured = (a0,a1,a2,a3,a0+a1,a0+a2)
                for actual,formula in zip(measured,formulas):
                    record('area_formulas',actual-formula)
                expected = expected_values(n,k,m)
                for key,actual in (('a0',a0),('a2',a2),('a1',a1),('adjacent',a0+a1)):
                    if expected[key] is not None:
                        counts[key] += 1
                        record('classified_values',actual-expected[key])
                # Recover arithmetic quantities from the independently measured
                # areas, and check consistency with the theorem's finite value sets.
                recovered = ((D(k)-n*a0)/2, (n*a2-k-l)/2,
                             D(t)-n*a1, D(m)-n*(a0+a1))
                for key,type_key,value,allowed in zip(
                        ('R','S','T','U'),('a0','a2','a1','adjacent'),recovered,
                        ((D(1)/3,D(1)),(D(2)/3,D(1)),(D(1)/3,D(1)),(D(1),D(2)))):
                    distance = min(abs(value-q) for q in allowed)
                    if expected[type_key] is None:
                        nonmember_distances[key] = min(nonmember_distances[key],distance)
                        assert distance > tolerance, (n,k,m,key,str(distance))
                    else:
                        assert distance < tolerance, (n,k,m,key,str(distance))
                if l == 4 and t == 1:
                    # Reflection of O in chord V_m V_end is V_m+V_end.
                    reflection = tuple(vertices[m][j]+vertices[end][j] for j in (0,1))
                    for j in (0,1):
                        record('E2_reflection',x[j]-reflection[j])
    return {'status':'PASS','decimal_precision':precision,'n_range':[5,n_max],
            'absolute_tolerance':str(tolerance),'counts':counts,
            'maximum_absolute_errors':{k:str(v) for k,v in max_errors.items()},
            'minimum_nonmember_distance_from_theorem_values':
                {k:str(v) for k,v in nonmember_distances.items()},
            'scope':'Numerical consistency only; no irrationality claim or proof by numerical search.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n-max',type=int,default=80)
    parser.add_argument('--precision',type=int,default=60)
    parser.add_argument('--tolerance',default='1e-50')
    args = parser.parse_args()
    assert args.n_max >= 5 and args.precision >= 60
    print(json.dumps(run(args.n_max,args.precision,D(args.tolerance)),indent=2))


if __name__ == '__main__':
    main()
