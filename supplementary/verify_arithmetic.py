#!/usr/bin/env python3
"""Regression tests for exact arithmetic and the new geometric calculations."""
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import unittest
import random
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from cyclotomic import Element, Field, profile, wedge
from explore_cells import Embedding, build_arrangement, enumerate_faces
from verify_geometry import pi, sin_cos


class ArithmeticTests(unittest.TestCase):
    def test_field_identities(self):
        for n in range(3,61):
            f = Field(n)
            self.assertEqual(sum(f.powers),0)
            a = 2+f.root(1)-3*f.root(2)
            self.assertEqual(a*a.inverse(),1)
            self.assertEqual(a.conjugate().conjugate(),a)
            for r in (1,n//2,n-1):
                self.assertEqual((1-f.root(r))*f.binomial_inverse(r),1)
            self.assertEqual(Element(f,[6,-3],9),Element(f,[-2,1],-3))
            self.assertEqual((F(7,11)*a).rational_ratio(a),F(7,11))
            self.assertIsNone(f.root(1).rational())
            self.assertIsNone(f.root(1).rational_ratio(f.one))
            with self.assertRaises(ZeroDivisionError):
                f.zero.inverse()
            with self.assertRaises(ZeroDivisionError):
                f.one.rational_ratio(f.zero)
            with self.assertRaises(ZeroDivisionError):
                f.binomial_inverse(0)

    def test_fast_arithmetic_against_naive_polynomials(self):
        rng = random.Random(20261023)
        def reduce_naively(coeff, modulus):
            coeff = list(coeff)
            degree = len(modulus)-1
            for k in range(len(coeff)-1,degree-1,-1):
                c = coeff[k]
                for j in range(degree):
                    coeff[k-degree+j] -= c*modulus[j]
            return tuple((coeff+[0]*degree)[:degree])
        def multiply_naively(a,b):
            coeff = [0]*(len(a)+len(b)-1)
            for i,x in enumerate(a):
                for j,y in enumerate(b):
                    coeff[i+j] += x*y
            return coeff
        for n in (5,8,12,25,30,49,60,64,81,89,95,97,100):
            f = Field(n)
            for magnitude in (3,10**30):
                aa = [rng.randint(-magnitude,magnitude) for _ in range(f.degree)]
                bb = [rng.randint(-magnitude,magnitude) for _ in range(f.degree)]
                a,b = Element(f,aa),Element(f,bb)
                expected = reduce_naively(multiply_naively(aa,bb),f.modulus)
                self.assertEqual((a*b).coeff,expected)
                conjugate = [sum(aa[j]*f.root(-j).coeff[k] for j in range(f.degree))
                             for k in range(f.degree)]
                self.assertEqual(a.conjugate().coeff,tuple(conjugate))
                for j in (-n-1,-1,0,1,n//2,n+1):
                    expected = reduce_naively(multiply_naively(aa,f.root(j).coeff),f.modulus)
                    self.assertEqual(a.rotate(j).coeff,expected)
            long = [rng.randrange(-100,101) for _ in range(3*n+5)]
            self.assertEqual(Element(f,long).coeff,reduce_naively(long,f.modulus))

    def test_dodecagon(self):
        f = Field(12)
        sqrt3 = f.root(1)+f.root(-1)
        self.assertEqual(sqrt3*sqrt3,3)
        normalizer = (12*(f.root(1)-f.root(-1))).inverse()
        for end,den in ((6,18),(8,24)):
            x = f.intersection(0,3,1,end)
            area = wedge(f.root(1)-1,x-1)*normalizer
            self.assertEqual(area,(2*sqrt3-3)*F(1,den))
            self.assertIsNone(area.rational())
        # Enumerate every crossing length-5 chord and check reflection orbits.
        crossing = {(a,b) for a in (1,2) for b in range(4,12)
                    if min((a-b)%12,(b-a)%12)==5}
        reflected = {tuple(sorted(((3-a)%12,(3-b)%12))) for a,b in ((1,6),(1,8))}
        self.assertEqual(crossing,{(1,6),(1,8)}|reflected)

    def test_profiles_against_trigonometric_formulas(self):
        with localcontext() as context:
            context.prec = 60
            p = pi()
            for n in list(range(5,21))+[47,97,100]:
                f = Field(2*n)
                embed = Embedding(f,60)
                def sn(j):
                    return sin_cos(D(j)*p/(2*n),p)[0]
                def cs(j):
                    return sin_cos(D(j)*p/(2*n),p)[1]
                ms = range(2,(n+1)//2) if n<=20 else (2,(n-1)//2)
                for m in ms:
                    ks = range(1,m) if n<=20 else (1,m-1)
                    for k in ks:
                        l,t = n-2*m,m-k
                        r = sn(2*k)*sn(l)/(sn(4)*sn(2*k+l))
                        s = sn(l)*cs(2*t-l)/(sn(4)*cs(2*t))
                        lateral = 2*sn(l)**2*sn(2*t)/(cs(2*t)*sn(4))
                        expected = (D(k)/n-2*r/n,D(t)/n-lateral/n,
                                    D(k+l)/n+2*s/n,D(t)/n-lateral/n)
                        for area,value in zip(profile(n,k,m,f),expected):
                            real,imag = embed.point(area)
                            self.assertLess(abs(real-value),D('1e-50'))
                            self.assertLess(abs(imag),D('1e-50'))
        self.assertEqual([a.rational() for a in profile(6,1,2)],
                         [F(1,18),F(1,9),F(13,18),F(1,9)])

    def test_rotational_quotient_against_full_expansion(self):
        from verify_cell_orbits import compare
        # Odd central faces, even fixed central vertices, concurrency, and the
        # exceptional rational configurations have different stabilizers.
        for n in (3,4,5,6,7,8,9,10,12,15,18):
            compare(n)

    def test_rotation_orbit_canonicalization(self):
        from explore_cell_orbits import Arrangement
        for n in (5,6,7,8,9,12):
            a = Arrangement(n)
            for j in range(n):
                self.assertEqual(a.intern(a.f.root(j)),j)
            z = a.f.intersection(0,2,1,3)
            vertex = a.intern(z)
            for j in range(n):
                self.assertEqual(a.intern(z*a.f.root(j)),a.shift(vertex,j))
            if n%2==0:
                self.assertEqual(a.intern(a.f.zero),n)

    def test_checkpoint_resume_and_stale_code(self):
        import verify_geometry
        root = Path(__file__).resolve().parent
        env = os.environ.copy()
        env['PYTHONPATH'] = str(Path(verify_geometry.__file__).resolve().parent)
        with tempfile.TemporaryDirectory(prefix='cell-resume-test-') as directory:
            directory = Path(directory)
            for name in ('explore_cell_orbits.py','explore_cells.py','cyclotomic.py'):
                shutil.copyfile(root/name,directory/name)
            command = [sys.executable,str(directory/'explore_cell_orbits.py'),
                       '--n','6','--checkpoint-dir',str(directory/'checkpoints')]
            first = subprocess.run(command,env=env,text=True,capture_output=True,check=True)
            second = subprocess.run(command,env=env,text=True,capture_output=True,check=True)
            self.assertEqual(json.loads(first.stdout),json.loads(second.stdout))
            with (directory/'explore_cell_orbits.py').open('a') as stream:
                stream.write('\n')
            stale = subprocess.run(command,env=env,text=True,capture_output=True)
            self.assertNotEqual(stale.returncode,0)
            self.assertIn('Stale checkpoint',stale.stderr)

    def test_cell_topology_and_concurrency(self):
        for n,cells,rational in ((3,1,1),(4,4,4),(5,11,0),(6,24,24),(8,80,0)):
            f,points,chords = build_arrangement(n)
            first = enumerate_faces(f,points,chords,60)
            self.assertEqual(first,enumerate_faces(f,points,chords,120))
            self.assertEqual(first['cells'],cells)
            self.assertEqual(first['rational_cells'],rational)
            self.assertEqual(first['exact_area_sum'],'1')
            if n%2 == 0:
                self.assertEqual(sum(p==0 for p in points),1)


if __name__ == '__main__':
    unittest.main()
