#!/usr/bin/env python3
"""Exact certificate for a non-null boundary point of a continuum color model.

Input decimal strings are first read as IEEE binary64 values, then interpreted
as exact rationals. PCHIP slopes and Hermite coefficients are computed in exact
rational arithmetic, not by rounding a floating-point interpolator. White
normalizers are outward rational enclosures of exact polynomial integrals.
All sign, root-count, and independence assertions below use interval bounds;
floating-point displays are not used to establish any assertion.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys


@dataclass(frozen=True)
class IV:
    lo: F
    hi: F

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError('Reversed interval')

    @staticmethod
    def point(x):
        x = F(x)
        return IV(x, x)

    def __add__(self, other):
        other = other if isinstance(other, IV) else IV.point(other)
        return IV(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return IV(-self.hi, -self.lo)

    def __sub__(self, other):
        other = other if isinstance(other, IV) else IV.point(other)
        return self + (-other)

    def __mul__(self, other):
        other = other if isinstance(other, IV) else IV.point(other)
        p = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return IV(min(p), max(p))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, IV) else IV.point(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError('Interval contains zero')
        return self * IV(1/other.hi, 1/other.lo)

    def sign(self):
        return 1 if self.lo > 0 else -1 if self.hi < 0 else 0


def sign(x):
    return (x > 0) - (x < 0)


def pchip(y):
    """Exact uniform-grid PCHIP, power coefficients in t=lambda-left."""
    d = [b-a for a, b in zip(y, y[1:])]
    m = [F(0)]*len(y)
    for k in range(1, len(y)-1):
        if d[k-1] and d[k] and sign(d[k-1]) == sign(d[k]):
            m[k] = 2*d[k-1]*d[k]/(d[k-1]+d[k])
    for slot, d0, d1 in [(0, d[0], d[1]), (-1, d[-1], d[-2])]:
        v = (3*d0-d1)/2
        if sign(v) != sign(d0):
            v = F(0)
        elif sign(d0) != sign(d1) and abs(v) > 3*abs(d0):
            v = 3*d0
        m[slot] = v
    return [[y[k], m[k], 3*d[k]-2*m[k]-m[k+1],
             -2*d[k]+m[k]+m[k+1]] for k in range(len(d))]


def outward(x, scale=10**40):
    floor = (x.numerator*scale)//x.denominator
    ceil = -((-x.numerator*scale)//x.denominator)
    return IV(F(floor, scale), F(ceil, scale))


def evaluate(a, t):
    t = t if isinstance(t, IV) else IV.point(t)
    v = a[-1]
    for c in a[-2::-1]:
        v = v*t+c
    return v


def bernstein(a):
    n = len(a)-1
    return [sum((a[k]*F(math.comb(j,k),math.comb(n,k))
                 for k in range(j+1)), IV.point(0)) for j in range(n+1)]


def certify_sign(b, wanted, depth=0):
    """Bernstein convex-hull bound, with exact midpoint subdivision."""
    if all(v.sign() == wanted for v in b):
        return 1
    if depth >= 30:
        raise AssertionError('Bernstein sign not certified after subdivision')
    rows = [b]
    while len(rows[-1]) > 1:
        last = rows[-1]
        rows.append([(last[i]+last[i+1])/2 for i in range(len(last)-1)])
    left = [row[0] for row in rows]
    right = [row[-1] for row in rows[::-1]]
    return certify_sign(left,wanted,depth+1)+certify_sign(right,wanted,depth+1)


def decimal_bound(x, places, upper=False):
    scale = 10**places
    n = -((-x.numerator*scale)//x.denominator) if upper else (x.numerator*scale)//x.denominator
    prefix = '-' if n < 0 else ''
    n = abs(n)
    return f'{prefix}{n//scale}.{n%scale:0{places}d}'


def describe(v, places=14):
    return {'lower': decimal_bound(v.lo, places),
            'upper': decimal_bound(v.hi, places, True),
            'lower_rational': str(v.lo), 'upper_rational': str(v.hi)}


def det3(v):
    a,b,c = v
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            -a[1]*(b[0]*c[2]-b[2]*c[0])
            +a[2]*(b[0]*c[1]-b[1]*c[0]))


def run(path):
    payload = path.read_bytes()
    rows = list(csv.reader(payload.decode().splitlines()))
    if not rows:
        raise ValueError('Empty data')
    try:
        float(rows[0][0])
    except ValueError:
        rows = rows[1:]
    data=[]
    for row in rows:
        if not row:
            continue
        if len(row) != 4:
            raise ValueError('Expected wavelength,L,M,S columns')
        values=[float(s) for s in row]
        if not all(math.isfinite(v) for v in values):
            raise ValueError('Nonfinite input')
        data.append([F.from_float(v) for v in values])
    if len(data) != 441 or [r[0] for r in data] != list(map(F,range(390,831))):
        raise ValueError('Certificate expects the complete 390--830 nm unit grid')
    if any(v < 0 for row in data for v in row[1:]):
        raise ValueError('Negative cone sensitivity')
    polys=[pchip([row[j+1] for row in data]) for j in range(3)]
    whites=[]
    for channel in polys:
        pieces=[outward(sum((a[k]/F(k+1) for k in range(4)),F(0))) for a in channel]
        whites.append(sum(pieces,IV.point(0)))
    if not all(w.lo > 0 for w in whites):
        raise AssertionError('White normalizer not strictly positive')
    balanced=[[[IV.point(a)/whites[j] for a in segment] for segment in polys[j]] for j in range(3)]
    # q=-c_L+2c_M-c_S; omitting positive factor 2/9 changes no signs.
    costs=[[(-balanced[0][i][k]+2*balanced[1][i][k]-balanced[2][i][k])
            for k in range(4)] for i in range(440)]
    transitions=[]
    local_roots=[]
    sign_counts={'positive_segments':0,'negative_segments':0,'transition_segments':0,
                 'bernstein_leaves':0}
    for i,a in enumerate(costs):
        s0,s1=evaluate(a,F(0)).sign(),evaluate(a,F(1)).sign()
        if not s0 or not s1:
            raise AssertionError(f'Unresolved endpoint sign in cell {390+i}')
        if s0 == s1:
            sign_counts['positive_segments' if s0 > 0 else 'negative_segments'] += 1
            sign_counts['bernstein_leaves'] += certify_sign(bernstein(a),s0)
            continue
        derivative=[a[k]*k for k in range(1,4)]
        wanted=s1
        certify_sign(bernstein(derivative),wanted)
        lo,hi=F(0),F(1)
        while hi-lo > F(1,10**10):
            mid=(lo+hi)/2
            sm=evaluate(a,mid).sign()
            if not sm:
                raise AssertionError('Normalizer enclosure too wide to isolate root')
            if sm == s0:
                lo=mid
            else:
                hi=mid
        root=IV(lo,hi)
        local_roots.append((i,root))
        L=evaluate(balanced[0][i],root)
        S=evaluate(balanced[2][i],root)
        x=L-S
        z=L+S
        if not x.sign() or z.lo <= 0:
            raise AssertionError('Crossing-vector independence sign unresolved')
        transitions.append({'cell_nm':[390+i,391+i],
                            'root_nm':describe(root+(390+i)),
                            'derivative_sign':wanted,
                            'cL_minus_cS_at_root':describe(x),
                            'cL_plus_cS_at_root':describe(z)})
        sign_counts['transition_segments']+=1
    if [t['cell_nm'][0] for t in transitions] != [486,599]:
        raise AssertionError('Unexpected transition cells')
    if [t['derivative_sign'] for t in transitions] != [1,-1]:
        raise AssertionError('Unexpected transition orientations')
    # Opposite x signs and positive z imply independent transformed vectors.
    if not (F(transitions[0]['cL_minus_cS_at_root']['upper_rational']) < 0 <
            F(transitions[1]['cL_minus_cS_at_root']['lower_rational'])):
        raise AssertionError('Crossing vectors lack the claimed opposite x signs')
    # The maximizing reflectance is the positive band between the roots.
    # Enclose its catches and coordinates, retaining interval arithmetic at
    # the two uncertain limits and outward rounding each full-cell integral.
    first,r0=local_roots[0]
    last,r1=local_roots[1]
    catches=[]
    for j in range(3):
        def primitive(i):
            return [IV.point(0)]+[a/F(k+1) for k,a in enumerate(balanced[j][i])]
        integral=evaluate(primitive(first),F(1))-evaluate(primitive(first),r0)
        integral+=evaluate(primitive(last),r1)
        for i in range(first+1,last):
            v=evaluate(primitive(i),F(1))
            integral+=IV(outward(v.lo).lo,outward(v.hi).hi)
        catches.append(integral)
    total=sum(catches,IV.point(0))
    xyz=[F(2,3)*(catches[0]-total/3),
         F(2,3)*(catches[1]-total/3),F(2,3)*total-1]
    if xyz[1].lo <= 0 or xyz[2].lo <= 0:
        raise AssertionError('Positive radial support or off-equator height not certified')
    witness_wavelengths=[440,530,610]
    determinant=det3([data[k-390][1:] for k in witness_wavelengths])
    if determinant == 0:
        raise AssertionError('Full-rank witness failed')
    return {'certified':True,'input_sha256':hashlib.sha256(payload).hexdigest(),
            'input_model':'Exact binary64 tabulated values; exact rational unit-grid PCHIP; continuum wavelength integration',
            'domain_nm':[390,830], 'normal_xy':[0,1],
            'cost':'-balanced_L+2*balanced_M-balanced_S',
            'white_integral_bounds':[describe(w,40) for w in whites],
            'sign_certificate':sign_counts,'roots':transitions,
            'crossing_vectors_independent':True,
            'support_point_xyz':[{k:v for k,v in describe(a).items() if not k.endswith('_rational')} for a in xyz],
            'positive_radial_support':True,'strictly_off_equator':True,
            'full_rank_witness':{'wavelengths_nm':witness_wavelengths,
                                 'raw_determinant_rational':str(determinant),
                                 'sign':sign(determinant)},
            'interpretation':'Exactly two simple support transitions, a positive band, and independent crossing vectors. With the accompanying convex-geometric argument these establish a unique horizontal outward normal and non-nullness under a metric with positive chromatic block and no chromatic-height cross terms. This is a continuum model certificate, not empirical psychophysics.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=Path(__file__).parent/'inputs/stockman_sharpe_2deg_390_830.csv')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    try:
        result=run(args.data)
    except Exception as exc:
        print(json.dumps({'certified':False,'error':str(exc)}),file=sys.stderr)
        return 1
    text=json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text,end='')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
