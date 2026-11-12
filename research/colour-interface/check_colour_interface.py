#!/usr/bin/env python3
"""Algebra checks for the colour-interface research note.

Python standard library only. Matrix checks use exact Fraction arithmetic.
The effect-body checks use the unitarily equivalent real Pauli pair
(sigma_x, sigma_z), rather than (R_A, Q_A). The active-carrier checks use
the original sector grading diag(I, -I) and its real reflection operator.
The one arbitrary-phase check uses complex floating-point arithmetic, with
an explicit 1e-12 tolerance. No empirical colour data or calibration is used.

These finite checks are regression examples; the note supplies the general
proofs. In particular, they do not identify colour with intersector phase,
prove physiological access, or interpret the final defect.
"""

from fractions import Fraction as F
import cmath
import itertools
import math
import unittest


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(s, a):
    return [[s * x for x in row] for row in a]


def multiply(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return [list(row) for row in zip(*a)]


def blocks(a, b, c, d):
    return [ar + br for ar, br in zip(a, b)] + [cr + dr for cr, dr in zip(c, d)]


def determinant2(a):
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def positive2(a):
    """Independent PSD test by all principal minors of a real symmetric 2x2."""
    if a != transpose(a):
        return False
    return a[0][0] >= 0 and a[1][1] >= 0 and determinant2(a) >= 0


def is_effect2(a):
    return positive2(a) and positive2(add(eye(2), scale(-1, a)))


def effect(x, y, z):
    """E = ((1+z)I + x sigma_x + y sigma_z)/2 in a real basis."""
    x, y, z = map(F, (x, y, z))
    return [[(1 + z + y) / 2, x / 2],
            [x / 2, (1 + z - y) / 2]]


def graph_compression(a):
    """W* a W for W=(I,I)^T/sqrt(2), without irrational arithmetic."""
    v = eye(2) + eye(2)
    return scale(F(1, 2), multiply(transpose(v), multiply(a, v)))


class ColourInterfaceChecks(unittest.TestCase):
    def test_effect_body_against_principal_minors(self):
        grid = [F(k, 4) for k in range(-6, 7)]
        for x, y, z in itertools.product(grid, repeat=3):
            bicone = abs(z) <= 1 and x*x + y*y <= (1 - abs(z))**2
            self.assertEqual(is_effect2(effect(x, y, z)), bicone, (x, y, z))

    def test_illustrative_points_have_reported_eigenvalues(self):
        examples = [(F(1, 2), -1, (F(-1, 4), F(1, 4))),
                    (F(1, 2), 1, (F(3, 4), F(5, 4))),
                    (F(13, 10), 0, (F(-3, 20), F(23, 20)))]
        minus, plus = [[1], [-1]], [[1], [1]]
        for r, z, (lo, hi) in examples:
            e = effect(r, 0, z)
            self.assertEqual(multiply(e, minus), scale(lo, minus))
            self.assertEqual(multiply(e, plus), scale(hi, plus))
            self.assertFalse(is_effect2(e))

    def test_minimal_uniform_enlargement_and_reference_points(self):
        for r, z, expected in [(F(1, 2), -1, F(3, 2)),
                               (F(1, 2), 1, F(3, 2)),
                               (F(13, 10), 0, F(13, 10))]:
            s = abs(z) + r
            self.assertEqual(s, expected)
            self.assertTrue(is_effect2(effect(r/s, 0, F(z)/s)))
            smaller = s - F(1, 100)
            self.assertFalse(is_effect2(effect(r/smaller, 0, F(z)/smaller)))
        self.assertEqual(effect(0, 0, F(-2, 3)), scale(F(1, 6), eye(2)))
        self.assertEqual(effect(0, 0, F(2, 3)), scale(F(5, 6), eye(2)))

    def test_pure_vector_interval_identity_exact_complex_components(self):
        # a=ar+i ai, b=br+i bi; all four real components are rational.
        for ar, ai, br, bi in itertools.product([F(-1), F(1, 3), F(2)], repeat=4):
            aa, bb = ar*ar + ai*ai, br*br + bi*bi
            x, y = 2*(ar*br + ai*bi), 2*(ar*bi - ai*br)
            z, t = aa + bb, aa - bb
            self.assertEqual(z*z - x*x - y*y, t*t)

    def test_mixed_interval_has_determinant_term(self):
        # M=[[a,c-id],[c+id,b]], with arbitrary rational Hermitian entries.
        for a, b, c, d in itertools.product([F(0), F(1, 2), F(1)], repeat=4):
            z, t, x, y = a+b, a-b, 2*c, 2*d
            det = a*b - c*c - d*d
            self.assertEqual(z*z - x*x - y*y, t*t + 4*det)
        # Equal incoherent sector weights are not on the pure-state null circle.
        a = b = F(1, 2)
        self.assertEqual((a-b)**2, 0)
        self.assertEqual((a+b)**2, 1)

    def test_bicone_conjugation_and_complement_symmetries(self):
        # In this real basis the omitted grading is sigma_y = i K.
        k = [[F(0), F(-1)], [F(1), F(0)]]
        e = effect(F(1, 4), F(1, 3), F(1, 5))
        jej = scale(-1, multiply(k, multiply(e, k)))
        self.assertEqual(jej, effect(F(-1, 4), F(-1, 3), F(1, 5)))
        complement = add(eye(2), scale(-1, e))
        self.assertEqual(complement, effect(F(-1, 4), F(-1, 3), F(-1, 5)))
        self.assertEqual(add(eye(2), scale(-1, jej)),
                         effect(F(1, 4), F(1, 3), F(-1, 5)))

    def test_active_graph_compression_is_scalar(self):
        i, zero = eye(2), scale(0, eye(2))
        reflection = blocks(zero, i, i, zero)
        for x, z in itertools.product([F(-1, 2), F(0), F(1, 3)], repeat=2):
            e = scale(F(1, 2), add(scale(1+z, eye(4)), scale(x, reflection)))
            self.assertEqual(graph_compression(e), scale((1+z+x)/2, i))

    def test_arbitrary_phase_statistics_are_input_independent(self):
        # Floating-point check of the original complex Q_A convention.
        tol = 1e-12
        x, y, z = .2, .3, -.1
        for theta in [0., .37, math.pi/2, math.pi]:
            expected = (1+z+x*math.cos(theta)+y*math.sin(theta))/2
            for p in [(1+0j, 0j), (.3+.4j, .5-.2j), (1j, 2-3j)]:
                norm = math.sqrt(sum(abs(v)**2 for v in p))
                p = [v/norm for v in p]
                v = [a/math.sqrt(2) for a in p] + [cmath.exp(1j*theta)*a/math.sqrt(2) for a in p]
                ev = [(1+z)*a/2 + (x-1j*y)*b/2 for a, b in zip(v[:2], v[2:])]
                ev += [(x+1j*y)*a/2 + (1+z)*b/2 for a, b in zip(v[:2], v[2:])]
                actual = sum(a.conjugate()*b for a, b in zip(v, ev))
                self.assertLess(abs(actual-expected), tol)

    def test_sector_blind_condition_allows_more_than_scalar_family(self):
        i, zero = eye(2), scale(0, eye(2))
        s = [[F(1), F(0)], [F(0), F(-1)]]
        j = blocks(i, zero, zero, scale(-1, i))
        e = scale(F(1, 2), blocks(i, s, s, i))
        self.assertEqual(multiply(e, e), e)  # Orthogonal projection, hence effect.
        self.assertEqual(e, transpose(e))
        self.assertEqual(add(e, multiply(j, multiply(e, j))), eye(4))
        # Every scalar structural-Pauli effect has scalar off-diagonal blocks.
        self.assertNotEqual(e[0][2], e[1][3])

    def test_defect_extension_invisible_statically_visible_after_mixing(self):
        # Physical dimension 2, phenomenal dimension 3, u is first-two inclusion.
        f = [[F(1, 4), F(0)], [F(0), F(3, 4)]]
        def extension(t):
            return [[f[0][0], 0, 0], [0, f[1][1], 0], [0, 0, F(t)]]
        e0, e1 = extension(0), extension(1)
        for psi in [[[1], [0], [0]], [[0], [1], [0]], [[F(3,5)], [F(4,5)], [0]]]:
            a = multiply(transpose(psi), multiply(e0, psi))
            b = multiply(transpose(psi), multiply(e1, psi))
            self.assertEqual(a, b)
        swap = [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
        psi = multiply(swap, [[1], [0], [0]])
        self.assertEqual(multiply(transpose(psi), multiply(e0, psi)), [[0]])
        self.assertEqual(multiply(transpose(psi), multiply(e1, psi)), [[1]])

    def test_internal_paired_readout_varies_with_input(self):
        f = [[F(1), F(0)], [F(0), F(0)]]
        zero = scale(0, eye(2))
        paired = blocks(f, zero, zero, f)
        self.assertEqual(graph_compression(paired), f)
        self.assertEqual(f[0][0], 1)
        self.assertEqual(f[1][1], 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
