#!/usr/bin/env python3

import sympy as sy
from fractions import Fraction

# Function that computes the Reynolds polynomials for a given total
# order d The input x_prime are the transformed reference coordinates
# gi * x, where gi, i=1..3 are the 2x2 rotation matrices associated
# with the Ro3 group.
def compute_reynolds(d, x_prime):
    polys = []
    Reynolds_polys = []
    for xpower in range(d+1):
        # Consider polynomial: x**xpower * y**ypower
        ypower = d - xpower
        polys.append(x**xpower * y**ypower)
        # Compute Reynolds operator for current polynomial
        Reynolds = 0
        for i in range(3):
            xi = x_prime[i][0]
            yi = x_prime[i][1]
            # Note that the (1/3) is part of the standard definition
            # of the Reynolds operator, but we consider any scalar
            # multiple of a given basis function to be equivalent.
            Reynolds += Fraction(1,3) * xi**xpower * yi**ypower
        Reynolds_polys.append(sy.simplify(Reynolds))

    return polys, Reynolds_polys

# Converts input vector of Reynolds polynomials defined in (x,y)
# coordinates to polar coordinates, and returns the results.
def to_polar(Reynolds_polys):
    Reynolds_polar = []
    r = sy.sympify('r')
    for q in Reynolds_polys:
        q_polar = q.subs([(x, r*sy.cos(theta)), (y, r*sy.sin(theta))])
        Reynolds_polar.append(sy.simplify(q_polar))
    return Reynolds_polar

####################
# Code starts here
####################

# This script computes the rotation matrices for the D3 region in R2,
# then uses them in conjunction with the "Reynolds operator" to
# compute the Ro3-invariant polynomials of a given degree.
theta = sy.sympify('theta')

# 2D rotation by angle theta about the x-axis
rot = sy.Matrix([[sy.cos(theta), -sy.sin(theta)],
                 [sy.sin(theta), sy.cos(theta)]])

# The group Ro3 consists of rotations by 2*pi/3, 4*pi/3, and 2*pi (identity)
# g is a list of sympy Matrices
g = [rot.subs(theta, Fraction(2,3)*sy.pi),
     rot.subs(theta, Fraction(4,3)*sy.pi),
     rot.subs(theta, 2*sy.pi)]

# Compute x_prime_i = g_i * x
# using symbolic matrix-vector multiplication
x, y = sy.sympify('x, y')
vec_x = sy.Matrix([[x], [y]])
x_prime = [g[i] * vec_x for i in range(3)]

# Compute Reynolds operator for linear basis functions.
# Note: there are no linear Ro3-invariant polynomials!
linear, Reynolds_linear = compute_reynolds(1, x_prime)
print(f'linear = {linear}')
print(f'Reynolds_linear = {Reynolds_linear}')

print('')

# Compute Reynolds operator for quadratic basis functions.
# Note: there is only 1 unique quadratic Ro3-invariant polynomial:
# P2 := r^2
quadratic, Reynolds_quadratic = compute_reynolds(2, x_prime)
Reynolds_quadratic_polar = to_polar(Reynolds_quadratic)
print(f'quadratic = {quadratic}')
print(f'Reynolds_quadratic = {Reynolds_quadratic}')
print(f'Reynolds_quadratic_polar = {Reynolds_quadratic_polar}')

print('')

# Compute Reynolds operator for cubic basis functions.
# Note: there are 2 unique quadratic Ro3-invariant polynomials:
# P3a := r^3 * cos(theta)
# P3b := r^3 * sin(theta)
cubic, Reynolds_cubic = compute_reynolds(3, x_prime)
Reynolds_cubic_polar = to_polar(Reynolds_cubic)
print(f'cubic = {cubic}')
print(f'Reynolds_cubic = {Reynolds_cubic}')
print(f'Reynolds_cubic_polar = {Reynolds_cubic_polar}')

print('')

# Compute Reynolds operator for quartic basis functions.
# Note: there is 1 unique quartic Ro3-invariant polynomial,
# P4 := P2**2 == r^4
# This can be seen as a consequence of the fact that the only way to
# get 4 by raising 2 and 3 to powers is by squaring 2.
quartic, Reynolds_quartic = compute_reynolds(4, x_prime)
Reynolds_quartic_polar = to_polar(Reynolds_quartic)
print(f'quartic = {quartic}')
print(f'Reynolds_quartic = {Reynolds_quartic}')
print(f'Reynolds_quartic_polar = {Reynolds_quartic_polar}')

print('')

# At 6th-order, the basis functions are computed by taking P2^r * P3a^s * P3b^t, where:
# r & s & t
# ---------
# 3 & 0 & 0
# 0 & 2 & 0
# 0 & 0 & 2
# 0 & 1 & 1
#
# However, a_d==3 for d==6, so we know that only three of these four
# 6th-order polynomials are unique. Using double-angle trig formulas,
# we can show that:
#
# P2^3 = r^6
# P3a^2 = [r^3 cos(3*theta)]^2 = r^6 cos^2(3*theta) = r^6 (1 + cos(6*theta))/2
# P3b^2 = [r^3 sin(3*theta)]^2 = r^6 sin^2(3*theta) = r^6 (1 - cos(6*theta))/2
# P3a * P3b = [r^3 cos(3*theta)] * [r^3 sin(3*theta)] = r^6 sin(6*theta)/2
#
# Thus, by selecting the basis:
#
# r^6 * {1, cos(6*theta), sin(6*theta)}
#
# we can represent any 6th-order Ro3-invariant polynomial.
sixth, Reynolds_sixth = compute_reynolds(6, x_prime)
Reynolds_sixth_polar = to_polar(Reynolds_sixth)
print(f'sixth = {sixth}')
print(f'Reynolds_sixth = {Reynolds_sixth}')
print(f'Reynolds_sixth_polar = {Reynolds_sixth_polar}')

print('')

# 10th-order.
# These are formed by taking P2^r * P3a^s * P3b^t, where:
# r & s & t
# ---------
# 5 & 0 & 0
# 2 & 2 & 0
# 2 & 0 & 2
# 2 & 1 & 1
# So the situation here is very similar to the 6th-order case, the only difference is
# the r column. Therefore we can select the basis:
#
# r^10 * {1, cos(6*theta), sin(6*theta)}
#
# to represent any 10th-order Ro3-invariant polynomial
tenth, Reynolds_tenth = compute_reynolds(10, x_prime)
Reynolds_tenth_polar = to_polar(Reynolds_tenth)
print(f'tenth = {tenth}')
print(f'Reynolds_tenth = {Reynolds_tenth}')
print(f'Reynolds_tenth_polar = {Reynolds_tenth_polar}')
