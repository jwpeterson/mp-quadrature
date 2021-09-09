#!/usr/bin/env python3

import sympy as sy
from fractions import Fraction
import numpy as np

# Compute (and expand) a linear combination from the given
# coefficients and monomials dicts.
def linear_comb(coeffs, monomials):
    comb = 0
    for key, val in coeffs.items():
        comb += val * monomials[key]
    return sy.expand(comb)

# Returns a polynomial string with the given coefficients and terms for printing
def poly_string(coeffs):
    poly_string = ''
    count = 0
    for key, val in coeffs.items():
        # For positive coeffs after the first, print a plus sign
        if count > 0 and val > 0:
            poly_string += ' + '

        # For negative coeffs after the first, we need a space
        if count > 0 and val < 0 and val != -1:
            poly_string += ' '

        # For -1 coefficients print only a minus sign
        if val == -1:
            poly_string += ' - '

        # Print "coeff*" unless coeff is 1 or -1
        if val != 1 and val != -1:
            poly_string += str(val) + '*'

        poly_string += 'Orb(' + key + ')'
        count += 1
    return poly_string

# Returns a string key like 'x2y3' based on the input powers and some
# conventions. We omit terms with 0 exponent, but we explicitly
# include 1 exponents when generating the string.
def key_from_powers(xpower, ypower):
    key = ''
    if xpower > 0:
        key += 'x' + str(xpower)
    if ypower > 0:
        key += 'y' + str(ypower)
    return key

# Check whether the orbit:
# {(a,b), (1-a-b,a), (b, 1-a-b)}
# gives the same quadrature contribution for some different polynomials

a, b = sy.sympify('a, b')
x = [(a,b), (1-a-b,a), (b, 1-a-b)]

# Compute orbits for monomials from degree 2 up to degree dmax. Initialize the
# dict with the constant "1" term. Note that Orb(1) = 3 technically but we
# are free to "normalize" any orbit by multiplying it by a scalar, so we just
# leave it as 1.
dmax = 6
monomial_orbits = {'1':1}
for d in range(2, dmax+1):
    for ypower in range(d+1):
        for i in range(3):
            # The xpower is the complement of the ypower
            xpower = d - ypower

            # Create dict key.
            key = key_from_powers(xpower, ypower)

            # Debugging
            # print(f'key={key}')

            # Compute monomial
            xi = x[i][0]
            yi = x[i][1]

            # KeyError: 'x3'
            # monomial_orbits[key] += (xi**xpower) * (yi**ypower)

            # Doing this in 1 line doesn't work for some reason, I think
            # it is the += operator? The error is:
            # SyntaxError: can't assign to function call
            # monomial_orbits.setdefault(key, 0) += (xi**xpower) * (yi**ypower)

            # So we use two lines
            monomial = monomial_orbits.setdefault(key, 0)
            monomial_orbits[key] += (xi**xpower) * (yi**ypower)

# Debugging
print('monomial_orbits =')
for key, val in monomial_orbits.items():
    print(f'  {key} = {sy.expand(val)}')

print('')
print('Results:')

# Quadratic
# Orb(x**2) + 2*Orb(x*y) - 1 = 0
quadratic_coeffs = {'x2':1, 'x1y1':2, '1':-1}
print('')
print('> Linear combination of quadratic monomial orbits which is not LI:')
print(f'{poly_string(quadratic_coeffs)} = {linear_comb(quadratic_coeffs, monomial_orbits)}')
print('')

# Cubic
# Orb(x**3) + Orb(x**2y) + Orb(xy**2) - Orb(x**2) = 0
cubic_coeffs = {'x3':1, 'x2y1':1, 'x1y2':1, 'x2':-1}
print('> Linear combination of cubic monomial orbits which is not LI:')
print(f'{poly_string(cubic_coeffs)} = {linear_comb(cubic_coeffs, monomial_orbits)}')
print('')

# Quintic
# 4*Orb(x**5) + 10*Orb(x**4 y + x**3 y**2) - 5*Orb(x**4) - 10*Orb(x**2 y) + 1
quintic_coeffs = {'x5':4, 'x4y1':10, 'x3y2':10, 'x4':-5, 'x2y1':-10, '1':1}
print('> Linear combination of quintic monomial orbits which is not LI:')
print(f'{poly_string(quintic_coeffs)} = {linear_comb(quintic_coeffs, monomial_orbits)}')
print('')

# Sixth-order
d = 6

# One more than "a_d" for this d. If we consider a linear combination
# of the first n basis functions, this is sufficient to show linear
# dependence.  For the d=6 case, we only have to consider up to x^3 *
# y^3 to show linear dependence.
n = 4

# Matrix to store coeffs
# We use a sympy Matrix since that has the nullspace() command
A_sym = sy.zeros(n, n)

for ypower in range(n):
    xpower = d - ypower
    key = key_from_powers(xpower, ypower)

    # Debugging
    # print(f'Coeffs of {key}:')

    # Build a Poly object out of the orbit with this key. We can then
    # call coeff_monomial() on the Poly object.
    poly = sy.Poly(monomial_orbits[key])

    # Debugging
    # print(f'poly = {poly}')

    # Extract powers of a^p * b^q in loop
    for bpower in range(n):
        apower = d - bpower
        coeff = poly.coeff_monomial(a**apower * b**bpower)

        # Debugging
        # print(f'coeff of a^{apower} * b^{bpower} = {coeff}')

        # We fill in the matrix of coefficients column by column.
        # This way, the first row of A is all the a^p coeffs, the
        # second row is all the a^{p-1} * b coeffs, etc.
        # Store in array
        A_sym[bpower, ypower] = coeff

# Print result
# print(f'A_sym = {A_sym}')

# Compute nullspace of this matrix. The result is:
# [1/2],
# [3/2],
# [3/2],
# [  1]
# print(f'A_sym.nullspace() = {A_sym.nullspace()}')

# The following degree=six coeffs are the (scaled) entries of the nullspace of the matrix A_sym,
# and we can indeed verify that the 6th-order terms are now cancelled by this linear combination.
# The remaining terms of the linear combination were then found by "inspecting" the remaining
# terms.
sixth_coeffs = {'x6':10, 'x5y1':30, 'x4y2':30, 'x3y3':20, 'x4y1':-30, 'x5':-18, 'x3':20, 'x2':-15, '1':3}
print('> Linear combination of sixth order monomial orbits which is not LI:')
print(f'{poly_string(sixth_coeffs)} = {linear_comb(sixth_coeffs, monomial_orbits)}')
print('')
