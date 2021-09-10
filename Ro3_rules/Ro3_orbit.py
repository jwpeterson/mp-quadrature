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

# Returns the number of basis functions for a given degree, d.
# These numbers are given by the Molien series for the Ro3 basis.
# Note that the _total_ size of the basis for a given d is the
# sum_d a_d.
def a_d(d):
    #    a_0, a_1, a_2, a_3, a_4
    a = [  1,   0,   1,   2,   1]

    # Handle pre-computed cases
    if d < len(a):
        return a[d]

    # Use formula to compute a_d for larger d
    a_ell = 0
    for ell in range(5, d+1):
        a.append(a[ell-1] + a[ell-3] - a[ell-4])

    # Debugging
    # print(f'a = {a}')

    # The last entry in a is the one we requested
    return a[-1]



# Check whether the orbit:
# {(a,b), (1-a-b,a), (b, 1-a-b)}
# gives the same quadrature contribution for some different polynomials

a, b = sy.sympify('a, b')
x = [(a,b), (1-a-b,a), (b, 1-a-b)]

# Compute orbits for monomials from degree 2 up to degree dmax. Initialize the
# dict with the constant "1" term. Note that Orb(1) = 3 technically but we
# are free to "normalize" any orbit by multiplying it by a scalar, so we just
# leave it as 1.
dmax = 7
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

# Quadratic (a_2 = 1)
# Orb(x**2) + 2*Orb(x*y) - 1 = 0
quadratic_coeffs = {'x2':1, 'x1y1':2, '1':-1}
print('')
print('> Linear combination of quadratic monomial orbits which is not LI:')
print(f'{poly_string(quadratic_coeffs)} = {linear_comb(quadratic_coeffs, monomial_orbits)}')
print('')

# Cubic (a_3 = 2)
# Orb(x**3) + Orb(x**2y) + Orb(xy**2) - Orb(x**2) = 0
cubic_coeffs = {'x3':1, 'x2y1':1, 'x1y2':1, 'x2':-1}
print('> Linear combination of cubic monomial orbits which is not LI:')
print(f'{poly_string(cubic_coeffs)} = {linear_comb(cubic_coeffs, monomial_orbits)}')
print('')

# Quintic (a_5 = 2)
# 4*Orb(x**5) + 10*Orb(x**4 y + x**3 y**2) - 5*Orb(x**4) - 10*Orb(x**2 y) + 1
quintic_coeffs = {'x5':4, 'x4y1':10, 'x3y2':10, 'x4':-5, 'x2y1':-10, '1':1}
print('> Linear combination of quintic monomial orbits which is not LI:')
print(f'{poly_string(quintic_coeffs)} = {linear_comb(quintic_coeffs, monomial_orbits)}')
print('')

# Sixth-order (a_6 = 3)
# The following degree=six coeffs are the (scaled) entries of the nullspace of the matrix A_sym,
# and we can indeed verify that the 6th-order terms are now cancelled by this linear combination.
# The remaining terms of the linear combination were then found by "inspecting" the remaining
# terms.
sixth_coeffs = {'x6':10, 'x5y1':30, 'x4y2':30, 'x3y3':20, 'x4y1':-30, 'x5':-18, 'x3':20, 'x2':-15, '1':3}
print('> Linear combination of sixth order monomial orbits which is not LI:')
print(f'{poly_string(sixth_coeffs)} = {linear_comb(sixth_coeffs, monomial_orbits)}')
print('')

# Seventh-order (a_7 = 2)
# Linear system of equations to cancel 6th-order terms
# A = sy.Matrix([[2, -1, 1, -7], [6, -4, 4, -14], [15, -10, 7, 0]])
# print(f'A={A}')
# print(f'A.rref()={A.rref()}')
# Result
# [1, 0, 0,    -7],
# [0, 1, 0, -56/3],
# [0, 0, 1, -35/3]]), (0, 1, 2))
# seventh_coeffs = {'x7':2, 'x6y1':7, 'x5y2':7}
seventh_coeffs = {'x7':18, 'x6y1':63, 'x5y2':63, 'x6':-63, 'x5y1':-168, 'x4y2':-105, 'x5':63, 'x4y1':105, 'x3':-35, 'x2':21, '1':-4}
print('> Linear combination of seventh order monomial orbits which is not LI:')
print(f'{poly_string(seventh_coeffs)} = {linear_comb(seventh_coeffs, monomial_orbits)}')
print('')



# Compute coefficients of linear combination that sum to zero.
d = 7

# Number of basis functions to consider for this d.
# .) If n = a_d(d) + 1, then we should get a non-trivial nullspace
#    for the linear system of equations that we build.
# .) If n = a_d(d), we should get only the trivial nullspace.
# For example, in the d=6 case, if we consider n=4 basis functions (up
# to x^3 * y^3) then we should get linear dependence.
n = a_d(d)
# n = a_d(d) + 1

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

# Returns a list of sympy matrices (vectors) representing the nullspace
# Compute nullspace of this matrix.
# For d=6, with n=4, for example, the result is:
# [1/2],
# [3/2],
# [3/2],
# [  1]
nullsp = A_sym.nullspace()

# from sympy import shape # does not work
# print(f'Shape of nullspace = {nullsp.shape}') # no
if not nullsp:
    print(f'd = {d} basis with {n} monomials is linearly independent.')
else:
    print(f'd = {d} basis with {n} monomials is linearly dependent with nullspace =\n{nullsp}')

# Test a_d() function
# print(f'a_0 = {a_d(0)} (should be 1)')
# print(f'a_6 = {a_d(6)} (should be 3)')
# print(f'a_7 = {a_d(7)} (should be 2)')
# print(f'a_12 = {a_d(12)} (should be 5)')
# print(f'a_17 = {a_d(17)} (should be 6)')
# print(f'a_24 = {a_d(24)} (should be 9)')
