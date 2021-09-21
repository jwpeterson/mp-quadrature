#!/usr/bin/env python3

import sympy as sy
from fractions import Fraction
import numpy as np
from collections import defaultdict

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

# Helper function used by a_d and s_d functions below
def a_d_array(d):
    #    a_0, a_1, a_2, a_3, a_4
    a = [  1,   0,   1,   2,   1]

    # Handle pre-computed cases
    if d < len(a):
        return a[0:d+1]

    # Use formula to compute a_d for larger d
    a_ell = 0
    for ell in range(5, d+1):
        a.append(a[ell-1] + a[ell-3] - a[ell-4])

    return a

# Returns the number of basis functions for a given degree, d.
# These numbers are given by the Molien series for the Ro3 basis.
# Note that the _total_ size of the basis for a given d is the
# sum_d a_d.
def a_d(d):
    # The last entry in the array returned by the a_d_array() function
    # is the one we requested.
    return a_d_array(d)[-1]

# Returns the sum of the a_d up to and including d
def s_d(d):
    a = a_d_array(d)
    return sum(a)

# I found this formula in my notes, it is much simpler to use than
# computing all of the a_d's
def s_d_simple(d):
    return int((d*d + 3*d + 6)/6)

# Given degree d and a the set of monomial_orbits of at least degree
# d, computes the nullspace for the orbits of degree d. If the
# nullspace is non-trivial, this gives us a linear combination of
# degree d polynomials which is linearly-dependent.
def compute_nullspace(d, monomial_orbits):
    # The input monomial_orbits dict contains polynomials a^p * b^q
    a, b = sy.sympify('a, b')

    # Number of basis functions to consider for this d.
    # .) If n = a_d(d) + 1, then we should get a non-trivial nullspace
    #    for the linear system of equations that we build.
    # .) If n = a_d(d), we should get only the trivial nullspace.
    # For example, in the d=6 case, if we consider n=4 basis functions (up
    # to x^3 * y^3) then we should get linear dependence.
    # n = a_d(d)
    n = a_d(d) + 1

    # Matrix to store coeffs
    # We use a sympy Matrix since that has the nullspace() command
    A_sym = sy.zeros(n, n)

    for ypower in range(n):
        xpower = d - ypower
        key = key_from_powers(xpower, ypower)

        # Debugging
        # print(f'monomial_orbits[{key}] = {monomial_orbits[key]}')

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
        # Compute linear combination of polynomials according to the
        # coeffs in nullsp.  We should find that it is at most degree
        # d-1, in other words, that the degree d terms have all
        # cancelled. Notes:
        # 1. We assume that the "nullsp" list has only a single entry,
        #    and we refer directly to it in the loop below.
        # 2. We convert to Poly when building the linear combination, since
        #    I was getting an error for the degree 2 case when trying to
        #    construct the Poly later.

        # Extract nullspace coeffs into list. Currently "nullsp" is a
        # list of Sympy column Matrix objects, which is not as
        # straightforward to work with.
        nullspace_coeffs = [val for val in nullsp[0].col(0)]
        print(f'd = {d} basis with {n} monomials is linearly dependent with nullspace = {nullspace_coeffs}')

        linear_comb = 0
        for ypower in range(n):
            xpower = d - ypower
            key = key_from_powers(xpower, ypower)
            linear_comb += nullspace_coeffs[ypower] * sy.Poly(monomial_orbits[key])

        # Debugging: print result
        # print(f'linear_comb = {linear_comb}')

        # Now check the degree of linear_comb, which should be <= d,
        # confirming that the polynomials tested are not linearly
        # independent.
        linear_comb_degree = sy.degree(linear_comb)
        # print(f'linear_comb_degree = {linear_comb_degree}')
        if not linear_comb_degree < d:
            raise RuntimeError(f'Expected linear combination to have degree < {d}, but got {linear_comb_degree} instead')



# Alternate version of the function above that computes a more
# complete nullspace by finding the coeffs of all monomial orbits, not
# just the highest degree ones.  Must be passed a monomial_orbits
# constructed by calling:
# compute_monomial_orbits(dmax, "nullspace")
# for some d.
def compute_nullspace2(dmax, monomial_orbits):
    # The input monomial_orbits dict contains polynomials a^p * b^q
    a, b = sy.sympify('a, b')

    # We are going to build an nrows x ncols system of equations, where ncols is the
    # number of entries in monomial_orbits, and nrows = (dmax+1)*(dmax+2)/2 is the number of
    # monomial terms in a 2D polynomial of degree dmax.
    nrows = int((dmax+1)*(dmax+2)/2)
    ncols = len(monomial_orbits)

    # Debugging
    # print(f'number of rows = {nrows}')
    # print(f'number of cols = {ncols}')

    # Matrix to store coeffs
    # We use a sympy Matrix since that has the nullspace() command
    A_symb = sy.zeros(nrows, ncols)

    # We are going to fill in column by column
    col = 0
    for key, val in monomial_orbits.items():
        # Debugging
        # print(f'col = {col}, key = {key}')

        # Special case for 1, can't create Poly from 1
        if val == 1:
            if col != 0:
                raise RuntimeError(f'Expected basis function 1 to come first')
            A_symb[0, col] = 1

        else:
            # Debugging
            # print(f'val = {val}')

            # Convert to poly so that we can call coeff_monomial() on it
            poly = sy.Poly(val)

            # Extract powers of a^{apower} * b^{bpower} in loop, where
            # apower + bpower = 0, 1, ..., dmax
            row = 0
            for d in range(dmax+1):
                for bpower in range(d+1):
                    apower = d - bpower

                    # Debugging
                    # print(f'row = {row}, apower = {apower}, bpower = {bpower}')

                    coeff = poly.coeff_monomial(a**apower * b**bpower)

                    # Debugging
                    # print(f'coeff of a^{apower} * b^{bpower} = {coeff}')

                    # Store this coeff in the matrix
                    A_symb[row, col] = coeff

                    # Go to next row
                    row += 1

        # Go to next col
        col += 1

    # Print result
    # print(f'A_symb = {A_symb}')

    # Returns a list of sympy matrices (vectors) representing the nullspace
    # Compute nullspace of this matrix.
    nullsp = A_symb.nullspace()

    # Debugging. For example in the dmax=4 case, we get:
    # [  0],   "1"
    # [1/2],   "x2"
    # [ -1],   "x3"
    # [ -1],   "x2y"
    # [1/2],   "x4"
    # [  1]])] "x3y"
    # print(f'nullsp = {nullsp}')

    if not nullsp:
        print(f'd = {d} basis with {n} monomials is linearly independent.')
    else:
        # Compute linear combination of polynomials according to the
        # coeffs in nullsp.  We should find that it the result is zero!
        # Note: We assume that the "nullsp" list has only a single entry.

        # Extract nullspace coeffs into list. Currently "nullsp" is a
        # list of Sympy column Matrix objects, which is not as
        # straightforward to work with.
        nullspace_coeffs = [val for val in nullsp[0].col(0)]

        # Construct dict for pretty printing and computing linear combinations
        ps = dict(zip(monomial_orbits.keys(), nullspace_coeffs))

        # Debugging
        # print(f'ps = {ps}')

        LC = linear_comb(ps, monomial_orbits)
        print(f'{poly_string(ps)} = {LC}')

        # Throw an error if LC != 0
        if LC != 0:
            raise RuntimeError(f'Expected linear combination to give zero polynomial, got {LC} instead')


# Compute orbits for monomials from degree 2 up to degree dmax. Initialize the
# dict with the constant "1" term. Note that Orb(1) = 3 technically but we
# are free to "normalize" any orbit by multiplying it by a scalar, so we just
# leave it as 1. The num_polys string can have any of the following values:
# "all" = computes all polynomials for every degree up to/including dmax
# "LI" = computes a_d polynomials for each degree up to/including dmax
# "nullspace" = computes a_d polynomials for all degrees < dmax, and a_d+1 polynomials of degree dmax
# Note: the polynomials are always computed in descending order of the x power.
def compute_monomial_orbits(dmax, num_polys):
    a, b = sy.sympify('a, b')
    x = [(a,b), (1-a-b,a), (b, 1-a-b)]
    # Create a defaultdict here since we want to "accumulate" into the
    # different keys below without worrying about whether an entry
    # already exists. Note that the default value for int (0) is
    # sufficient for us to use here, even though we will actually be
    # storing sympy objects in the dict.
    monomial_orbits = defaultdict(int)
    monomial_orbits['1'] = 1
    for d in range(2, dmax+1):

        # Determine the number of polynomials we are computing orbits for
        np = None
        if num_polys == "all":
            np = d+1
        elif num_polys == "LI":
            np = a_d(d)
        elif num_polys == "nullspace":
            np = a_d(d)
            if d == dmax:
                np += 1
        else:
            raise RuntimeError(f"num_polys must be one of 'all', 'LI', or 'nullspace'")

        # Debugging:
        # print(f'd = {d}, np = {np}')

        for ypower in range(np):
            for i in range(3):
                # The xpower is the complement of the ypower
                xpower = d - ypower

                # Create dict key.
                key = key_from_powers(xpower, ypower)

                # Shorthand notation for current Ro3 orbit values to evaluate
                xi = x[i][0]
                yi = x[i][1]

                # Accumulate orbit contribution
                monomial_orbits[key] += (xi**xpower) * (yi**ypower)

    # Debugging
    # print('monomial_orbits =')
    # for key, val in monomial_orbits.items():
    #     print(f'  {key} = {sy.expand(val)}')

    return monomial_orbits



"""
This code can be used to find linear combinations of monomial orbits
which are linearly dependent. For a given monomial
x^m * y^n, with m+n=d, we compute

Orb(x^m * y^n) := sum_{q=1}^3 (x_q)^m * (y_q)^n

where:

(x_q, y_q) \in {(a,b), (1-a-b,a), (b, 1-a-b)}

is the set of Ro3 points generated from the arbitrary input (a,b).
This gives a polynomial of total degree d in a and b. We can then find
a linear combination of a_d + 1 such polynomials which is linearly
dependent, and this confirms that we need only consider a_d polynomials
for each d.
"""

# Compute orbits for monomials from degree 2 up to degree dmax.
dmax = 10
monomial_orbits = compute_monomial_orbits(dmax, "all")

print('--------------------------------------------------------------------------------')
print('I. Linear combinations of Ro3-invariant polynomials that sum to zero (are not LI)')
print('--------------------------------------------------------------------------------')

# Quadratic (a_2 = 1)
# Orb(x**2) + 2*Orb(x*y) - 1 = 0
quadratic_coeffs = {'x2':1, 'x1y1':2, '1':-1}
print('> Linear combination of quadratic monomial orbits which is not LI:')
print(f'{poly_string(quadratic_coeffs)} = {linear_comb(quadratic_coeffs, monomial_orbits)}')
print('')

# Cubic (a_3 = 2)
# Orb(x**3) + Orb(x**2y) + Orb(xy**2) - Orb(x**2) = 0
cubic_coeffs = {'x3':1, 'x2y1':1, 'x1y2':1, 'x2':-1}
print('> Linear combination of cubic monomial orbits which is not LI:')
print(f'{poly_string(cubic_coeffs)} = {linear_comb(cubic_coeffs, monomial_orbits)}')
print('')

# Quartic (a_4 = 1)
print('> Linear combination of quartic monomial orbits which is not LI:')
dmax = 4
monomial_orbits_nullspace = compute_monomial_orbits(dmax, "nullspace")
compute_nullspace2(dmax, monomial_orbits_nullspace)
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
seventh_coeffs = {'x7':18, 'x6y1':63, 'x5y2':63, 'x6':-63, 'x5y1':-168, 'x4y2':-105, 'x5':63, 'x4y1':105, 'x3':-35, 'x2':21, '1':-4}
print('> Linear combination of seventh order monomial orbits which is not LI:')
print(f'{poly_string(seventh_coeffs)} = {linear_comb(seventh_coeffs, monomial_orbits)}')
print('')

# Compute nullspaces for degrees up to dmax
# print('--------------------------------------------------------------------------------')
# print(f'II. Nullspaces for degrees up to {dmax}')
# print('--------------------------------------------------------------------------------')
# for d in range(2,dmax+1):
#     compute_nullspace(d, monomial_orbits)


# Test a_d() function
# print(f'a_0 = {a_d(0)} (should be 1)')
# print(f'a_1 = {a_d(1)} (should be 0)')
# print(f'a_2 = {a_d(2)} (should be 1)')
# print(f'a_3 = {a_d(3)} (should be 2)')
# print(f'a_4 = {a_d(4)} (should be 1)')
# print(f'a_5 = {a_d(5)} (should be 2)')
# print(f'a_6 = {a_d(6)} (should be 3)')
# print(f'a_7 = {a_d(7)} (should be 2)')
# print(f'a_12 = {a_d(12)} (should be 5)')
# print(f'a_17 = {a_d(17)} (should be 6)')
# print(f'a_24 = {a_d(24)} (should be 9)')

# Test s_d() function. Apparently s_d = (d*d + 3*d + 6)/6 as well?
# for d in range(dmax+1):
#     print(f's_{d} = {s_d(d)}, s_d_simple({d}) = {s_d_simple(d)}')
