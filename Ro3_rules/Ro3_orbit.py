#!/usr/bin/env python3

import sympy as sy
from fractions import Fraction

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

# Check whether the orbit:
# {(a,b), (1-a-b,a), (b, 1-a-b)}
# gives the same quadrature contribution for some different polynomials

a, b = sy.sympify('a, b')
x = [(a,b), (1-a-b,a), (b, 1-a-b)]

# Compute orbits for monomials from degree 2 up to degree dmax. Initialize the
# dict with the constant "1" term. Note that Orb(1) = 3 technically but we
# are free to "normalize" any orbit by multiplying it by a scalar, so we just
# leave it as 1.
dmax = 5
monomial_orbits = {'1':1}
for d in range(2, dmax+1):
    for ypower in range(d+1):
        for i in range(3):
            # The xpower is the complement of the ypower
            xpower = d - ypower

            # Create dict key. We omit terms with 0 exponent, but we
            # explicitly include 1 exponents when generating the string.
            key = ''
            if xpower > 0:
                key += 'x' + str(xpower)
            if ypower > 0:
                key += 'y' + str(ypower)

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
print('.) Linear combination of quadratic monomial orbits which is not LI:')
print(f'{poly_string(quadratic_coeffs)} = {linear_comb(quadratic_coeffs, monomial_orbits)}')
print('')

# Cubic
# Orb(x**3) + Orb(x**2y) + Orb(xy**2) - Orb(x**2) = 0
cubic_coeffs = {'x3':1, 'x2y1':1, 'x1y2':1, 'x2':-1}
print('.) Linear combination of cubic monomial orbits which is not LI:')
print(f'{poly_string(cubic_coeffs)} = {linear_comb(cubic_coeffs, monomial_orbits)}')
print('')

# Quintic
# 4*Orb(x**5) + 10*Orb(x**4 y + x**3 y**2) - 5*Orb(x**4) - 10*Orb(x**2 y) + 1
quintic_coeffs = {'x5':4, 'x4y1':10, 'x3y2':10, 'x4':-5, 'x2y1':-10, '1':1}
print('.) Linear combination of quintic monomial orbits which is not LI:')
print(f'{poly_string(quintic_coeffs)} = {linear_comb(quintic_coeffs, monomial_orbits)}')
print('')
