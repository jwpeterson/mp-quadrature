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

        poly_string += key
        count += 1
    return poly_string

# Check whether the orbit:
# {(a,b), (1-a-b,a), (b, 1-a-b)}
# gives the same quadrature contribution for some different polynomials

a, b = sy.sympify('a, b')
x = [(a,b), (1-a-b,a), (b, 1-a-b)]

# Compute monomials from degree 2 up to degree dmax. Initialize the
# dict with the constant "1" term.
dmax = 5
monomials = {'1':1}
for d in range(2, dmax+1):
    for ypower in range(d+1):
        for i in range(3):
            # The xpower is the complement of the ypower
            xpower = d - ypower

            # Create dict key. We omit terms with 0 exponent, but we
            # explicitly include 1 exponents.
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
            # monomials[key] += (xi**xpower) * (yi**ypower)

            # Doing this in 1 line doesn't work for some reason, I think
            # it is the += operator? The error is:
            # SyntaxError: can't assign to function call
            # monomials.setdefault(key, 0) += (xi**xpower) * (yi**ypower)

            # So we use two lines
            monomial = monomials.setdefault(key, 0)
            monomials[key] += (xi**xpower) * (yi**ypower)

# Debugging
print('monomials =')
for key, val in monomials.items():
    print(f'  {key} = {sy.expand(val)}')

# Cubic
# x**3 + x**2y + xy**2 - x**2 = 0
cubic_coeffs = {'x3':1, 'x2y1':1, 'x1y2':1, 'x2':-1}
print('')
print('Results:')
print('1.) Linear combination of cubic polynomials which is not LI:')
print(f'{poly_string(cubic_coeffs)} = {linear_comb(cubic_coeffs, monomials)}')
print('')

# Quintic
# 4 x**5 + 10 (x**4 y + x**3 y**2) - 5 x**4 - 10 x**2 y + 1
quintic_coeffs = {'x5':4, 'x4y1':10, 'x3y2':10, 'x4':-5, 'x2y1':-10, '1':1}
print('2.) Linear combination of quintic polynomials which is not LI:')
print(f'{poly_string(quintic_coeffs)} = {linear_comb(quintic_coeffs, monomials)}')
print('')
