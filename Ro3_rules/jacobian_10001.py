#!/usr/bin/env python
import sympy
from sympy.abc import a, w, x, y, z
from fractions import Fraction
from sympy.matrices import *
import sys

ninth = Fraction(1,9)
one27 = Fraction(1,27)
J = Matrix(4, 4,
           [1, 3, 0, 0,
            ninth, x**2 + y**2 + z**2, 2 * w * (x - z), 2 * w * (y - z),
            one27, x**3 + y**3 + z**3, w*(3*x**2 - 3*z**2), w*(3*y**2 - 3*z**2),
            one27, x**2 * y + z**2 * x + y**2 * z, w*(2*x*y + z**2 - 2*x*z - y**2), w*(2*y*z + x**2 - 2*x*z - y**2)])
print('Jacobian = {}'.format(J))

# Compute determinant, do basic simplification.
det = J.det()
print('---')
print('det(J) = {}'.format(det))

det_xy = sympy.simplify(det.subs([(z,1-x-y)]))
print('---')
print('With z=1-x-y: {}'.format(det))

det_xy = sympy.simplify(det_xy / (2 * w**2))
print('---')
print('To simplify the problem, divide by 2w^2 (we note that w=0 makes Jacobian singular):\n{}'.format(det_xy))

# Amazingly, this works really well and makes the result much simpler to understand. Is this even true?
# Result: (3*x**2 + 3*x*y - 3*x + 3*y**2 - 3*y + 1)**3/9
# Note: according to Wofram-alpha, the value in parentheses has only one real root at (x=1/3, y=1/3). It
# is an elliptic paraboloid with values >= 0 for all input values (x,y).
# Also when y=x, the value in parens is (1-3*x)**2 = (3x-1)**2, which is similar to/consistent
# with what we obtained for the reduced Jacobian's determinant, (-4/9)*wg*(3*x-1)^4.
det_xy = sympy.factor(det_xy)
print('---')
print('sympy.factor:\n{}'.format(det_xy))

# Let's do another test to make sure this determinant is correct? If
# check==0 then I think we can be relatively confident that the factor
# approach did not lead us astray.
check = sympy.simplify(det.subs([(z,1-x-y)]) - 2 * w**2 * (3*x**2 + 3*x*y - 3*x + 3*y**2 - 3*y + 1)**3/9)
print('---')
print('check:\n{}'.format(check))
