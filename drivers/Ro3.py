import numpy as np
from scipy.optimize import minimize
import scipy.special

# Integral of x^a y^b over the reference triangle
def exact_tri(a, b):
    return 1. / (b+1) * scipy.special.beta(a+1, b+2)

# Polynomial x**a y**b powers to be tested
polys = [(0,0),
         (2,0),                  # 2nd
         (3,0),  (2,1),          # 3rd
         (4,0),                  # 4th
         (5,0),  (4,1),          # 5th
         (6,0),  (5,1),  (4,2),  # 6th
         (7,0),  (6,1)]          # 7th

# This function computes the residual, then "f" is just 0.5 * dot(r,r)
def res(x):
    # Residual vector
    r = np.zeros(len(x))

    # FIXME: Currently assumes no centroid point.
    for i in xrange(0, len(polys)):
        xpower = polys[i][0]
        ypower = polys[i][1]

        for q in xrange(0, len(x), 3):
            wq = x[q]
            xq = x[q+1]
            yq = x[q+2]
            zq = 1 - xq - yq
            spatial = np.power(xq, xpower) * np.power(yq, ypower) + \
                      np.power(zq, xpower) * np.power(xq, ypower) + \
                      np.power(yq, xpower) * np.power(zq, ypower)
            # Accumulate ith residual contribution
            r[i] += wq * spatial

        # Subtract exact integral value
        r[i] -= exact_tri(xpower, ypower)
    return r

def myfunc(x):
    r = res(x)
    return 0.5 * np.dot(r,r)

# Test exact_tri()
# for p in polys:
#     print('exact_tri({},{})={}'.format(p[0], p[1], exact_tri(p[0],p[1])))

# Bounds for Ro3 weight parameters are (0, 1./6) since they come in sets of 3.
wmax = 1./6

# Vector of initial guesses. rand() returns values in (0,1) and then
# we scale the weights so they are feasible.
x0 = np.random.rand(1,12)
for q in xrange(0, len(x0)):
    x0[q] *= wmax

# print('initial guess={}'.format(x0))

# The 'ftol' parameter of minimize is related to the 'factr' parameter
# of the L-BFGS-B algorithm according to the formula:
# ftol = factr * numpy.finfo(float).eps
result = minimize(myfunc, x0, method='L-BFGS-B', \
               bounds=[(0,wmax), (0,1), (0,1), \
                       (0,wmax), (0,1), (0,1), \
                       (0,wmax), (0,1), (0,1), \
                       (0,wmax), (0,1), (0,1)],
               options={'disp': True,
                        'ftol' : 1.e-30,
                        'gtol' : 1.e-9,
                        'maxls' : 30,
                        'maxcor' : 20})

# The solution is x_i = 1 for all i.
print('Solution = {}'.format(result.x))

# Check sum of weights
sum_of_weights = 0.
for q in xrange(0, len(result.x), 3):
    sum_of_weights += result.x[q]
print('sum_of_weights={}'.format(sum_of_weights))

# Check that for each (x,y) pair, 1-x-y >= 0. Nothing a priori
# enforces this condition, so it might not be satisfied.
for q in xrange(0, len(result.x), 3):
    xq = result.x[q+1]
    yq = result.x[q+2]
    zq = 1-xq-yq
    if (zq < 0.):
        print('Points found outside element: q={}, zq=1-xq-yq={}'.format(q,1-xq-yq))

# Compute final residual
r = res(result.x)
print('Final residual norm={}'.format(np.linalg.norm(r)))
