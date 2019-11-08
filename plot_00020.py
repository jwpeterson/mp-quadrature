import numpy as np
import matplotlib.pyplot as plt
import sys

from matplotlib import rcParams
rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 12
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True

"""
For configuration (0,0,0,2,0), plots the curve:

(w1(alpha), w2(alpha))

for r2-eps < alpha < r2+eps, where r2 = (9-sqrt(21))/30 and eps ~ .022
"""

# It's convenient to have Python functions for these polynomials.
def f(x):
    return 6*x**2 - 4*x + 0.5

def g(x):
    return x * (15*x**2 - 9*x + 1)

r1 = (9 + np.sqrt(21))/30
r2 = (9 - np.sqrt(21))/30
alphas = np.linspace(1.e-3, 0.17048)
print('alphas={}'.format(alphas))
w1 = np.zeros(len(alphas))
w2 = np.zeros(len(alphas))

for i in xrange(len(alphas)):
    alpha = alphas[i]
    # Compute sigma based on alpha
    sigma = f(alpha) / g(alpha)
    # Compute roots of "characteristic" eqn.
    roots = np.roots([15*sigma, -(9*sigma + 6), (sigma + 4), -0.5])
    print('roots={}'.format(roots))
    # Find root which is closest to r1
    x1_new = min(roots, key=lambda x : abs(x - r1))

    # Exit if we failed to find a real root
    if not np.isreal(x1_new):
        print('Imaginary root found for alpha={}, try smaller neighborhood around r2.'.format(alpha))
        sys.exit(1)

    # Now solve for the weights w1, w2. Result should be the same
    # regardless of whether we use the ratio of f's or g's to compute.
    w1_over_w2 = -f(alpha) / f(x1_new)
    # w1_over_w2 = -g(alpha) / g(x1_new)
    w2_over_w1 = 1. / w1_over_w2

    # It seems that w1/w2 = constant, the plot in (w1, w2) space is a straight line.
    # print('w1_over_w2={}'.format(w1_over_w2))
    # print('w2_over_w1={}'.format(w2_over_w1))
    w1[i] = 1. / 6 / (1 + w2_over_w1)
    w2[i] = 1. / 6 / (1 + w1_over_w2)

    # Print current result
    # print('---')
    # print('w1={}'.format(w1[i]))
    # print('x1={}'.format(x1_new))
    # print('w2={}'.format(w2[i]))
    # print('x2={}'.format(alpha))

# Print results
# print('w1={}, w2={}'.format(w1, w2))

# Plot results
fig = plt.figure()
ax1 = fig.add_subplot(111)

# Plot y=x line of symmetry
ax1.plot([0,.2],[0,.2], color='lightgray', linestyle='--', linewidth=1)

# Plot line x=0
ax1.plot([0,0],[0,.2], color='lightgray', linestyle='--', linewidth=1)
# Plot line y=0
ax1.plot([0,.2],[0,0], color='lightgray', linestyle='--', linewidth=1)

# Plot main solutions (and symmetric solutions
ax1.plot(w1, w2, color='black', marker=None)
# ax1.plot(w2, w1, color='salmon', marker=None)

# Plot other solutions.

# The PI solution: should be on the curve since we are expanding about
# the PI solution.
ax1.plot([1./12 - np.sqrt(21)/168],
         [1./12 + np.sqrt(21)/168], color='red', linestyle='', marker='o')

# The PB/NB solution pair
ax1.plot([(39 - np.sqrt(21))/240],
         [(1 + np.sqrt(21))/240], color='blue', linestyle='', marker='o')
ax1.plot([(1 - np.sqrt(21))/240],
         [(39 + np.sqrt(21))/240], color='green', linestyle='', marker='o')

# The PB solution without partner
ax1.plot([1./60],
         [3./20], color='khaki', linestyle='', marker='o')

#ax1.set_xlim([0, 0.1])
#ax1.set_ylim([0, 0.1])
ax1.set_xlabel(r'$w_1$')
ax1.set_ylabel(r'$w_2$')
ax1.axis('equal')
plt.savefig('plot_00020.pdf', format='pdf')
