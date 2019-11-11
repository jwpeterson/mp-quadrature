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
# The curve goes from:
# (x1,x2) = ((9 + np.sqrt(21))/30, 0) to
# (x1,x2) = (1/2, 1/6)
# But w1/w2 -> \infty at alpha=1/6
alphas1 = np.linspace(1.e-6, 0.15)
alphas2 = np.linspace(0.15, 1./6 - 1.e-6)
alphas = np.concatenate((alphas1,alphas2))
# print('alphas={}'.format(alphas))
w1 = np.zeros(len(alphas))
w2 = np.zeros(len(alphas))
x1 = np.zeros(len(alphas))

end = len(alphas)
for i in xrange(len(alphas)):
    alpha = alphas[i]
    # Compute sigma based on alpha
    sigma = f(alpha) / g(alpha)
    # Compute roots of "characteristic" eqn.
    roots = np.roots([15*sigma, -(9*sigma + 6), (sigma + 4), -0.5])
    # print('roots={}'.format(roots))
    # Find root which is closest to r1
    candidate_x1 = min(roots, key=lambda x : abs(x - r1))

    # Don't accept candidate result if root is imaginary or outside the reference element.
    if (not np.isreal(candidate_x1)) or (candidate_x1 > 0.5):
        print('Invalid root = {} found for i={}, alpha={}, try smaller neighborhood around r2.'.format(candidate_x1, i, alpha))
        # Record the last valid index and break out of the loop
        end = i
        break

    # If we made it here, we accept the candidate solution
    x1[i] = candidate_x1

    # Now solve for the weights w1, w2. Result should be the same
    # regardless of whether we use the ratio of f's or g's to compute.
    w1_over_w2 = -f(alpha) / f(x1[i])
    # w1_over_w2 = -g(alpha) / g(x1[i])
    w2_over_w1 = 1. / w1_over_w2

    # It seems that w1/w2 = constant, the plot in (w1, w2) space is a straight line.
    # print('w1_over_w2={}'.format(w1_over_w2))
    # print('w2_over_w1={}'.format(w2_over_w1))
    w1[i] = 1. / 6 / (1 + w2_over_w1)
    w2[i] = 1. / 6 / (1 + w1_over_w2)

    # Print current result
    # print('---')
    # print('w1={}'.format(w1[i]))
    # print('x1={}'.format(x1[i]))
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

# Plot (x1(alpha), alpha)
fig = plt.figure()
ax1 = fig.add_subplot(111)
# Plot line y=1/6
ax1.plot([0.44,0.55],[1./6,1./6], color='lightgray', linestyle='--', linewidth=1)
# Plot line y=0
ax1.plot([0.44,0.55],[0.,0.], color='lightgray', linestyle='--', linewidth=1)
# Plot line x=0.5
# ax1.plot([0.5,0.5],[0,1./6], color='lightgray', linestyle='--', linewidth=1)
ax1.plot(x1[0:end], alphas[0:end], color='black', marker=None)
ax1.plot([r1], [0], color='black', linestyle='', marker='o')
ax1.plot([0.5], [1./6], color='black', linestyle='', marker='o')
ax1.text(r1+.001, 0.+.002, r'$\alpha = 0$')
ax1.text(0.5-.003, 1./6-.015, r'$\alpha = \frac{1}{6}$')
ax1.set_xlabel(r'$x_1$')
ax1.set_ylabel(r'$x_2$')
ax1.set_xlim([0.443, 0.505])
# ax1.axis('equal')
plt.savefig('plot_00020_x1.pdf', format='pdf')