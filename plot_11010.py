import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 12
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True

"""
For configuration (1,1,0,1,0), plots the curve:

(wc(alpha), wv(alpha))

for 0 < alpha < 1/2, where wc is the centroid weight and wv is the
vertex weight.
"""
alpha1 = np.linspace(.01, (9-np.sqrt(21))/30)
alpha2 = np.linspace((9-np.sqrt(21))/30, 0.2)
alpha3 = np.linspace(0.2, (9+np.sqrt(21))/30)
alpha4 = np.linspace((9+np.sqrt(21))/30, 0.5)

alpha = np.concatenate((alpha1,alpha2,alpha3,alpha4))
wm = 1. / 120 / alpha / (3*alpha - 1)**2
wv = (5. * alpha - 1) / (120 * alpha)
wc = 1./2 - 3*wm - 3*wv

# print('alpha={}'.format(alpha))
# print('wm={}'.format(wm))
# print('wv={}'.format(wv))
# print('wc={}'.format(wc))

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(wc, wv, color='black')

# Add grid lines. Do this first so they are "under" the plot. I did not
# find a way to do this with ax1.grid so I used manual lines.
# ax1.grid(b=True, which='major', axis='both', color='lightgray', linestyle='--', linewidth=1)
ax1.plot([0,0],[-1,1], color='lightgray', linestyle='--', linewidth=1)
ax1.plot([-1,1],[0,0], color='lightgray', linestyle='--', linewidth=1)

# Fill in upper right part of the graph, these are PB rules
plt.fill([0, 0.5, 0.5, 0], [0, 0, 0.05, 0.05], color='lightgray')

# Mark "special" points, either endpoints or where one of the weights goes to 0.
ax1.plot([0.0, -0.28125, 0.0, 9./40],
         [-0.014927398729, 0.0, 0.0232607320623, 1./40], color='black', linestyle='', marker='o')

# Put labels on the plot
ax1.text(0.0, -0.014927398729, r'$\alpha = \frac{9-\sqrt{21}}{30}$')
ax1.text(-0.28125, 0.0+.001, r'$\alpha = \frac{1}{5}$')
ax1.text(0.0-.06, 0.0232607320623+.0025, r'$\alpha = \frac{9+\sqrt{21}}{30}$')
ax1.text(9./40, 1./40, r'$\alpha = \frac{1}{2}$')
ax1.annotate(r"$\alpha\rightarrow 0$", xy=(0.22, -0.05), xytext=(0.14, -0.035),
             arrowprops=dict(arrowstyle="->"))
ax1.annotate(r"$\alpha\rightarrow \frac{1}{3}$", xy=(-0.5, 0.012), xytext=(-0.4, 0.01),
             arrowprops=dict(arrowstyle="->"))
# ax1.text(0.3, 0.01, 'PB')
ax1.text(0.1, 0.02, 'PB')
ax1.text(-0.28125-.02, 0.0-.005, 'NI')

# Set limits, labels, etc. and make plot
ax1.set_xlim([-.5, .5])
ax1.set_ylim([-0.05, 0.05])
ax1.set_xlabel(r'$w_c$')
ax1.set_ylabel(r'$w_v$')
plt.savefig('plot_11010.pdf', format='pdf')
