#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, KDTree, voronoi_plot_2d
from math import sqrt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.colorbar as cbar
import sys

# Returns True if the point (x,y) is inside the equilateral triangle,
# otherwise returns False. This is a version of the check in
# src/geom/face_tri3.C in libmesh which has been simplified to the
# case of an equilateral triangle.
def contains_point(x,y):
    # Treat points which are outside the triangle by less than tol as
    # being inside the triangle. We generally want to err on the side
    # of keeping more regions, since even if we accidentally get some
    # Voronoi regions outside the triangle, we are going to cover them
    # up using plt.fill() calls anyway.
    tol = 1.e-3
    dot00 = 1
    dot01 = 0.5
    dot02 = x
    dot11 = 1
    dot12 = 0.5*(x + sqrt(3.)*y)

    # Compute barycentric coords
    inv_denom = 1 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom

    # Check if point is in triangle
    return (u > -tol) and (v > -tol) and (u + v < 1 + tol)

# Convert a numerical value 'f' to a hex color using an 'f2rgb'
# object as returned by cm.ScalarMappable
# https://stackoverflow.com/questions/27206300/assign-color-to-value
def f2hex(f2rgb, f):
    rgb = f2rgb.to_rgba(f)[:3]
    return '#%02x%02x%02x' % tuple([255*fc for fc in rgb])

################################################################################

# Quadrature point data generated by calling
# ./drivers/tri_rule -w -i inputs/some_file.in

# d=3

# 4 QPs
# filename = 'quad_2d_p03_NI_equilateral.csv'
# filename = 'quad_2d_p03_CP_equilateral.csv'

# 6 QPs
# filename = 'quad_2d_p03_Ro3_00020_06QP_equilateral.csv'
# filename = 'quad_2d_p03_Ro3_00020_06QP_second_equilateral.csv'
# filename = 'quad_2d_p03_Ro3_00020_06QP_third_equilateral.csv'
filename = 'quad_2d_p03_Ro3_00020_06QP_PI_equilateral.csv' # solved for by hand!
# filename = 'quad_2d_p03_Ro3_00110_06QP_equilateral.csv'
# filename = 'quad_2d_p03_Ro3_00110_06QP_Case2a_equilateral.csv'
# filename = 'quad_2d_p03_Ro3_00110_06QP_Case2b_equilateral.csv' # negative weights

# 7 QPs
# filename = 'quad_2d_p03_Ro3_11010_07QP_equilateral.csv'
# filename = 'quad_2d_p03_Ro3_11010_07QP_second_equilateral.csv'
# filename = 'quad_2d_p03_Ro3_11010_07QP_third_equilateral.csv'
# filename = 'quad_2d_p03_Ro3_11100_07QP_equilateral.csv'

# d=4
# filename = 'quad_2d_p04_06QP_equilateral.csv'
# filename = 'quad_2d_p04_Ro3_0011_06QP_equilateral.csv'

# d=5. These are actually the same. In the Ro3-invariant case, the
# two general orbits converge to "median" orbits and the result matches
# the D3-invariant case.
# filename = 'quad_2d_p05_Ro3_1002_07QP_equilateral.csv'
# filename = 'quad_2d_p05_07QP_equilateral.csv'

# d=6
# D3-invariant rules.
# filename = 'quad_2d_p06_dunavant_12QP_equilateral.csv'
 #filename = 'quad_2d_p06_hompack1_12QP_equilateral.csv'
# Ro3-invariant rules
# filename = 'quad_2d_p06_Ro3_0022_12QP_equilateral.csv'

# d=7
# filename = 'quad_2d_p07_Ro3_1013_13QP_equilateral.csv'
# filename = 'quad_2d_p07_Ro3_0032_15QP_equilateral.csv'
# filename = 'quad_2d_p07_Ro3_1122_16QP_equilateral.csv'
# d=7, Gatermann's Ro3-invariant  rule
# filename = 'quad_2d_p07_gatermann_12QP_equilateral.csv'
# d=7, D3-invariant rules
# filename = 'quad_2d_p07_mine_15QP_equilateral.csv'
# filename = 'quad_2d_p07_hompack1_15QP_equilateral.csv'
# filename = 'quad_2d_p07_dunavant_13QP_equilateral.csv'
# filename = 'quad_2d_p07_zhang_underdetermined_15QP_equilateral.csv'

# d=8
# filename = 'quad_2d_p08_Ro3_1014_16QP_equilateral.csv'
# filename = 'quad_2d_p08_Ro3_0114_18QP_equilateral.csv'
# filename = 'quad_2d_p08_Ro3_1123_19QP_equilateral.csv'
# d=8, D3-invariant rule
# filename = 'quad_2d_p08_16QP_equilateral.csv'

# d=9
# filename = 'quad_2d_p09_Ro3_1006_19QP_equilateral.csv' # same as D3-invariant rule.
# filename = 'quad_2d_p09_Ro3_0106_21QP_equilateral.csv'
# filename = 'quad_2d_p09_Ro3_1115_22QP_equilateral.csv'
# filename = 'quad_2d_p09_Ro3_0134_24QP_equilateral.csv'
# d=9, D3-invariant rule
# filename = 'quad_2d_p09_19QP_equilateral.csv'

# d=10
# filename = 'quad_2d_p10_Ro3_0026_24QP_equilateral.csv'
# filename = 'quad_2d_p10_Ro3_0026_24QP_second_equilateral.csv'
# filename = 'quad_2d_p10_Ro3_0026_24QP_third_equilateral.csv'
# filename = 'quad_2d_p10_Ro3_0026_24QP_fourth_equilateral.csv'
# filename = 'quad_2d_p10_Ro3_0026_24QP_fifth_equilateral.csv'
# filename = 'quad_2d_p10_Ro3_0026_24QP_sixth_equilateral.csv'

# filename = 'quad_2d_p10_Ro3_1116_25QP_equilateral.csv'
# filename = 'quad_2d_p10_Ro3_1035_25QP_equilateral.csv'
# filename = 'quad_2d_p10_Ro3_1035_25QP_second_equilateral.csv'
# filename = 'quad_2d_p10_Ro3_0135_27QP_equilateral.csv'

# d=10, D3-invariant rules (found before or previously known)
# filename = 'quad_2d_p10_hompack1_25QP_equilateral.csv'
# filename = 'quad_2d_p10_hompack2_25QP_equilateral.csv'
# filename = 'quad_2d_p10_dunavant_25QP_equilateral.csv'
# filename = 'quad_2d_p10_zhang_underdetermined_25QP_equilateral.csv'

# d=11, Ro3-invariant rules
# filename = 'quad_2d_p11_Ro3_0018_27QP_equilateral.csv'
# filename = 'quad_2d_p11_Ro3_0018_27QP_second_equilateral.csv'
# filename = 'quad_2d_p11_Ro3_1027_28QP_equilateral.csv'
# filename = 'quad_2d_p11_Ro3_0127_30QP_equilateral.csv'
# filename = 'quad_2d_p11_Ro3_0046_30QP_equilateral.csv'
# filename = 'quad_2d_p11_Ro3_1136_31QP_equilateral.csv'
# filename = 'quad_2d_p11_Ro3_1136_31QP_second_equilateral.csv'

# d=11, D3-invariant rules
# filename = 'quad_2d_p11_hompack1_30QP_equilateral.csv'
# filename = 'quad_2d_p11_hompack2_30QP_equilateral.csv'
# filename = 'quad_2d_p11_mine_30QP_equilateral.csv'
# filename = 'quad_2d_p11_zhang_underdetermined_28QP_equilateral.csv'

# d=12, D3-invariant rules
# filename = 'quad_2d_p12_dunavant_33QP_equilateral.csv'
# filename = 'quad_2d_p12_mine_33QP_equilateral.csv'

# d=13, D3-invariant rules
# filename = 'quad_2d_p13_zhang_37QP_equilateral.csv'

# d=14, Ro3-invariant rules
# filename = 'quad_2d_p14_Ro3_1-0-3-11_43QP_equilateral.csv'
# filename = 'quad_2d_p14_Ro3_0-1-3-11_45QP_equilateral.csv'

# d=14, D3-invariant rule. This rule appears to have 6 points on the
# boundary, but they are just very close and not actually on it.
# filename = 'quad_2d_p14_dunavant_42QP_equilateral.csv'

# d=15, D3-invariant rules. The Wandzura rule is rank-deficient.
# filename = 'quad_2d_p15_mine_49QP_equilateral.csv'
# filename = 'quad_2d_p15_mine_51QP_equilateral.csv'
# filename = 'quad_2d_p15_mine_52QP_equilateral.csv'
# filename = 'quad_2d_p15_zhang_witherden_49QP_equilateral.csv'
# filename = 'quad_2d_p15_wandzura_54QP_equilateral.csv'

# d=16, D3-invariant rules.
# filename = 'quad_2d_p16_mine_55QP_equilateral.csv'
# filename = 'quad_2d_p16_zhang_underdetermined_55QP_equilateral.csv'

# d=17, D3-invariant rules.
# filename = 'quad_2d_p17_mine_63QP_equilateral.csv'
# filename = 'quad_2d_p17_mine_63QP_second_equilateral.csv'
# filename = 'quad_2d_p17_dunavant_underdetermined_61QP_equilateral.csv'
# filename = 'quad_2d_p17_zhang_xiao_underdetermined_60QP_equilateral.csv'

################################################################################

# Read CSV data into an array
data = np.genfromtxt(filename, delimiter=',')

# Extract first two columns of data (assuming there is a third column of weights
points = data[:, [0, 1]]
wts = data[:, 2]

# print(points)
# print('wts={}'.format(wts))

# Create a colormap. This is mostly a matter of preference?
cmap = cm.get_cmap('RdYlGn')
# cmap = cm.get_cmap(plt.cm.jet)

# Perceptually uniform sequential colormaps
# cmap = cm.get_cmap(plt.cm.viridis) # Can't see black QPs on purple background
# cmap = cm.get_cmap(plt.cm.plasma)
# cmap = cm.get_cmap(plt.cm.inferno)
# cmap = cm.get_cmap(plt.cm.magma)

# Create a numerical value -> color mapping based on the quadrature weights.
# We use f2rgb when calling 'f2hex', see below.
norm = colors.Normalize(vmin=np.min(wts), vmax=np.max(wts))
f2rgb = cm.ScalarMappable(norm=norm, cmap=cmap)

# To create a "bounded" Voronoi region, we first reflect the original
# points across each edge of the equilateral triangle, then we take
# all those points and reflect them across the edges of the "macro"
# element. There may be better ways of doing it, this is just and idea
# I came across at https://stackoverflow.com/questions/28665491

# Slope "m" and y-intercept "b" of the edges of the equilateral triangle
# in the formula y=mx + b
s3 = sqrt(3.)

# First set of reflections across equilatral triangle edges.
mb1 = [[0, 0], [-s3, s3], [s3, 0]]

# Second set of reflections across "macro element" edges.
mb2 = [[0, s3/2], [s3, -s3], [-s3, 0]]

# First we reflect points across each edge of the original triangle.
# Then we take that *entire set of points and reflect it across each
# edge of the macro triangle.
combined_mb = [mb1, mb2]

# Save a copy of the original set of points before we start adding a bunch
# of reflections.
orig_points = points

for refset in combined_mb:
    # Keep track of transformed points. Numpy arrays are not really
    # designed to be appended to, so we keep track of a list of
    # transformed points to be added and then concatenate them all at the
    # end.
    transformed_points = []
    for mbpair in refset:
        m = mbpair[0]
        b = mbpair[1]
        for i in xrange(len(points)):
            # (x,y) position of original point
            x = points[i, 0]
            y = points[i, 1]

            # Transformed positions
            u = ((1 - m**2)*x + 2*m*y - 2*m*b) / (m**2 + 1)
            v = ((m**2 - 1)*y + 2*m*x + 2*b) / (m**2 + 1)

            # Don't include points which are right on the edge when
            # transforming points.
            if sqrt((x-u)**2 + (y-v)**2) >= 1.e-6:
                transformed_points.append([u, v])

    # Concatenate the lists of original and transformed points.
    points = np.concatenate((points, transformed_points), axis=0)
    # print('shape(points)={}'.format(np.shape(points)))

# Generate the Voroni diagram
vor = Voronoi(points)

# Debugging: no filtering, just plot Voronoi diagram as-is and exit.
# fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange',
#                       line_width=2, line_alpha=0.6, point_size=2)
# plt.savefig('voronoi.pdf')
# sys.exit()

# Filter regions. Here we drop any regions which have a point at
# infinity or those which have all points outside the equilateral
# triangle. Again, this filtering isn't perfect, but it's OK because
# we are going to cover up the regions which "stick out" later, and it
# does make generating the plots a bit faster when you don't have so
# many regions to plot.
regions = []
for region in vor.regions:
    has_infinite = False
    all_points_outside = True

    for index in region:
        # Discard any region with point(s) at infinity
        if index == -1:
            has_infinite = True
            break
        else:
            # Check whether at least 1 point is inside the equalateral triangle.
            x = vor.vertices[index, 0]
            y = vor.vertices[index, 1]
            if contains_point(x,y):
                all_points_outside = False
                break

            # Debugging: keep every region
            # all_points_outside = False

    if region != [] and not has_infinite and not all_points_outside:
        regions.append(region)

# Overwrite the original set of regions. Note that this doesn't perfectly get
# rid of *all* the regions we want to get rid of. But that's OK because we are
# just going to cover up the remaining ones with filled in polygons in the end
# anyway...
vor.regions = regions

# The top of the triangle is at y = sqrt(3)/2 ~ 0.8660254037844386.
# If we want the axes to be square, we need to use the same limits for
# both x and y. I found that does not really look as good, however.
border=0.05
plt.xlim(0-border, 1+border)
plt.ylim(0-border, s3/2+border)

# Plot edges which have at least one vertex in/on the equilateral triangle.
# Note: by making these lines white we are kind of cheating, since it allows
# us to hide the ones that are outside the region, but this is easier than
# trying to find where the Voronoi edges intersect the triangle and somehow
# cutting them...
for simplex in vor.ridge_vertices:
    simplex = np.asarray(simplex)

    # Skip any edges with vertices at infinity (-1 indices). We have
    # not filtered these out yet, even though we did filter out
    # regions with points at infinity.
    if np.any(simplex < 0):
        continue

    # Extract (x,y) positions to pass into contains_point()
    x0 = vor.vertices[simplex, 0][0]
    x1 = vor.vertices[simplex, 0][1]
    y0 = vor.vertices[simplex, 1][0]
    y1 = vor.vertices[simplex, 1][1]
    if (contains_point(x0, y0) or contains_point(x1, y1)):
        plt.plot([x0, x1], [y0, y1], '-', color='white', linewidth=1)

# Build a KDTree from the original points. We will use this
# later to figure out quickly which region maps to which QP.
tree = KDTree(orig_points)

# Here's code to fill in the non-infinite regions! If you don't specify a
# color to the 'fill' command, it chooses a random sequence of colors.
# The regions are not numbered in the same order as the QPs, so we use
# the KDTree to find the closest QP to the centroid of each region, and
# assign that region a corresponding color.
for region in vor.regions:
    # We should have got rid of all the infinite regions by now.
    if -1 in region:
        raise RuntimeError('Infinite region.')

    polygon = []
    centroid = np.zeros((1, 2))
    for i in region:
        polygon.append(vor.vertices[i])
        # print('type of vor.vertices[i] is {}'.format(type(vor.vertices[[i]])))
        # print('vor.vertices[i]={}'.format(vor.vertices[i]))
        centroid = np.add(centroid, vor.vertices[i])
        # centroid[0] = centroid[0] + vor.vertices[i][0]
        # centroid[1] = centroid[1] + vor.vertices[i][1]

    # Scale by number of vertices in the region
    centroid = np.divide(centroid, float(len(region)))
    # print('centroid = {}'.format(centroid))

    # Look for the closest corresponding point in the KDTree
    result = tree.query(centroid)
    dist = np.asscalar(result[0])
    index = np.asscalar(result[1])
    # print('kdtree: dist={}, index={}, pts[index]=({},{}), wts[index]={}'.
    #       format(dist,
    #              index,
    #              orig_points[index][0],
    #              orig_points[index][1],
    #              wts[index]))

    # print(polygon)
    plt.fill(*zip(*polygon), color=f2hex(f2rgb, wts[index]))

# Fill the outside of the triangle with white regions to cover
# up the Voronoi elements that extend outside the domain.
# I don't really understand exactly how this works, but if
# you have a list of 2 element numpy arrays representing the
# vertices of the polygon, and then you call zip() on that, this
# is the input expected by plt.fill()

# Fill on left side of the triangle.
polygon=[np.array([-1,-s3]), np.array([1,s3]), np.array([-1,s3])]
plt.fill(*zip(*polygon), color='white')

# Fill below the triangle
polygon=[np.array([-1,0]), np.array([2,0]), np.array([2,-1]), np.array([-1,-1])]
plt.fill(*zip(*polygon), color='white')

# Fill to the right of the triangle
polygon=[np.array([2,-s3]), np.array([2,s3]), np.array([0,s3])]
plt.fill(*zip(*polygon), color='white')

# Plot lines for the edges of the equilateral triangle
ax = plt.gca()
ax.plot([0, 1, 0.5, 0], [0, 0, sqrt(3.)/2., 0], color='k', linestyle='-', linewidth=2)

# Plot the original points.
plt.plot(orig_points[:, 0], orig_points[:, 1], 'ko')

# Add a colorbar to the plot. I have no idea how this syntax works, I
# found it at https://stackoverflow.com/questions/17660071
cax, _ = cbar.make_axes(ax, shrink=0.6)
cb2 = cbar.ColorbarBase(cax, cmap=cmap, norm=norm)

# Save a PDF. Generate an output filename from the input filename.
output_filename = filename.rsplit( ".", 1)[0] + '.pdf'
plt.savefig(output_filename)
