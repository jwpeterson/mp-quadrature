import numpy
import matplotlib.pyplot as plt

def latin_hypercube_2d_uniform(n):
    # For each i, the interval (lower_limits[i], upper_limits[i]) defines a "bin"
    # The total number of two-dimensional (x,y) sample points will be n**2
    # Note: this also assumes the range of both parameters is [0,1], but that
    # is fine because we can always map them to some other range after the random
    # selection is made.
    lower_limits = numpy.arange(0,n) / float(n)
    upper_limits = numpy.arange(1,n+1) / float(n)
    # Now pick two (note the 'size' arg) random numbers from each bin.
    points = numpy.random.uniform(low=lower_limits, high=upper_limits, size=[2,n]).T
    # print('points={}'.format(points))
    # Shuffle only the second column of points so that when treated as (x,y)
    # ordered pairs we don't just get the "diagonal" set of bins. We are still
    # guaranteed to have only one entry in each y bin.
    numpy.random.shuffle(points[:, 1])
    return points

# Same thing, but takes the number of dimensions as an argument too.
def latin_hypercube_uniform(n_samples, n_dim):
    # Again, assume all ranges are [0,1]
    lower_limits = numpy.arange(0,n_samples) / float(n_samples)
    upper_limits = numpy.arange(1,n_samples+1) / float(n_samples)
    # print('lower_limits={}'.format(lower_limits))
    # print('upper_limits={}'.format(upper_limits))

    # Now pick n_dim (note the 'size' arg) random numbers from each bin.
    points = numpy.random.uniform(low=lower_limits, high=upper_limits, size=[n_dim,n_samples]).T
    # print('before shuffle points={}'.format(points))

    # Shuffle all columns after the first column.
    numpy.random.shuffle(points[:, 1:])
    # print('after shuffle points={}'.format(points))

    # Each row of "points" is a sample in n_dim-dimensional space
    return points

n = 10
# p = latin_hypercube_2d_uniform(n)
p = latin_hypercube_uniform(10, 3)
print('p={}'.format(p))

# plt.figure(figsize=[5,5])
# plt.xlim([0,1])
# plt.ylim([0,1])
# plt.scatter(p[:, 0], p[:, 1], c='r')
# plt.show()
