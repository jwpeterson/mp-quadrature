#include <cstdlib> // std::abort
#include <assert.h>
#include <sstream>
#include "grundmann_moller.h"
#include "compose_all.h"

GrundmannMoller::GrundmannMoller()
{
}



void GrundmannMoller::rule(unsigned s)
{
  // Currently we only support 3D rules:
  unsigned dim = 3;

  // s is the so-called "index" of the rule, it can be useful to use in place of d
  unsigned degree = 2*s + 1;

  // The number of points in the rule
  unsigned n_pts = (s+4)*(s+3)*(s+2)*(s+1) / 24;

  // Allocate space for points and weights
  x.resize(n_pts);
  w.resize(n_pts);

  // (-1)^i -> This one flips sign at each iteration of the i-loop below.
  int one_pm=1;

  // Where we store all the integer point compositions/permutations
  std::vector<std::vector<unsigned> > permutations;

  // Index into the vector where we should start adding the next round of points/weights
  std::size_t offset=0;

  // Implement the GM formula 4.1 on page 286 of the paper
  for (unsigned i=0; i<=s; ++i)
    {
      // Get all the ordered compositions (and their permutations)
      // of |beta| = s-i into dim+1=4 parts
      compose_all(s-i, dim+1, permutations);

      for (unsigned p=0; p<permutations.size(); ++p)
        {
          // We use the first dim=3 entries of each permutation to
          // construct an integration point.
          for (unsigned int j=0; j<3; ++j)
            {
              x[offset+p](j) = 2*permutations[p][j] + 1;
              x[offset+p](j) /= degree + dim - 2*i;
            }
        }

      // This for loop needs to run for dim, degree, or dim+degree-i iterations,
      // whichever is largest.
      const unsigned weight_loop_index =
        std::max(dim, std::max(degree, degree+dim-i))+1;

      // Compute the weight for this i
      mpq_class weight = one_pm;

      for (unsigned int j=1; j<weight_loop_index; ++j)
        {
          if (j <= degree) // Accumulate (d+n-2i)^d term
            weight *= (degree + dim - 2*i);

          if (j <= 2*s) // Accumulate 2^{-2s}
            weight /= 2;

          if (j <= i) // Accumulate (i!)^{-1}
            weight /= j;

          if (j <= degree+dim-i) // Accumulate ( (d+n-i)! )^{-1}
            weight /= j;
        }

      // Debugging
      // std::cout << "weight=" << weight << std::endl;

      // This is the weight for each of the points computed previously
      for (unsigned int j=0; j<permutations.size(); ++j)
        w[offset+j] = weight;

      // Change sign for next iteration
      one_pm = -one_pm;

      // Update offset for the next set of points
      offset += permutations.size();
    }
}




