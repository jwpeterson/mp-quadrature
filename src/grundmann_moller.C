#include <cstdlib> // std::abort
#include <assert.h>
#include <sstream>
#include "grundmann_moller.h"

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

  // Used to form a string for constructing the multi-precision
  // numbers since this should work with real and rational number
  // representations.
  std::ostringstream number_stream;

  // Implement the GM formula 4.1 on page 286 of the paper
  for (unsigned i=0; i<=s; ++i)
    {
      // Get all the ordered compositions (and their permutations)
      // of |beta| = s-i into dim+1=4 parts
      this->compose_all(s-i, dim+1, permutations);

      for (unsigned p=0; p<permutations.size(); ++p)
        {
          // We use the first dim=3 entries of each permutation to
          // construct an integration point.
          for (unsigned int j=0; j<3; ++j)
            {
              // Reset the stream
              number_stream.str("");

              // Create a string with numerator/denominator
              number_stream << 2*permutations[p][j] + 1 << "/" << degree + dim - 2*i;

              // Debugging
              // std::cout << "number_stream.str()=" << number_stream.str() << std::endl;

              // Assign the multi-precision value using a string.
              // Only an mpq_class can be constructed from this type
              // of string, the mpfr_class constructor throws an
              // exception...
              mpq_class rational(number_stream.str());

              // Debugging
              // std::cout << "rational = " << rational << std::endl;

              // But it is valid to assign a rational number into a real number.
              x[offset+p](j) = rational;
            }
        }

      // This for loop needs to run for dim, degree, or dim+degree-i iterations,
      // whichever is largest.
      const unsigned weight_loop_index =
        std::max(dim, std::max(degree, degree+dim-i))+1;

      // Compute the weight for this i
      mpfr_class weight = one_pm;

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

      // This is the weight for each of the points computed previously
      for (unsigned int j=0; j<permutations.size(); ++j)
        w[offset+j] = weight;

      // Change sign for next iteration
      one_pm = -one_pm;

      // Update offset for the next set of points
      offset += permutations.size();
    }
}




void GrundmannMoller::compose_all(unsigned s,
                                  unsigned p,
                                  std::vector<std::vector<unsigned> >& result)
{
  // Clear out results remaining from previous calls
  result.clear();

  // Allocate storage for a workspace.  The workspace will periodically
  // be copied into the result container.
  std::vector<unsigned int> workspace(p);

  // The first result is always (s,0,...,0)
  workspace[0] = s;
  result.push_back(workspace);

  // the value of the first non-zero entry
  unsigned int head_value=s;

  // When head_index=-1, it refers to "off the front" of the array.  Therefore,
  // this needs to be a regular int rather than unsigned.  I initially tried to
  // do this with head_index unsigned and an else statement below, but then there
  // is the special case: (1,0,...,0) which does not work correctly.
  int head_index = -1;

  // At the end, all the entries will be in the final slot of workspace
  while (workspace.back() != s)
    {
      // If the previous head value is still larger than 1, reset the index
      // to "off the front" of the array
      if (head_value > 1)
        head_index = -1;

      // Either move the index onto the front of the array or on to
      // the next value.
      head_index++;

      // Get current value of the head entry
      head_value = workspace[head_index];

      // Put a zero into the head_index of the array.  If head_index==0,
      // this will be overwritten in the next line with head_value-1.
      workspace[head_index] = 0;

      // The initial entry gets the current head value, minus 1.
      // If head_value > 1, the next loop iteration will start back
      // at workspace[0] again.
      assert (head_value > 0);
      workspace[0] = head_value - 1;

      // Increment the head+1 value
      workspace[head_index+1] += 1;

      // Save this composition in the results
      result.push_back(workspace);
    }
}
