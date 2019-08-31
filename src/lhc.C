// Local includes
#include "lhc.h"

// C++ includes
#include <algorithm> // std::random_shuffle
#include <stdlib.h> // random

std::vector<std::vector<double>>
lhc(unsigned int n_dim, unsigned int n_tests)
{
  std::vector<std::vector<double>> random_numbers(n_dim);
  double bin_size = 1. / double(n_tests);

  for (unsigned int i=0; i<random_numbers.size(); ++i)
    {
      random_numbers[i].resize(n_tests);

      // Get random number in the current bin.
      for (unsigned int j=0; j<random_numbers[i].size(); ++j)
        {
          // Lower bound of current bin
          double lower = j * bin_size;

          random_numbers[i][j] = bin_size * double(random()) / RAND_MAX + lower;
        }
    }

  // random_shuffle all rows after the first one. This way we don't have
  // a sample from the *same* bin for each random variable.
  for (unsigned int i=0; i<random_numbers.size(); ++i)
    std::random_shuffle(random_numbers[i].begin(), random_numbers[i].end());

  return random_numbers;
}
