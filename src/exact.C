#include "common_definitions.h"
#include "exact.h"

mpq_class exact_tri(unsigned a, unsigned b)
{
  mpq_class analytical = 1;

  unsigned
    larger_power = std::max(a, b),
    smaller_power = std::min(a, b);

  // Cancel the larger of the two numerator terms with the
  // denominator, and fill in the remaining entries.
  std::vector<unsigned>
    numerator(smaller_power > 1 ? smaller_power-1 : 0),
    denominator(2+smaller_power);

  // Fill up the vectors with sequences starting at the right values.
  iota(numerator.begin(), numerator.end(), 2);
  iota(denominator.begin(), denominator.end(), larger_power+1);

  // The denominator is guaranteed to have more terms...
  for (unsigned i=0; i<denominator.size(); ++i)
    {
      if (i < numerator.size())
        analytical *= numerator[i];
      analytical /= denominator[i];
    }

  return analytical;
}
