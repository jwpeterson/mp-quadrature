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



mpq_class exact_tet(unsigned a, unsigned b, unsigned c)
{
  mpq_class analytical = 1;

  // Sort the a, b, c values
  unsigned sorted_powers[3] = {a, b, c};
  std::sort(sorted_powers, sorted_powers+3);

  // Cancel the largest power with the denominator, fill in the
  // entries for the remaining numerator terms and the denominator.
  std::vector<unsigned>
    numerator_1(sorted_powers[0] > 1 ? sorted_powers[0]-1 : 0),
    numerator_2(sorted_powers[1] > 1 ? sorted_powers[1]-1 : 0),
    denominator(3 + sorted_powers[0] + sorted_powers[1]);

  // Fill up the vectors with sequences starting at the right values.
  iota(numerator_1.begin(), numerator_1.end(), 2);
  iota(numerator_2.begin(), numerator_2.end(), 2);
  iota(denominator.begin(), denominator.end(), sorted_powers[2]+1);

  // The denominator is guaranteed to have the most terms...
  for (unsigned i=0; i<denominator.size(); ++i)
    {
      if (i < numerator_1.size())
        analytical *= numerator_1[i];

      if (i < numerator_2.size())
        analytical *= numerator_2[i];

      analytical /= denominator[i];
    }

  return analytical;
}



mpq_class exact_pyr(unsigned p, unsigned q)
{
  // 0 if p or q is odd
  if (p%2 != 0 || q%2 != 0)
    return 0.;

  return mpq_class(4) / mpq_class(p+q+3) / mpq_class(q+1) / mpq_class(p+1);
}
