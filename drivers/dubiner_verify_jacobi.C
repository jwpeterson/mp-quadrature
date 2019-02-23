#include "common_definitions.h"
#include "dubiner.h"
#include "conical.h"

// This test driver verifies the orthogonality of the Dubiner
// polynomials, both symbolic and numeric versions.
int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  Dubiner dubiner;
  dubiner.compare();

  return 0;
}
