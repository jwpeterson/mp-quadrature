#ifndef __common_definitions__
#define __common_definitions__

#include <string>

// For mpfr_class.  It may be possible to forward declare this, but
// I'm not sure how.
#include "gmpfrxx.h"

typedef double Real;

// Returns a string which contains the 32 decimal digits of x in the form:
// -9.4489927222288222340758013830322e-01L
std::string fix_string(const mpfr_class& x);

// Helper function to raise an arbitrary precision parameter to an integer power
mpfr_class pow(const mpfr_class& x, unsigned power);

#endif
