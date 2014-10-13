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

// iota fills up a range with the entries (value, value+1, value+2, ...)
// and is a handy way to initialize an array.
template <typename ForwardIter, typename T>
void iota (ForwardIter first, ForwardIter last, T value)
{
  while (first != last)
    {
      *first = value++;
      ++first;
    }
}

#endif
