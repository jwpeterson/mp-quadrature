#include <sstream>

#include "common_definitions.h"

std::string fix_string(const mpfr_class& x)
{
  // Go to some extra trouble to print e-01L instead of e-1
  std::ostringstream number_stream;
  number_stream.precision(32);
  number_stream.setf(std::ios_base::scientific);
  number_stream << x;
  std::string number = number_stream.str();

  // We almost always want to append an L (long double literal) but not always
  bool append_L = true;

  // Find the exponent 'e'.
  size_t e_pos = number.find('e');

  // If 'e' is not found...
  if (e_pos == std::string::npos)
    {
      if (number == "0")
        {
          // if number == "0", append a decimal point
          number.append(".");
          append_L = false;
        }
      else
        {
          // otherwise, append "e+00"
          number.append("e+00");
        }
    }
  else
    {
      // If 'e' is found, insert a leading zero to be consistent
      // with the Gauss rules in libmesh.  The number of digits
      // used in the exponent is implementation-dependent.  To be
      // totally general, we should check to see how many zeros
      // there are, and only insert exactly as many as we need...
      number.insert(e_pos+2, "0");
    }

  // The L suffix in C/C++ makes the compiler treat a number literal as 'long double'
  if (append_L)
    number.append("L");

  return number;
}
