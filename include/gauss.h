#ifndef __gauss_rule_h__
#define __gauss_rule_h__

#include <cmath>
#include <vector>

#include "gmp.h"
#include "mpfr.h"
#include "gmpfrxx.h"

#include "common_definitions.h" // Real



void gauss_rule(unsigned int n, // Number of points (not order!)
                std::vector<mpfr_class>& x,  // Quadrature points
                std::vector<mpfr_class>& w); // Quadrature weights

#endif
