#ifndef EXACT_H
#define EXACT_H

#include "gmpfrxx.h"

// Computes the exact integral of the polynomial x^{a} * y^{b}
// Maple commands:
// assume(a, integer, a, nonnegative);
// assume(b, integer, b, nonnegative);
// simplify(int(int(x^a * y^b, y=0..1-x), x=0..1));
// Result is:
// 1/(b+1)*Beta(a+1,b+2)
// where Beta(p,q) := (p-1)! (q-1)! / (p+q-1)! is the beta function, as defined for positive integers.
// After simplification, the result is:
// a! b! / (a + b + 2)!
mpq_class exact_tri(unsigned a, unsigned b);

#endif
