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



// Maple has a lot of trouble computing this integral exactly. Try for yourself:
// assume(a, integer, a, nonnegative);
// assume(b, integer, b, nonnegative);
// assume(c, integer, c, nonnegative);
// int(int(int(x^a * y^b * z^c, z=0..1-x-y), y=0..1-x), x=0..1);
//
// The main trick is to perform a change of variables to get the
// second integral into the form of the beta function integral...
//
// 1.) Integrate in z: result is I = int(int((x^a * y^b * (1-x-y)^(c+1))/(c+1), y=0..1-x), x=0..1)
// 2.) Introduce change of variables eta = y/(1-x)
//     Then I = Beta(b+1, c+2) / (c+1) * int(x^a * (1-x)^(b+c+2), x=0..1)
// 3.) I can now be computed using the Beta function integral again to get:
//     I = Beta(b+1, c+2) / (c+1) * Beta(a+1, b+c+3)
// 4.) Substituting in the Beta function definition,
//     Beta(p,q) := (p-1)! (q-1)! / (p+q-1)!,
//     and canceling terms, we finally obtain:
//     I = a! b! c! / (a + b + c + 3)!
mpq_class exact_tet(unsigned a, unsigned b, unsigned c);



// The exact integral of x^p over the Pyramid reference element.
mpq_class exact_pyr(unsigned p);

#endif
