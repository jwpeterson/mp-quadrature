#ifndef GENERATOR_H
#define GENERATOR_H

#include "common_definitions.h"
#include "gmpfrxx.h"

/**
 * Simple object that stores all the information about a Generator.
 */
class Generator
{
public:
  // There are three different types of generators
  enum GeneratorType
    {
      CENTROID  =  0,
      MEDIAN    =  1,
      ARBITRARY =  2,
      INVALID   = 99
    };

  Generator(GeneratorType type_in = INVALID,
            const mpfr_class & w_in = mpfr_class(0.),
            const mpfr_class & a_in = mpfr_class(0.),
            const mpfr_class & b_in = mpfr_class(0.)) :
    type(type_in),
    w(w_in),
    a(a_in),
    b(b_in)
  {}

  // A generator is defined to be valid when its required parameters
  // are non-zero
  bool is_valid()
  {
    if (type == CENTROID && w != mpfr_class(0.))
      return true;

    if (type == MEDIAN && w != mpfr_class(0.) && a != mpfr_class(0.))
      return true;

    if (type == ARBITRARY && w != mpfr_class(0.) && a != mpfr_class(0.) && b != mpfr_class(0.))
      return true;

    // If we made it here, we're not valid
    std::cout << "Generator type " << type << " not valid: weight = " << w << ", a = " << a << ", b = " << b << std::endl;
    return false;
  }

  virtual ~Generator() {}

  GeneratorType & get_type() { return type; }
  mpfr_class & get_w() { return w; }
  mpfr_class & get_a() { return a; }
  mpfr_class & get_b() { return b; }

  // Generate the points and weights associated with this Generator
  // and store them in the passed-in vectors.
  void generate_points_and_weights(std::vector<Point<mpfr_class> > & generated_points,
                                   std::vector<mpfr_class> & generated_weights);

protected:
  GeneratorType type;
  mpfr_class w;
  mpfr_class a;
  mpfr_class b;
};

#endif
