#include "generator.h"

bool Generator::is_valid() const
{
  if (type == CENTROID && w != mpfr_class(0.))
    return true;

  if (type == MEDIAN && w != mpfr_class(0.) && a != mpfr_class(0.))
    return true;

  if (type == ARBITRARY && w != mpfr_class(0.) && a != mpfr_class(0.) && b != mpfr_class(0.))
    return true;

  if (type == SINGLEPOINT && w != mpfr_class(0.) && a != mpfr_class(0.) && b != mpfr_class(0.))
    return true;

  if (type == RO3 && w != mpfr_class(0.) && a != mpfr_class(0.))
    return true;

  // If we made it here, we're not valid
  std::cout << "Generator type " << type
            << " not valid: weight = " << w
            << ", a = " << a
            << ", b = " << b
            << std::endl;

  return false;
}



void Generator::generate_points_and_weights(std::vector<Point<mpfr_class> > & generated_points,
                                            std::vector<mpfr_class> & generated_weights) const
{
  // Refuse to generate points and weights if the Generator is not valid
  if (!this->is_valid())
    {
      std::cerr << "Error: cannot generate points/weights for invalid Generator!" << std::endl;
      std::abort();
    }

  if (type == CENTROID)
    {
      mpfr_class one_third = mpfr_class(1.)/mpfr_class(3.);
      generated_points.resize(1);
      generated_points[0] = Point<mpfr_class>(one_third, one_third);

      // Note: the clear() call here is important: if the existing
      // vector is larger than the new size, the new w value will not
      // be used!
      generated_weights.clear();
      generated_weights.resize(1, w);
    }
  else if (type == MEDIAN)
    {
      mpfr_class c = mpfr_class(1.) - 2.*a;
      generated_points.resize(3);
      generated_points[0] = Point<mpfr_class>(a, a);
      generated_points[1] = Point<mpfr_class>(a, c);
      generated_points[2] = Point<mpfr_class>(c, a);

      generated_weights.clear();
      generated_weights.resize(3, w);
    }
  else if (type == ARBITRARY)
    {
      mpfr_class c = mpfr_class(1.) - a - b;
      generated_points.resize(6);
      generated_points[0] = Point<mpfr_class>(a, b);
      generated_points[1] = Point<mpfr_class>(b, a);
      generated_points[2] = Point<mpfr_class>(a, c);
      generated_points[3] = Point<mpfr_class>(c, a);
      generated_points[4] = Point<mpfr_class>(b, c);
      generated_points[5] = Point<mpfr_class>(c, b);

      generated_weights.clear();
      generated_weights.resize(6, w);
    }
  else if (type == SINGLEPOINT)
    {
      generated_points.resize(1);
      generated_points[0] = Point<mpfr_class>(a, b);

      generated_weights.clear();
      generated_weights.resize(1, w);
    }
  else if (type == RO3)
    {
      mpfr_class c = mpfr_class(1.) - a - b;
      generated_points.resize(3);
      generated_points[0] = Point<mpfr_class>(a, b);
      generated_points[1] = Point<mpfr_class>(c, a);
      generated_points[2] = Point<mpfr_class>(b, c);

      generated_weights.clear();
      generated_weights.resize(3, w);
    }
  else
    {
      std::cerr << "We'll never get here!" << std::endl;
      std::abort();
    }
}
