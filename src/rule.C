#include "rule.h"
#include "exact.h"

unsigned Rule::n_generators() const
{
  return generators.size();
}



void Rule::push_back(const Generator & generator)
{
  generators.push_back(generator);
}



// Writeable reference to the ith generator, if it exists
Generator& Rule::operator[](unsigned i)
{
  if (i >= generators.size())
    {
      std::cerr << "Invalid index " << i << " requested!" << std::endl;
      std::abort();
    }

  return generators[i];
}



void Rule::generate_points_and_weights(std::vector<Point<mpfr_class> > & generated_points,
                                       std::vector<mpfr_class> & generated_weights) const
{
  // Vectors for individual generators
  std::vector<Point<mpfr_class> > single_generator_points;
  std::vector<mpfr_class> single_generator_weights;

  for (unsigned i=0; i<generators.size(); ++i)
    {
      generators[i].generate_points_and_weights(single_generator_points,
                                                single_generator_weights);

      // Append the points and weights to the end of the user's vectors
      generated_points.insert(/*pos=*/generated_points.end(),
                              /*first=*/single_generator_points.begin(),
                              /*last=*/single_generator_points.end());

      generated_weights.insert(/*pos=*/generated_weights.end(),
                               /*first=*/single_generator_weights.begin(),
                               /*last=*/single_generator_weights.end());
    }
}



bool Rule::verify(unsigned int total_degree) const
{
  // Generate the points and weights
  std::vector<Point<mpfr_class> > generated_points;
  std::vector<mpfr_class> generated_weights;
  this->generate_points_and_weights(generated_points, generated_weights);

  // Originally from conical_product_2D.C
  for (unsigned x_power=0; x_power<total_degree; ++x_power)
    for (unsigned y_power=0; y_power<total_degree; ++y_power)
      {
        // Only try to integrate polynomials we can integrate exactly
        if (x_power + y_power > total_degree)
          continue;

        // Conical product rule uses 0-based array access!
        mpfr_class sum = 0.;
        for (unsigned n_qp=0; n_qp<generated_weights.size(); ++n_qp)
          sum += generated_weights[n_qp] *
            pow(generated_points[n_qp](0), x_power) *
            pow(generated_points[n_qp](1), y_power);

        // std::cout << "quadrature = " << sum << std::endl;

        // Compute the analytical integral value.
        mpfr_class analytical = exact_tri(x_power, y_power);

        // std::cout << "analytical = " << analytical << std::endl;

        // Compute the absolute error:
        mpfr_class abs_err = abs(sum-analytical);

        // Print message
        std::cout << "Computing integral of: x^" << x_power << " y^" << y_power
                  << ", abs_err = " << abs_err << std::endl;

        // Abort if error is too large.
        if (abs_err > mpfr_class(1.e-30))
          {
            // std::cerr << "Quadrature error too large, "
            //           << "possible problem with points and weights!"
            //           << std::endl;

            // std::abort();
            return false;
          }
      } // end for (x_power, y_power)

  // If we made it here, we're verified.
  return true;
}
