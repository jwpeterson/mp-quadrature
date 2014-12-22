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
  Conical conical;

  const unsigned max_dubiner_degree = 4;

  // 1.) Compute integral(phi(i)) for each of the Dubiner
  // polynomials.  They should all be zero except the first one,
  // since they are all orthogonal to constants.
  for (unsigned dubiner_degree=/*0*/max_dubiner_degree; dubiner_degree <= max_dubiner_degree; ++dubiner_degree)
    {
      // The number of polynomials in the Dubiner basis of degree 'dubiner_degree'
      const unsigned Np = (dubiner_degree+1)*(dubiner_degree+2)/2;

      // The points and weights from a conical product quadrature rule
      // of high enough order to integrate a Dubiner polynomial of
      // degree dubiner_degree.
      conical.rule2D(dubiner_degree);
      const std::vector<Point<mpfr_class> > & conical_rule_points = conical.get_points();
      const std::vector<mpfr_class> & conical_rule_weights = conical.get_weights();

      // Nq is the number of points in the quadrature rule
      const unsigned Nq = conical_rule_points.size();

      // // Verification: print the Conical product rules and weights
      // {
      //   std::cout << "\n";
      //   for (unsigned int i=0; i<Nq; ++i)
      //     {
      //       std::cout << "_points[" << std::setw(2) << i << "](0)=" << fix_string(conical_rule_points[i](0)) << ";\n";
      //       std::cout << "_points[" << std::setw(2) << i << "](1)=" << fix_string(conical_rule_points[i](1)) << ";\n";
      //     }
      //   for (unsigned int i=0; i<Nq; ++i)
      //     std::cout << "_weights[" << std::setw(2) << i << "]=" << fix_string(conical_rule_weights[i]) << ";\n";
      // }

      std::vector<mpfr_class> quadrature_sums(Np);
      for (unsigned q=0; q<Nq; ++q)
        {
          // Evaluate all the Dubiner polynomials at the current qp
          // using the "numerical" Dubiner polynomials.
          std::vector<mpfr_class> current_vals;
          dubiner.p_numeric(dubiner_degree,
                            /*xi=*/  conical_rule_points[q](0),
                            /*eta=*/ conical_rule_points[q](1),
                            current_vals);

          // Evaluate all the Dubiner polynomials at the current qp
          // using the "symbolic" Dubiner polynomials.
          // dubiner.p(dubiner_degree,
          //           /*xi=*/  conical_rule_points[q](0),
          //           /*eta=*/ conical_rule_points[q](1),
          //           current_vals);

          // Compute the phi(i) integrals using conical product rule quadrature.
          for (unsigned i=0; i<Np; ++i)
            quadrature_sums[i] += conical_rule_weights[q] * current_vals[i];
        }

      // Print the results
      std::cout << "\n int phi(i)*1" << std::endl;
      for (unsigned i=0; i<Np; ++i)
        std::cout << quadrature_sums[i] << std::endl;
    }



  // 2.) Compute integral(phi(i)**2) for each of the Dubiner
  // polynomials.  These are the "orthogonality coefficients"...
  for (unsigned dubiner_degree=/*0*/max_dubiner_degree; dubiner_degree <= max_dubiner_degree; ++dubiner_degree)
    {
      // The number of polynomials in the Dubiner basis of degree 'dubiner_degree'
      const unsigned Np = (dubiner_degree+1)*(dubiner_degree+2)/2;

      // The points and weights from a conical product quadrature rule
      // of high enough order to integrate phi(i)**2, which has degree 2*dubiner_degree
      conical.rule2D(2*dubiner_degree);
      const std::vector<Point<mpfr_class> > & conical_rule_points = conical.get_points();
      const std::vector<mpfr_class> & conical_rule_weights = conical.get_weights();

      // Nq is the number of points in the quadrature rule
      const unsigned Nq = conical_rule_points.size();

      std::vector<mpfr_class> quadrature_sums(Np);
      for (unsigned q=0; q<Nq; ++q)
        {
          std::vector<mpfr_class> current_vals;

          // Evaluate all the Dubiner polynomials at the current qp
          // using the "numerical" Dubiner polynomials.
          // dubiner.p_numeric(dubiner_degree,
          //                   /*xi=*/  conical_rule_points[q](0),
          //                   /*eta=*/ conical_rule_points[q](1),
          //                   current_vals);

          // Evaluate all the Dubiner polynomials at the current qp
          // using the "symbolic" Dubiner polynomials.
          dubiner.p(dubiner_degree,
                    /*xi=*/  conical_rule_points[q](0),
                    /*eta=*/ conical_rule_points[q](1),
                    current_vals);

          // Square each of the current_vals to represent phi(i)**2
          for (unsigned i=0; i<Np; ++i)
            current_vals[i] = current_vals[i]*current_vals[i];

          // Accumulate the phi(i)**2 integrals using conical product rule quadrature.
          for (unsigned i=0; i<Np; ++i)
            quadrature_sums[i] += conical_rule_weights[q] * current_vals[i];
        }

      // Print the results
      std::cout << "\n int phi(i)**2" << std::endl;
      for (unsigned i=0; i<Np; ++i)
        std::cout << quadrature_sums[i] << std::endl;
    }

  return 0;
}
