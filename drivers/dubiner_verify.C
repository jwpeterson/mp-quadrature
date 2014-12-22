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

  enum DubinerImplementation {SYMBOLIC=0, NUMERIC};

  // Pick an implementation to verify
  DubinerImplementation dubiner_implementation = NUMERIC;

  // Compute integral(phi(i)*phi(j)) for each of the Dubiner polynomials.
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

      // There are Np^2 quadrature sums to compute in this case
      Matrix<mpfr_class> quadrature_sums(Np, Np);
      for (unsigned q=0; q<Nq; ++q)
        {
          std::vector<mpfr_class> current_vals;

          // Evaluate all the Dubiner polynomials at the current qp
          // using the "numerical" Dubiner polynomials.
          if (dubiner_implementation == NUMERIC)
            dubiner.p_numeric(dubiner_degree,
                              /*xi=*/  conical_rule_points[q](0),
                              /*eta=*/ conical_rule_points[q](1),
                              current_vals);

          // Evaluate all the Dubiner polynomials at the current qp
          // using the "symbolic" Dubiner polynomials.
          else
            dubiner.p(dubiner_degree,
                      /*xi=*/  conical_rule_points[q](0),
                      /*eta=*/ conical_rule_points[q](1),
                      current_vals);

          if (current_vals.size() != Np)
            {
              std::cerr << "current_vals.size() = " << current_vals.size() << " does not match Np = " << Np << std::endl;
              std::abort();
            }

          // std::cout << "computing mass matrix entries." << std::endl;

          // Make a matrix of the phi(i) * phi(j) values
          Matrix<mpfr_class> mass(Np, Np);
          for (unsigned i=0; i<Np; ++i)
            for (unsigned j=0; j<Np; ++j)
              {
                // std::cout << i << "," << j << std::endl;
                mass(i,j) = current_vals[i]*current_vals[j];
              }

          // std::cout << "accumulating quadrature sums." << std::endl;

          // Accumulate the phi(i)*phi(j) integrals using conical product rule quadrature.
          for (unsigned i=0; i<Np; ++i)
            for (unsigned j=0; j<Np; ++j)
              quadrature_sums(i,j) += conical_rule_weights[q] * mass(i,j);
        }

      // Print the results
      std::cout << "\n int phi(i)*phi(j)" << std::endl;
      quadrature_sums.print();

      // Print an error message if an off-diagonal entry is too large
     for (unsigned i=0; i<Np; ++i)
       {
         for (unsigned j=0; j<Np; ++j)
           {
             if ((i != j) && (abs(quadrature_sums(i,j)) > 1.e-30))
               {
                 std::cerr << "Matrix entry " << i << "," << j << " should be zero, but is " << quadrature_sums(i,j) << std::endl;
                 std::abort();
               }
           }
       }
    }

  return 0;
}
