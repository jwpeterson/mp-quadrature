#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <getopt.h> // getopt_long()

#include "gmpfrxx.h"
#include "common_definitions.h"
#include "generator.h"
#include "rule.h"
#include "read_rule.h"
#include "compose_all.h"
#include "exact.h"
#include "dubiner.h"

// Function explaining command line options that gets called when the
// program is run with --help.
void usage();

// Output the generated points and weights (which are assumed to be on
// a reference right triangle) on an equilateral triangle with
// vertices (0,0), (1,0), and (1/2,sqrt(3)/2).  An output filename is
// constructed from input_filename.
void output_equilateral(const std::vector<Point<mpfr_class> > & generated_points,
                        const std::vector<mpfr_class> & generated_weights,
                        const std::string & input_filename);

// Compute the discrete residual function r_N for the current rule
void compute_rN(const std::vector<Point<mpfr_class> > & generated_points,
                const std::vector<mpfr_class> & generated_weights);

// This is code for reading the input/output files used by Zhang's quadrule program
int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  // Filename must now be specified on the command line
  std::string filename = "";

  // options descriptor - this associates several "long" and "short"
  // command line options.  The last element of the longopts array has
  // to be filled with zeros.
  static struct option longopts[] =
    {
      {"input-file",      required_argument, NULL, 'i'},
      {"help",            no_argument,       NULL, 'h'},
      { NULL,             0,                 NULL,  0 }
    };

  // Parse command line options using getopt_long()
  int ch = -1;
  while ((ch = getopt_long(argc, argv, "i:h", longopts, NULL)) != -1)
    {
      switch (ch)
        {
        case 'i':
          {
            if (optarg != NULL)
              filename = std::string(optarg);
            break;
          }

        case 'h':
          usage();
          return 0;

        default:
          usage();
        }
    } // end while

  // Check that the user specified a filename
  if (filename == "")
    {
      std::cerr << "Error, filename was not specified!" << std::endl;
      usage();
      std::abort();
    }

  // Create a Rule object and fill it up from the file
  Rule rule;
  read_rule(filename, rule);

  // Check that the rules are valid and print information about them
  for (unsigned i=0; i<rule.n_generators(); ++i)
    if (!rule[i].is_valid())
      {
        std::cout << "Rule type is: " << rule[i].get_type() << std::endl;
        std::cout << "Error reading parameters for generator " << i << std::endl;
      }

  // Generate the points and weights vectors for this Rule
  std::vector<Point<mpfr_class> > generated_points;
  std::vector<mpfr_class> generated_weights;
  rule.generate_points_and_weights(generated_points, generated_weights);

  // Check that the rules are properly scaled.  If not, we can scale
  // them now...  This may be the case for "foo.in" files, since our
  // foo.out files have all been scaled correctly.
  mpfr_class weight_sum = 0.;
  for (unsigned i=0; i<generated_weights.size(); ++i)
    weight_sum += generated_weights[i];

  if ( abs(weight_sum - mpfr_class(1.)/mpfr_class(2.)) > 1.e-30 )
    {
      // Scale weights by (1/2 * weight_sum)
      for (unsigned i=0; i<generated_weights.size(); ++i)
        generated_weights[i] /= (mpfr_class(2.) * weight_sum);
    }

  // Output points on the equilateral triangle
  // output_equilateral(generated_points, generated_weights, filename);

  // Compute the discrete residual function r_N
  compute_rN(generated_points, generated_weights);

  return 0;
}



void usage()
{
  std::cout << "\n";
  std::cout << "This program reads the generators for 2D quadrature rules on the triangle.\n";
  std::cout << "\n";
  std::cout << "Valid command line options are:\n";
  std::cout << "--input-file, -i filename = Read generators from file 'filename'\n";
  std::cout << "--help, -h                = Print this message.\n";
  std::cout << "\n";
}




void output_equilateral(const std::vector<Point<mpfr_class> > & generated_points,
                        const std::vector<mpfr_class> & generated_weights,
                        const std::string & input_filename)
{
  // Write to a file that is named similarly to the input filename
  std::string output_filename = input_filename;

  // Erase leading path
  size_t slash_pos = output_filename.find_last_of('/');
  if (slash_pos != std::string::npos)
    output_filename.erase(0, slash_pos+1);

  // Erase file extension
  size_t dot_pos = output_filename.find_last_of('.');
  output_filename.erase(dot_pos);
  output_filename += "_equilateral.dat";
  std::cout << "output_filename=" << output_filename << std::endl;

  // Open a stream for output
  std::ofstream out(output_filename.c_str());

  // Set precision flags on the output file
  out.precision(32);
  out.setf(std::ios_base::scientific);

  for (unsigned i=0; i<generated_points.size(); ++i)
    {
      mpfr_class
        xi  = generated_points[i](0),
        eta = generated_points[i](1);

      mpfr_class
        equilateral_x = xi + mpfr_class(0.5)*eta,
        equilateral_y = sqrt(mpfr_class(3.))*0.5*eta;

      out << equilateral_x << ", " << equilateral_y << std::endl;
    }
}



void compute_rN(const std::vector<Point<mpfr_class> > & generated_points,
                const std::vector<mpfr_class> & generated_weights)
{
  // Verify integration of Dubiner polynomials.  Note: they should all
  // be zero except the first one, which should be 1.0
  Dubiner dubiner;

  // The degree of Dubiner polynomials to compute the residual with
  const unsigned max_dubiner_degree = 30;

  for (unsigned dubiner_degree=0; dubiner_degree <= max_dubiner_degree; ++dubiner_degree)
    {
      std::vector<mpfr_class> E, current_vals;

      for (unsigned i=0; i<generated_points.size(); ++i)
        {
          dubiner.p_numeric(dubiner_degree,
                            /*xi=*/ generated_points[i](0),
                            /*eta=*/ generated_points[i](1),
                            current_vals);

          // After the first loop iteration, this resize() should do nothing
          E.resize(current_vals.size());

          for (unsigned j=0; j<current_vals.size(); ++j)
            E[j] += generated_weights[i] * current_vals[j];
        }

      // E[] now contains all of the E(phi_i) except the first
      // entry, which needs to have 1/2 subtracted.  This is because the
      // true values of all the integrals is zero except for the first
      // one.
      E[0] -= 0.5;

      // Build the H1-projection matrix
      Matrix<mpfr_class> A;
      dubiner.build_H1_projection_matrix(dubiner_degree, A);

      // Make sure they are the same size
      if (E.size() != A.n_rows())
        {
          std::cerr << "E vector size does not match matrix size!" << std::endl;
          std::abort();
        }

      // The LU solve will destroy A, therefore make a copy
      // of A first since we also need it to compute the H1-norm of r.
      Matrix<mpfr_class> A_copy(A);

      // Solve Ar = E for r.
      std::vector<mpfr_class> r;
      A.lu_solve(r, E);

      // Compute the H1-norm (squared) of the residual.  This is simply
      // r^T * A * r.
      mpfr_class residual_h1_norm_squared = 0.;
      for (unsigned i=0; i<r.size(); ++i)
        for (unsigned j=0; j<r.size(); ++j)
          residual_h1_norm_squared += r[i] * A_copy(i,j) * r[j];

      // Print norm, squared
      // std::cout << "||r_N||_1^2 = "
      //           << fix_string(residual_h1_norm_squared, /*append_L=*/false) << std::endl;

      // Print norm
      std::cout << std::setw(2) << dubiner_degree
                << ", "
                << fix_string(sqrt(residual_h1_norm_squared), /*append_L=*/false) << std::endl;
    }
}
