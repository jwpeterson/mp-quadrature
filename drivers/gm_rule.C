#include <iomanip>
#include <algorithm> // std::sort
#include <cstdlib> // std::abort
#include <getopt.h> // getopt_long()
#include <cmath> // log10

#include "grundmann_moller.h"

// An unoptimized pow function for mpq_class objects
mpq_class pow(const mpq_class & base, unsigned power)
{
  mpq_class result(1);
  for (unsigned i=0; i<power; ++i)
    result *= base;
  return result;
}

// Formats a rational number with decimal points and long double
// string literals for generating C++ code.
std::string format_rational(const mpq_class & number)
{
  std::string number_string = number.get_str();

  // Find the forward slash, if there is one
  size_t slash_pos = number_string.find('/');

  if (slash_pos != std::string::npos)
    {
      // Insert ".L" before the "/"
      number_string.insert(slash_pos, ".L");

      // And insert a ".L" for the denominator
      number_string += ".L";
    }

  return number_string;
}

// Function explaining command line options that gets called when the
// program is run with --help.
void usage();


// Function which verifies a given quadrature rule defined by the
// points x and weights w.  The third argument is a vector defining
// the order in which to evaluate the rule.
template <class T>
void verify(unsigned rule_index,
            const std::vector<Point<T> >& x,
            const std::vector<T>& w,
            const std::vector<size_t>& indices);


int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // If false, print the points and weights in the order they are
  // computed, if true print them in the weight-sorted order.
  // Eventually make this selectable on the command line.
  bool print_sorted = true;

  // If false, will not compute all the verification integrals.  This
  // can be nice for high-order rules which take an extremely long
  // time to verify.  Eventually make this selectable on the command
  // line.
  bool do_verification = true;

  // Read rule index, s, from command line.  Any integer s>=0 is valid.
  unsigned s = 0;

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  unsigned b = 256;

  // options descriptor - this associates several "long" and "short"
  // command line options.  The last element of the longopts array has
  // to be filled with zeros.
  static struct option longopts[] =
    {
      {"no-verify",      no_argument,       NULL, 'x'},
      {"rule-index",     required_argument, NULL, 's'},
      {"print-unsorted", no_argument,       NULL, 'p'},
      {"binary-digits",  required_argument, NULL, 'b'},
      {"help",           no_argument,       NULL, 'h'},
      { NULL,            0,                 NULL,  0 }
    };

  // Parse command line options using getopt_long()
  int ch = -1;
  while ((ch = getopt_long(argc, argv, "xs:pb:h", longopts, NULL)) != -1)
    {
      switch (ch)
        {
        case 'x':
          do_verification = false;
          break;

        case 's':
          s = atoi(optarg);
          break;

        case 'p':
          print_sorted = false;
          break;

        case 'b':
          b = atoi(optarg);
          break;

        case 'h':
          usage();
          return 0;

        default:
          usage();
        }
    } // end while

  // Set the default precision for real-valued numbers
  // std::cout << "Using " << b << " binary digits for mpfr_class objects." << std::endl;
  mpfr_set_default_prec(b);

  GrundmannMoller gm;
  gm.rule(s);

  // Get a reference to the points and weights arrays
  const std::vector<Point<mpq_class> >& x = gm.get_points();
  const std::vector<mpq_class>& w = gm.get_weights();

  // Initialize the indices array.  This will be used to access into
  // the x and w arrays from now on...
  std::vector<size_t> indices(w.size());
  iota(indices.begin(), indices.end(), 0);

  // If requested, indirect sort the points and weights vectors
  if (print_sorted)
    indirect_sort(w.begin(), w.end(), indices);

  std::cout << "Index " << s
            << " rule has degree " << 2*s+1
            << ", and " << x.size() << " point(s)." << std::endl;

  mpq_class sum_weights = 0;

  // Generate C++ code to resize the _points and _weights vectors appropriately.
  std::cout << "_points.resize(" << x.size() << ");" << std::endl;
  std::cout << "_weights.resize(" << x.size() << ");\n" << std::endl;

  // Determine how many spaces to use to make everything line up in the printout
  unsigned print_width = 1 + static_cast<unsigned>(std::log10(static_cast<double>(x.size())));

  // Pre-compute the widest entry of each column.  We don't have to sort to do this...
  size_t
    widest_x = 0,
    widest_y = 0,
    widest_z = 0,
    widest_w = 0;
  for (unsigned i=0; i<x.size(); ++i)
    {
      widest_x = std::max(widest_x, format_rational(x[i](0)).size());
      widest_y = std::max(widest_y, format_rational(x[i](1)).size());
      widest_z = std::max(widest_z, format_rational(x[i](2)).size());
      widest_w = std::max(widest_w, format_rational(w[i]).size());
    }

  // Print the points and weights
  for (unsigned i=0; i<x.size(); ++i)
    {
      // Choose the sorted or "natural" ordering
      unsigned mapped_index = indices[i];

      // Print points in C++ code that can be used in libmesh
      std::cout << "_points[" << std::setw(print_width) << i << "](0) = " << std::setw(widest_x) << format_rational(x[mapped_index](0)) << ";  "
                << "_points[" << std::setw(print_width) << i << "](1) = " << std::setw(widest_y) << format_rational(x[mapped_index](1)) << ";  "
                << "_points[" << std::setw(print_width) << i << "](2) = " << std::setw(widest_z) << format_rational(x[mapped_index](2)) << ";  "
                << "_weights[" << std::setw(print_width) << i << "] = "   << std::setw(widest_w) << format_rational(w[mapped_index]) << ";"
                << std::endl;

      sum_weights += w[i];
    }

  // std::cout << "Sum of weights = " << sum_weights << std::endl;
  mpq_class error_weights = abs(sum_weights - mpq_class("1/6"));
  std::cout << "Error in sum of weights = " << error_weights << std::endl;

  if (do_verification)
    {
      // We can run the verification for the rational values, but the
      // error is always identically zero (though it's good to verify
      // this).  The more interesting thing is what happens when we
      // convert the rational points and weights to real-valued
      // numbers with a fixed precision and see how the ordering
      // affects the maximum error.
      // std::cout << "\nVerifying rule..." << std::endl;
      // verify(s, x, w, indices);
      // std::cout << "... Verification complete!" << std::endl;

      // Convert the rational point/weight values to finite precision
      // real-valued numbers with b binary digits.  Compare the error in
      // the quadrature approximations for the sorted vs. unsorted rules.
      std::vector<mpfr_class> w_real(w.begin(), w.end());
      std::vector<Point<mpfr_class> > x_real(x.size());
      for (unsigned i=0; i<x_real.size(); ++i)
        for (unsigned j=0; j<3; ++j)
          x_real[i](j) = x[i](j);

      std::cout << "\nVerifying real-valued rule..." << std::endl;
      verify(s, x_real, w_real, indices);
      std::cout << "... Verification complete!" << std::endl;
    }

  return 0;
}



void usage()
{
  std::cout << "\n";
  std::cout << "This program generates the points and weights for a 3D Grundmann-Moller quadrature rule.\n";
  std::cout << "\n";
  std::cout << "Valid command line options are:\n";
  std::cout << "--rule-index, -s #    = Build GM rule with index s, degree 2*s+1.\n";
  std::cout << "--no-verify, -x       = Skip verifying that the rule exactly integrates polynomials with degree<=2*s+1.\n";
  std::cout << "--print-unsorted, -p  = Print the points and weights in unsorted order.\n";
  std::cout << "--binary-digits, -b # = Default number of binary digits to use for mpfr_class variables.\n";
  std::cout << "--help, -h            = Print this message.\n";
  std::cout << "\n";
}



template <class T>
void verify(unsigned rule_index,
            const std::vector<Point<T> >& x,
            const std::vector<T>& w,
            const std::vector<size_t>& indices)
{
  // Keep track of the max error encountered
  T max_error = 0;

  // Mostly copied from conical_product_3D.C, the verification code
  // for tets should probably be factored out into a common
  // function...
  unsigned max_order = 2*rule_index + 1;
  for (unsigned x_power=0; x_power<max_order; ++x_power)
    for (unsigned y_power=0; y_power<max_order; ++y_power)
      for (unsigned z_power=0; z_power<max_order; ++z_power)
        {
          // Only try to integrate polynomials we can integrate exactly
          if (x_power + y_power + z_power > max_order)
            continue;

          // Compute the quadrature result
          T sum = 0;
          for (unsigned n_qp=0; n_qp<x.size(); ++n_qp)
            {
              unsigned mapped_index = indices[n_qp];

              sum +=
                w[mapped_index] *
                pow(x[mapped_index](0), x_power) *
                pow(x[mapped_index](1), y_power) *
                pow(x[mapped_index](2), z_power);
            }

          // std::cout << "quadrature = " << sum << std::endl;

          // Compute the exact solution as described above, divide
          // term-by-term to avoid huge numbers/overflows (even though
          // this is GMP).
          T analytical = 1;
          {
            // Sort the a, b, c values
            unsigned sorted_powers[3] = {x_power, y_power, z_power};
            std::sort(sorted_powers, sorted_powers+3);

            // Cancel the largest power with the denominator, fill in the
            // entries for the remaining numerator terms and the denominator.
            std::vector<unsigned>
              numerator_1(sorted_powers[0] > 1 ? sorted_powers[0]-1 : 0),
              numerator_2(sorted_powers[1] > 1 ? sorted_powers[1]-1 : 0),
              denominator(3 + sorted_powers[0] + sorted_powers[1]);

            // Fill up the vectors with sequences starting at the right values.
            iota(numerator_1.begin(), numerator_1.end(), 2);
            iota(numerator_2.begin(), numerator_2.end(), 2);
            iota(denominator.begin(), denominator.end(), sorted_powers[2]+1);

            // The denominator is guaranteed to have the most terms...
            for (unsigned i=0; i<denominator.size(); ++i)
              {
                if (i < numerator_1.size())
                  analytical *= numerator_1[i];

                if (i < numerator_2.size())
                  analytical *= numerator_2[i];

                analytical /= denominator[i];
              }
          }

          // std::cout << "analytical = " << analytical << std::endl;

          // Compute the absolute error:
          T abs_err = abs(sum-analytical);

          max_error = max(abs_err, max_error);

          // Print message.  In 3D, this is just way too much output.
          // std::cout << "Computing integral of: x^" << x_power << " y^" << y_power << " z^" << z_power
          //           << ", abs_err = " << abs_err << std::endl;

          // Abort if error is too large.
          // if (abs_err > T(1.e-30))
          //   {
          //     std::cerr << "Quadrature error too large, possible problem with points and weights!" << std::endl;
          //     std::abort();
          //   }
        } // end triple-nested for loop

  std::cout << "max_error = " << max_error << std::endl;
}
