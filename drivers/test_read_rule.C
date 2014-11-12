#include <iostream>
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


// Function explaining command line options that gets called when the
// program is run with --help.
void usage();



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
    {
      if (rule[i].is_valid())
        {
          std::cout << "w[gen=" << i << "]=" << rule[i].get_w()
                    << ", a[gen=" << i << "]=" << rule[i].get_a()
                    << ", b[gen=" << i << "]=" << rule[i].get_b()
                    << std::endl;
        }
      else
        {
          std::cout << "Rule type is: " << rule[i].get_type() << std::endl;
          std::cout << "Error reading parameters for generator " << i << std::endl;
        }
    }

  // Test generating points and weights vectors for this Rule
  std::vector<Point<mpfr_class> > generated_points;
  std::vector<mpfr_class> generated_weights;
  rule.generate_points_and_weights(generated_points, generated_weights);
  for (unsigned i=0; i<generated_points.size(); ++i)
    {
      std::cout << "Point " << i
                << ": (" << generated_points[i](0)
                << ", " << generated_points[i](1)
                << "), weight: "
                << generated_weights[i]
                << std::endl;
    }


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


