// nlopt includes. The preprocessor define is set on the compile line
// with a -D flag.
#ifndef HAVE_NLOPT

int main()
{
  std::cout << "This driver requires nlopt. Install nlopt and then rerun configure, "
            << "using the --with-nlopt-include and --with-nlopt-lib flags to specify "
            << "its location." << std::endl;
  return 0;
}

#else

// mp-quadrature includes
#include "solver_data.h"
#include "newton.h"
#include "gradient_descent.h"
#include "nlcg.h"
#include "ro3.h"

// nlopt includes
#include <nlopt.h>

// C/C++ includes
#include <iostream>
#include <math.h>
#include <stdlib.h> // random()
#include <getopt.h> // getopt_long()

// Objective function f : R^n -> R and, if grad != nullptr, grad(f).
double myfunc(unsigned n, const double * x, double * grad, void * my_func_data);

// Constraints are of the form fc(x) <= 0, where fc is required to have
// the same parameters as the objective function. In our case, there is
// one inequality constraint associated with each ro3 orbit of the form:
// 1 - x - y >= 0,
// which we multiply by -1 to write as:
// x + y - 1 <= 0
// So in this case, the incoming "data" pointer just needs to tell us
// the index of the x dof for the orbit being constrained.
double myconstraint(unsigned n, const double *x, double *grad, void *data);

// Convert a numerical nlopt_result successful return code to a human readable string.
std::string nlopt_result_to_string(nlopt_result res);

// Call this function when the program is run with the wrong arguments.
void usage();

int main(int argc, char ** argv)
{
  // Convenience: Print exact integral values for some polynomials before returning.
  // std::cout << "exact_00 = " << exact_tri(0,0) << std::endl; // 1/2
  // std::cout << "exact_20 = " << exact_tri(2,0) << std::endl; // 1/12
  // std::cout << "exact_30 = " << exact_tri(3,0) << std::endl; // 1/20
  // std::cout << "exact_21 = " << exact_tri(2,1) << std::endl; // 1/60
  // std::cout << "exact_40 = " << exact_tri(4,0) << std::endl; // 1/30
  // std::cout << "exact_50 = " << exact_tri(5,0) << std::endl; // 1/42
  // std::cout << "exact_41 = " << exact_tri(4,1) << std::endl; // 1/210
  // return 0;

  // You can't trust all these digits from doubles, but you can with
  // mpfr_class objects.
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // Process command line arguments
  static struct option longopts[] =
    {
      {"degree",      required_argument, NULL, 'd'},
      {"n-centroid",  required_argument, NULL, 'c'},
      {"n-vertex",    required_argument, NULL, 'v'},
      {"n-edge",      required_argument, NULL, 'e'},
      {"n-median",    required_argument, NULL, 'm'},
      {"n-general",   required_argument, NULL, 'g'},
      {"rule-string", required_argument, NULL, 'r'},
      {"help",        no_argument,       NULL, 'h'},
      { NULL,         0,                 NULL,  0 }
    };

  // To be set from the command line
  unsigned int d = 0;

  // Number of centroid, vertex, edge, median, and general orbits.
  unsigned int nc = 0;
  unsigned int nv = 0;
  unsigned int ne = 0;
  unsigned int nm = 0;
  unsigned int ng = 0;

  // A string in the format "nc,nv,ne,nm,ng" specifying the number
  // of each kind of orbit. There must be exactly 5 comma-separated
  // values or else an error is thrown.
  std::string rule_string = "";

  // A colon following an option means it has an argument.
  // If there's an unrecognized argument, getopt_long()
  // prints a message, so we don't handle it in the cases below.
  int ch = -1;
  while ((ch = getopt_long(argc, argv, "hd:c:v:e:m:g:r:", longopts, NULL)) != -1)
    {
      switch (ch)
        {
        case 'd': { d = atoi(optarg); break; }
        case 'c': { nc = atoi(optarg); break; }
        case 'v': { nv = atoi(optarg); break; }
        case 'e': { ne = atoi(optarg); break; }
        case 'm': { nm = atoi(optarg); break; }
        case 'g': { ng = atoi(optarg); break; }
        case 'r': { if (optarg) rule_string = std::string(optarg); break; }
        case 'h': { usage(); return 0; }
        }
    } // end while

  if (d == 0 || nc > 1 || nv > 1)
    {
      std::cout << "Error, invalid command line arguments provided!" << std::endl;
      usage();
      return 1;
    }

  // Parse the rule string if present
  if (!rule_string.empty())
    {
      // Debugging
      // std::cout << "rule_string = " << rule_string << std::endl;

      std::stringstream rule_ss(rule_string);
      std::string val;
      std::vector<unsigned int> vals;
      while (std::getline(rule_ss, val, ','))
        {
          char * endptr;
          unsigned int tmp =
            std::strtol(val.c_str(), &endptr, /*base=*/10);

          if (val.c_str() != endptr)
            vals.push_back(tmp);
        }

      // Debugging
      // std::cout << "rule numerals = " << std::endl;
      // for (unsigned int i=0; i<vals.size(); ++i)
      //   std::cout << vals[i] << std::endl;

      // Make sure exactly the right number of values were parsed.
      if (vals.size() != 5)
        {
          std::cout << "Error, invalid command line arguments provided!" << std::endl;
          usage();
          return 1;
        }

      // Assign values, overwriting anything that may have already
      // been there.
      nc = vals[0];
      nv = vals[1];
      ne = vals[2];
      nm = vals[3];
      ng = vals[4];
    }

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  // It's faster to comment this out while debugging.
  mpfr_set_default_prec(256);

  // Use the current time since epoch as a seed.
  time_t seed = time(nullptr);
  // long seed = 1567251384;
  std::cout << "seed=" << seed << std::endl;
  srandom(seed);

  // Build ro3 rule object
  Ro3 r(d, nc, nv, ne, nm, ng);

  // d==2, dim==2. Best known PI rule has 3 QP.
  // Note: This case is somewhat tricky despite also being the
  // simplest non-trivial case!  The issue is that the Jacobian is
  // singular at the root, so Newton iterations converge very slowly,
  // and you don't obtain the full requested accuracy in the final
  // solution.  I am a bit surprised that lu_solve() does not throw a
  // singular matrix exception, I think the problem is that the matrix
  // is only approximately singular, and we only throw if it's
  // *exactly* singular.
  // -d2 -c0 -v0 -e1 -g0 # 3 QP <-- Solution confirmed

  // d==3, dim==4. Best known PI rule is a conical product rule with 4 QP.
  // -d3 -c1 -v0 -e0 -m0 -g1 # 4 QP <-- No solution (-ve wt soln only)
  // -d3 -c0 -v0 -e2 -m0 -g0 # 6 QP <-- No solution (I think we can prove this.)
  // -d3 -c0 -v0 -e1 -m1 -g0 # 6 QP <-- ???
  // -d3 -c0 -v0 -e0 -m2 -g0 # 6 QP <-- _Many_ solutions found?! (We can show the possibility of infinitely many solutions.)
  // -d3 -c1 -v1 -e1 -m0 -g0 # 7 QP <-- New (?) solution
  // -d3 -c1 -v1 -e0 -m1 -g0 # 7 QP <-- _Many_ solutions found

  // d==4, dim==5
  // o The best known deg=4 PI rule is a D3-invariant rule with 6 QP.
  // o We also found a degree=4 rule with 6 QPs having ne=ng=1.
  // o For the case with nc=1, ne=2 there does not seem to be a solution, but
  //   the same (global) minimum objective function value of ~2.57337-08
  // is found by the algorithm regardless of starting point.
  // -d4 -c0 -v0 -e0 -m1 -g1 # 6 QP <-- No solution (Should be D3-invariant in libmesh)
  // -d4 -c0 -v0 -e1 -m0 -g1 # 6 QP <-- New (?) soln found!
  // -d4 -c1 -v1 -e0 -m0 -g1 # 7 QP <-- No solution
  // -d4 -c1 -v0 -e2 -m0 -g0 # 7 QP <-- No solution
  // -d4 -c0 -v1 -e2 -m0 -g0 # 9 QP <-- No solution

  // d==5, dim==7. Best known PI rule has 7 QP.
  // -d5 -c1 -v0 -e0 -m0 -g2 # 7 QP <-- Also a D3 rule in libmesh
  // -d5 -c1 -v0 -e0 -m3 -g0 # 10 QP <-- A D3 rule

  // d=6, dim==10. Best known PI rule has 12 QP.
  // -d6 -c1 -v0 -e0 -m0 -g3 # 10 QP <-- No solution
  // -d6 -c0 -v1 -e0 -m0 -g3 # 12 QP <-- No solution
  // -d6 -c0 -v0 -e2 -m0 -g2 # 12 QP <-- New (?) solution
  // -d6 -c1 -v1 -e1 -m0 -g2 # 13 QP <-- No solution
  // -d6 -c1 -v0 -e3 -m0 -g1 # 13 QP <-- No solution
  // -d6 -c0 -v0 -e5 -m0 -g0 # 15 QP <-- No solution
  // -d6 -c1 -v1 -e4 -m0 -g0 # 16 QP <-- No solution

  // d==7, dim=12. Best known PI rule is Gatermann's rule.
  // -d7 -c0 -v0 -e0 -m0 -g4 # 12 QP <-- Gatermann solution
  // -d7 -c1 -v0 -e1 -m0 -g3 # 13 QP <-- New (?) solution
  // -d7 -c0 -v1 -e1 -m0 -g3 # 15 QP <-- No solution
  // -d7 -c0 -v0 -e3 -m0 -g2 # 15 QP <-- New (?) solution
  // -d7 -c1 -v1 -e2 -m0 -g2 # 16 QP <-- New (?) solution
  // -d7 -c1 -v0 -e4 -m0 -g1 # 16 QP <-- No solution
  // -d7 -c0 -v1 -e4 -m0 -g1 # 18 QP <-- No solution
  // -d7 -c1 -v1 -e5 -m0 -g0 # 18 QP <-- No solution

  // d==8, dim=15, best PI degree 8 rule in libmesh has 16 QPs
  // -d8 -c0 -v0 -e0 -m0 -g5 # 15 QP <-- No solution
  // -d8 -c1 -v0 -e1 -m0 -g4 # 16 QP <-- New (?) solution
  // -d8 -c0 -v1 -e1 -m0 -g4 # 18 QP <-- New (?) solution
  // -d8 -c0 -v0 -e3 -m0 -g3 # 18 QP <-- No solution
  // -d8 -c1 -v1 -e2 -m0 -g3 # 19 QP <-- New (?) solution
  // -d8 -c1 -v0 -e4 -m0 -g2 # 19 QP <-- No solution
  // -d8 -c0 -v1 -e4 -m0 -g2 # 21 QP <-- No solution
  // -d8 -c0 -v0 -e6 -m0 -g1 # 21 QP <-- No solution
  // -d8 -c1 -v1 -e5 -m0 -g1 # 22 QP <-- No solution
  // -d8 -c1 -v0 -e7 -m0 -g0 # 22 QP <-- No solution
  // -d8 -c0 -v1 -e7 -m0 -g0 # 24 QP <-- No solution

  // d==9, dim=19, best PI degree 9 rule in libmesh has 19 QPs
  // -d9 -c1 -v0 -e0 -m0 -g6 # 19 QP <-- One solution found, same as D3!
  // -d9 -c0 -v1 -e0 -m0 -g6 # 21 QP <-- New (?) solution
  // -d9 -c0 -v0 -e2 -m0 -g5 # 21 QP <-- No solution
  // -d9 -c1 -v1 -e1 -m0 -g5 # 22 QP <-- New (?) solution
  // -d9 -c1 -v0 -e3 -m0 -g4 # 22 QP <-- No solution
  // -d9 -c0 -v1 -e3 -m0 -g4 # 24 QP <-- New (?) solution
  // -d9 -c0 -v0 -e5 -m0 -g3 # 24 QP <-- No solution
  // -d9 -c1 -v1 -e4 -m0 -g3 # 25 QP <-- No solution
  // -d9 -c1 -v0 -e6 -m0 -g2 # 25 QP <-- No solution
  // -d9 -c0 -v1 -e6 -m0 -g2 # 27 QP <-- No solution
  // -d9 -c0 -v0 -e8 -m0 -g1 # 27 QP <-- No solution
  // -d9 -c1 -v1 -e7 -m0 -g1 # 28 QP <-- No solution
  // -d9 -c1 -v0 -e9 -m0 -g0 # 28 QP <-- No solution
  // -d9 -c0 -v1 -e9 -m0 -g0 # 30 QP <-- No solution

  // d==10, dim=22, best PI degree 10 rule in libmesh has 25 QPs
  // -d10 -c1 -v0 -e0 -m0 -g7 # 22 QP <-- No solution
  // -d10 -c0 -v1 -e0 -m0 -g7 # 24 QP <-- No solution
  // -d10 -c0 -v0 -e2 -m0 -g6 # 24 QP <-- SIX New (?) solutions
  // -d10 -c1 -v1 -e1 -m0 -g6 # 25 QP <-- New (?) solution
  // -d10 -c1 -v0 -e3 -m0 -g5 # 25 QP <-- TWO New (?) solutions
  // -d10 -c0 -v1 -e3 -m0 -g5 # 27 QP <-- New (?) solution
  // -d10 -c0 -v0 -e5 -m0 -g4 # 27 QP <-- No solution
  // -d10 -c1 -v1 -e4 -m0 -g4 # 28 QP <-- No solutions found
  // -d10 -c1 -v0 -e6 -m0 -g3 # 28 QP
  // -d10 -c1 -v0 -e7 -m0 -g2 # 28 QP
  // -d10 -c0 -v1 -e7 -m0 -g2 # 30 QP
  // -d10 -c0 -v0 -e8 -m0 -g2 # 30 QP
  // -d10 -c1 -v0 -e9 -m0 -g1 # 31 QP
  // -d10 -c0 -v1 -e9 -m0 -g1 # 33 QP
  // -d10 -c1 -v1 -e10 -m0 -g0 # 34 QP

  // d==11, dim=26, best PI degree 11 rule in libmesh has 30 QPs
  // -d11 -c0 -v0 -e1 -m0 -g8 # 27 QP <-- Two New solutions
  // -d11 -c1 -v1 -e0 -m0 -g8 # 28 QP <-- No solutions found
  // -d11 -c1 -v0 -e2 -m0 -g7 # 28 QP <-- One New solution
  // -d11 -c0 -v1 -e2 -m0 -g7 # 30 QP <-- One New solution
  // -d11 -c0 -v0 -e4 -m0 -g6 # 30 QP <-- One New solution
  // -d11 -c1 -v1 -e3 -m0 -g6 # 31 QP <-- Two New solutions
  // -d11 -c1 -v0 -e5 -m0 -g5 # 31 QP <-- No solutions found
  // -d11 -c0 -v1 -e5 -m0 -g5 # 33 QP <-- No solutions found
  // -d11 -c0 -v0 -e7 -m0 -g4 # 33 QP <-- No solutions found
  // -d11 -c1 -v1 -e6 -m0 -g4 # 34 QP <-- No solutions found
  // -d11 -c1 -v0 -e8 -m0 -g3 # 34 QP
  // -d11 -c0 -v1 -e8 -m0 -g3 # 36 QP
  // -d11 -c0 -v0 -e10 -m0 -g2 # 36 QP
  // -d11 -c1 -v1 -e9  -m0 -g2 # 37 QP
  // -d11 -c1 -v0 -e11 -m0 -g1 # 37 QP
  // -d11 -c0 -v1 -e11 -m0 -g1 # 39 QP
  // -d11 -c0 -v0 -e13 -m0 -g0 # 39 QP
  // -d11 -c1 -v1 -e12 -m0 -g0 # 40 QP

  // d==12, dim=31, best PI degree 12 rule in libmesh has 33 QPs
  // -d12 -c1 -v0 -e0 -m0 -g10 # 31 QP <-- No solutions found
  // -d12 -c0 -v1 -e0 -m0 -g10 # 33 QP <-- No solutions found
  // -d12 -c0 -v0 -e2 -m0 -g9 # 33 QP <-- No solutions found
  // -d12 -c0 -v0 -e1 -m1 -g9 # 33 QP <-- No solutions found
  // -d12 -c0 -v0 -e0 -m2 -g9 # 33 QP <-- No solutions found
  // -d12 -c1 -v1 -e1 -m0 -g9 # 34 QP <-- No solutions found
  // -d12 -c1 -v0 -e3 -m0 -g8 # 34 QP <-- No solutions found
  // -d12 -c0 -v1 -e3 -m0 -g8 # 36 QP <-- No solutions found
  // -d12 -c0 -v0 -e5 -m0 -g7 # 36 QP <-- No solutions found
  // -d12 -c1 -v1 -e4 -m0 -g7 # 37 QP <-- No solutions found
  // ...

  // d==13, dim=35, best PI degree 13 rule in libmesh has 37 QPs
  // -d13 -c0 -v0 -e1 -m0 -g11 # 36 QP <-- No solutions found
  // -d13 -c1 -v1 -e0 -m0 -g11 # 37 QP <-- No solutions found
  // -d13 -c1 -v0 -e2 -m0 -g10 # 37 QP <-- No solutions found
  // -d13 -c1 -v0 -e1 -m1 -g10 # 37 QP <-- No solutions found
  // -d13 -c1 -v0 -e0 -m2 -g10 # 37 QP <-- No solutions found
  // -d13 -c0 -v1 -e2 -m0 -g10 # 39 QP
  // -d13 -c0 -v0 -e4 -m0 -g9 # 39 QP
  // -d13 -c1 -v1 -e3 -m0 -g9 # 40 QP
  // -d13 -c1 -v0 -e5 -m0 -g8 # 40 QP
  // -d13 -c0 -v1 -e5 -m0 -g8 # 42 QP
  // ...

  // d==14, dim=40, best PI degree 14 rule in libmesh has 42 QPs.
  // When plotted, Dunavant's 42 QP rule *appears* to have points on the
  // boundary, but they are in fact only very *close* to the boundary
  // and not on it, for example one QP is located at position
  // (1.196e-1, 1.0984e-3) in the equilateral reference triangle.
  // -d14 -c1 -v0 -e0 -m0 -g13 # 40 QP, <-- No solutions found
  // -d14 -c0 -v0 -e2 -m0 -g12 # 42 QP, <-- No solutions found
  // -d14 -c0 -v0 -e1 -m1 -g12 # 42 QP, <-- No solutions found
  // -d14 -c0 -v0 -e0 -m2 -g12 # 42 QP, <-- No solutions found
  // -d14 -c1 -v1 -e1 -m0 -g12 # 43 QP, <-- No solutions found
  // -d14 -c1 -v0 -e3 -m0 -g11 # 43 QP, <-- One new solution found
  // -d14 -c0 -v1 -e3 -m0 -g11 # 45 QP, <-- One new solution found
  // -d14 -c0 -v0 -e5 -m0 -g10 # 45 QP, <-- No solutions found
  // -d14 -c1 -v1 -e4 -m0 -g10 # 46 QP, <-- No solutions found
  // ...

  // d==15, dim=46, best PI rule in libmesh has 49 QPs
  // -d15 -c1 -v0 -e0 -m0 -g15 # 46 QP, <-- No solutions found
  // -d15 -c0 -v1 -e0 -m0 -g15 # 48 QP, <-- No solutions found
  // -d15 -c1 -v1 -e1 -m0 -g14 # 49 QP, <-- No solutions found
  // -d15 -c1 -v1 -e0 -m1 -g14 # 49 QP, <-- No solutions found
  // -d15 -c1 -v0 -e3 -m0 -g13 # 49 QP, <-- No solutions found
  // -d15 -c0 -v1 -e3 -m0 -g13 # 51 QP, <-- No solutions found
  // -d15 -c0 -v0 -e5 -m0 -g12 # 51 QP
  // -d15 -c1 -v1 -e4 -m0 -g12 # 52 QP
  // -d15 -c1 -v0 -e6 -m0 -g11 # 52 QP
  // -d15 -c0 -v1 -e6 -m0 -g11 # 54 QP
  // -d15 -c1 -v1 -e7 -m0 -g10 # 55 QP

  // d==16, dim=51, best PI rule in libmesh has 55 QPs
  // -d16 -r 0,0,0,0,17 # 51 QP, <-- instance-1
  // -d16 -r 1,0,1,0,16 # 52 QP, <-- instance-2
  // -d16 -r 1,0,0,1,16 # 52 QP, <-- instance-3
  // -d16 -r 0,1,1,0,16 # 54 QP, <-- instance-4
  // -d16 -r 0,1,0,1,16 # 54 QP, <-- instance-5
  // -d16 -r 0,0,0,3,15 # 54 QP, <-- instance-6
  // -d16 -r 0,0,1,2,15 # 54 QP, <-- instance-7
  // -d16 -r 0,0,2,1,15 # 54 QP, <-- instance-8
  // -d16 -r 0,0,3,0,15 # 54 QP
  // -d16 -r 1,1,0,2,15 # 55 QP
  // -d16 -r 1,1,1,1,15 # 55 QP
  // -d16 -r 1,1,2,0,15 # 55 QP
  // -d16 -r 1,0,4,0,14 # 55 QP
  // -d16 -r 1,0,3,1,14 # 55 QP
  // -d16 -r 1,0,2,2,14 # 55 QP
  // -d16 -r 1,0,1,3,14 # 55 QP
  // -d16 -r 1,0,0,4,14 # 55 QP

  // d==17, dim=57, best PI rule in libmesh has 63 QPs
  // -d17 -r 0,0,0,0,19 # 57 QP
  // -d17 -r 1,0,1,0,18 # 58 QP
  // -d17 -r 1,0,0,1,18 # 58 QP
  // -d17 -r 0,1,1,0,18 # 60 QP
  // -d17 -r 0,1,0,1,18 # 60 QP

  // d==18, dim=64, there is no PI degree 18 rule in libmesh, next highest has 73 pts.
  // -d18 -r 1,0,0,0,21 # 64 QP
  // -d18 -r 0,1,0,0,21 # 66 QP
  //
  // -d18 -r 0,0,2,0,20 # 66 QP
  // -d18 -r 0,0,1,1,20 # 66 QP
  // -d18 -r 0,0,0,2,20 # 66 QP
  // -d18 -r 1,1,1,0,20 # 67 QP
  // -d18 -r 1,1,0,1,20 # 67 QP
  //
  // -d18 -r 1,0,3,0,19 # 67 QP
  // -d18 -r 1,0,2,1,19 # 67 QP
  // -d18 -r 1,0,1,2,19 # 67 QP
  // -d18 -r 1,0,0,3,19 # 67 QP
  // -d18 -r 0,1,3,0,19 # 69 QP
  // -d18 -r 0,1,2,1,19 # 69 QP
  // -d18 -r 0,1,1,2,19 # 69 QP
  // -d18 -r 0,1,0,3,19 # 69 QP
  //
  // -d18 -r 0,0,5,0,18 # 69 QP
  // -d18 -r 0,0,4,1,18 # 69 QP
  // -d18 -r 0,0,3,2,18 # 69 QP
  // -d18 -r 0,0,2,3,18 # 69 QP
  // -d18 -r 0,0,1,4,18 # 69 QP
  // -d18 -r 0,0,0,5,18 # 69 QP
  // -d18 -r 1,1,4,0,18 # 70 QP
  // -d18 -r 1,1,3,1,18 # 70 QP
  // -d18 -r 1,1,2,2,18 # 70 QP
  // -d18 -r 1,1,1,3,18 # 70 QP
  // -d18 -r 1,1,0,4,18 # 70 QP
  //
  // -d18 -r 1,0,6,0,17 # 70 QP
  // -d18 -r 1,0,5,1,17 # 70 QP
  // -d18 -r 1,0,4,2,17 # 70 QP
  // -d18 -r 1,0,3,3,17 # 70 QP
  // -d18 -r 1,0,2,4,17 # 70 QP
  // -d18 -r 1,0,1,5,17 # 70 QP
  // -d18 -r 1,0,0,6,17 # 70 QP
  // -d18 -r 0,1,6,0,17 # 72 QP
  // -d18 -r 0,1,5,1,17 # 72 QP
  // -d18 -r 0,1,4,2,17 # 72 QP
  // -d18 -r 0,1,3,3,17 # 72 QP
  // -d18 -r 0,1,2,4,17 # 72 QP
  // -d18 -r 0,1,1,5,17 # 72 QP
  // -d18 -r 0,1,0,6,17 # 72 QP
  //
  // -d18 -r 0,0,8,0,16 # 72 QP
  // -d18 -r 0,0,7,1,16 # 72 QP
  // -d18 -r 0,0,6,2,16 # 72 QP
  // -d18 -r 0,0,5,3,16 # 72 QP
  // -d18 -r 0,0,4,4,16 # 72 QP
  // -d18 -r 0,0,3,5,16 # 72 QP
  // -d18 -r 0,0,2,6,16 # 72 QP
  // -d18 -r 0,0,1,7,16 # 72 QP
  // -d18 -r 0,0,0,8,16 # 72 QP
  // -d18 -r 1,1,7,0,16 # 73 QP
  // -d18 -r 1,1,6,1,16 # 73 QP
  // -d18 -r 1,1,5,2,16 # 73 QP
  // -d18 -r 1,1,4,3,16 # 73 QP
  // -d18 -r 1,1,3,4,16 # 73 QP
  // -d18 -r 1,1,2,5,16 # 73 QP
  // -d18 -r 1,1,1,6,16 # 73 QP
  // -d18 -r 1,1,0,7,16 # 73 QP
  //
  // -d18 -r 1,0,9,0,15 # 73 QP
  // -d18 -r 0,1,9,0,15 # 75 QP
  // ...

  // Print information.
  // std::cout << "Rule has " << r.n_qp() << " quadrature points." << std::endl;
  //
  // std::cout << "CENTROID: ["
  //           << r.begin(Ro3::CENTROID)
  //           << ","
  //           << r.end(Ro3::CENTROID)
  //           << ")" << std::endl;
  //
  // std::cout << "VERTEX: ["
  //           << r.begin(Ro3::VERTEX)
  //           << ","
  //           << r.end(Ro3::VERTEX)
  //           << ")" << std::endl;
  //
  // std::cout << "EDGE: ["
  //           << r.begin(Ro3::EDGE)
  //           << ","
  //           << r.end(Ro3::EDGE)
  //           << ")" << std::endl;
  //
  // std::cout << "GENERAL: ["
  //           << r.begin(Ro3::GENERAL)
  //           << ","
  //           << r.end(Ro3::GENERAL)
  //           << ")" << std::endl;
  //
  // return 0;

  SolverData solver_data(r);
  solver_data.verbose = true;
  solver_data.maxits = 50;
  solver_data.do_backtracking = true;
  solver_data.residual_reduction_required = true;

  // Tolerances used by the optimization routines (and the 'inner'
  // optimization routine, if using AUGLAG). Recall that these are
  // objective function values, so they are proportional to the "true"
  // residual squared, i.e. an ftol of 1.e-16 corresponds to a true
  // residual of about 1.e-8. If maxeval is negative, then this
  // criterion is disabled.
  //
  // To give some idea about the relationship between maxeval and time
  // spent optimizing, we ran a few seeds of a large (-d18) test case
  // and recorded how long they took. In each case, the 'convergence'
  // reason was NLOPT_MAXEVAL_REACHED. I am not sure whether
  // specifying maxeval on the "local" optimizer effectively squares
  // the number of evaluations that are done, or if the value is
  // somehow cumulative.  Note that there seems to be a lot of
  // variability possible even for a single value of maxeval. In fact,
  // there can be overlap in the times required for 5k and 10k
  // evals. Clearly if you make the maxeval too small you will not
  // obtain any starting points close to anything like the global
  // minimum. I think the main point here is to use a reasonably
  // large value which is not infinite, so you don't spend too much
  // time in local minima where you can't achieve the required ftol.
  //
  // maxeval   Time(s) for "nlopt_test -d18 -c1 -v0 -e0 -g21"
  // 10000  &  83.72 & 36.37 & 45.09 & 81.59 & 38.72 & 142.19
  //  5000  &  19.86 & 47.02 & 19.28 & 41.47 & 18.37
  //  2500  &  9.60  & 10.32 & 9.58  & 9.46  & 21.34
  //  1000  &  3.85  & 9.23  & 9.49  & 9.37  & 6.01
  double ftol_rel = 1.e-30;
  double xtol_rel = 1.e-8;
  int maxeval = 20000;

  // "GN" = global derivative-free optimization
  // nlopt_algorithm alg = NLOPT_GN_DIRECT_L;

  // "LD" = local gradient based optimization

  // .) MMA = method of moving asymptotes
  // o gradient-based
  // o globally-convergent
  // o includes nonlinear constraints but not equality constraints
  // nlopt_algorithm alg = NLOPT_LD_MMA;

  // .) sequential quadratic programming
  // nlopt_algorithm alg = NLOPT_LD_SLSQP;

  // Preconditioned truncated Newton
  // nlopt_algorithm alg = NLOPT_LD_TNEWTON_PRECOND_RESTART;

  // Low storage BFGS
  // nlopt_algorithm alg = NLOPT_LD_LBFGS;

  // You must use a local/subsidiary optimization algorithm with AUGLAG,
  // this is set by calling nlopt_set_local_optimizer().
  nlopt_algorithm alg = NLOPT_AUGLAG;

  // MLSL also requires a local optimizer. MLSL is a multistart
  // algorithm that performs a sequence of local optimizations, and
  // uses a clustering heuristic to avoid repeated searches of the
  // same local optima. Unfortunately, I let it run for a long time
  // (over night) and it never produced any output even on the d==7
  // problem which has a known solution, so I'm not really sure what
  // it's doing or if it's even working?
  // nlopt_algorithm alg = NLOPT_G_MLSL;
  // nlopt_algorithm alg = NLOPT_G_MLSL_LDS;

  // Improved Stochastic Ranking Evolution Strategy.
  // "This method supports arbitrary nonlinear inequality and equality
  // constraints in addition to the bound constraints"
  // nlopt_algorithm alg = NLOPT_GN_ISRES;

  // The problem dimension depends only on "d".
  unsigned int dim = r.dim();

  // Set minimization algorithm and dimensionality
  nlopt_opt opt = nlopt_create(alg, dim);

  if (alg == NLOPT_AUGLAG ||
      alg == NLOPT_G_MLSL)
    {

      ////////////////////////////////////////////////////////////////////////////////
      // Global optimization
      ////////////////////////////////////////////////////////////////////////////////

      // The locally biased variant seemed to work well for d==7 case
      // given a large maxeval.  You have to be careful when using
      // these, at least the first two *ignore* your initial guess, so
      // if you run multiple times with the same maxeval, you will
      // just get the same result every time... Another issue that
      // makes these unsuitable for our use is they will often choose
      // the same values for more than (w,x,y) triple, which generally
      // leads to a singular Jacobian (although for some reason it can
      // produce a small residual).
      //
      // All of the algorithms perform a "deterministic search based
      // on systematic division of the search domain into smaller and
      // smaller hyperrectangles". The "L" versions of the algorithms
      // are specifically stated to "do well for small problems with a
      // single global minimizer and only a few local minimizers"
      // which I don't think describes our problem particularly
      // well. A possible use of the algorithms is to "generate
      // initial iterates for other sampling methods", which is
      // something that seems relevant.
      // nlopt_algorithm local_alg = NLOPT_GN_DIRECT;
      // nlopt_algorithm local_alg = NLOPT_GN_DIRECT_L; // "locally biased" variant
      // nlopt_algorithm local_alg = NLOPT_GN_DIRECT_L_RAND; // randomized variant of DIRECT-L
      // nlopt_algorithm local_alg = NLOPT_GN_DIRECT_NOSCAL; // Unscaled variant
      // nlopt_algorithm local_alg = NLOPT_GN_DIRECT_L_NOSCAL; // Unscaled variant
      // nlopt_algorithm local_alg = NLOPT_GN_DIRECT_L_RAND_NOSCAL; // Unscaled variant
      // nlopt_algorithm local_alg = NLOPT_GN_ORIG_DIRECT; // implementation based on original Fortran
      // nlopt_algorithm local_alg = NLOPT_GN_ORIG_DIRECT_L; // implementation based on original Fortran

      // Controlled random search. This one can also use the
      // nl_set_population() flag. In my experience this one seemed to
      // frequently produce "solutions" that contained many 0s or 1s,
      // i.e. that were pegged at the lower or upper constraint, for
      // some reason.  Possibly it is not designed to be used as an
      // inner unconstrained optimization algorithm wrapped by AUGLAG,
      // or perhaps there is some bug in the implementation. In
      // contrast, the evolutionary algorithm seemed not to do this
      // for whatever reason.
      // nlopt_algorithm local_alg = NLOPT_GN_CRS2_LM;

      ////////////////////////////////////////////////////////////////////////////////
      // Evolutionary algorithms
      ////////////////////////////////////////////////////////////////////////////////

      // The population size for ISRES defaults to 20*(n+1) in n dimensions,
      // this can be changed by calling nlopt_set_population(). Both of these seem to provide
      // somewhat reasonable looking initial guesses (i.e. not full of
      // 0s and 1s like CRS2) and therefore have been my go-to for
      // initial guess selection selection up to this point.
      // For a -d15 problem, ISRES takes about 11.5 seconds while ESCH takes about 14, although
      // some of that time was spent in newton_min.
      nlopt_algorithm local_alg = NLOPT_GN_ISRES;
      // nlopt_algorithm local_alg = NLOPT_GN_ESCH;

      // Global optimization algorithms that use derivative
      // information.  Requires a bound constrained problem. Original
      // source code No longer seems to be available. Requires linking
      // against the NLOPT C++ library, when I try to run it, it just
      // returns error codes, and that could be why...
      // nlopt_algorithm local_alg = NLOPT_GD_STOGO;
      // nlopt_algorithm local_alg = NLOPT_GD_STOGO_RAND;

      ////////////////////////////////////////////////////////////////////////////////
      // Derivative-free algorithms
      // I doubt that any derivative-free optimization algorithm will
      // ever be as good as the derivative-based ones, at least for a
      // given maxeval. But I figured I might as well try them all
      // anyway. Most ran surprisingly (to me) slow compared to the
      // derivative based algorithms.
      ////////////////////////////////////////////////////////////////////////////////

      // nlopt_algorithm local_alg = NLOPT_LN_COBYLA; // 12/33 converged, slow
      // nlopt_algorithm local_alg = NLOPT_LN_BOBYQA; // 1/22 converged, slow
      // nlopt_algorithm local_alg = NLOPT_LN_NEWUOA; // solutions violate constraints?
      // nlopt_algorithm local_alg = NLOPT_LN_NEWUOA_BOUND; // did not work?
      // nlopt_algorithm local_alg = NLOPT_LN_PRAXIS; // 2/40 converged
      // nlopt_algorithm local_alg = NLOPT_LN_NELDERMEAD; // none converged
      // nlopt_algorithm local_alg = NLOPT_LN_SBPLX; // 3/28 converged

      ////////////////////////////////////////////////////////////////////////////////
      // Derivative-based algorithms
      ////////////////////////////////////////////////////////////////////////////////

      // Must set a maxeval for this otherwise it will run forever.
      // Generally does not work as well as SLSQP even for the
      // relatively simple degree=7 problem.
      // nlopt_algorithm local_alg = NLOPT_LD_MMA;

      // Instead of constructing local MMA approximations, constructs
      // simple quadratic approximations. *Must* set a maxeval for this
      // or else it will run forever. This is the only optimization algorithm
      // in nlopt which supports a user-defined Hessian approximation/preconditioner.
      // See: nlopt_set_precond_min_objective(). The function prototype for the
      // preconditioner is
      // void pre(unsigned n, const double *x, const double *v,
      //          double *vpre, void *f_data);
      // and you are meant to compute and return the action of H, vpre = H * v.
      // nlopt_algorithm local_alg = NLOPT_LD_CCSAQ;

      // Low storage BFGS. Runs really fast. Does not seem to respect
      // the ftol we set, as it returns NLOPT_FTOL_REACHED even with
      // the minimum is much larger than the tolerance? Does actually
      // find some roots, but that may be because it runs so fast that
      // I got more data for it... Also when it did converge, it
      // seemed to be my Newton solver that was actually obtaining the
      // root, and not really the optimization program.
      // nlopt_algorithm local_alg = NLOPT_LD_LBFGS; // 10/89 converged

      // Preconditioned truncated Newton. For some reason it seems to
      // find many solutions which are pegged at multiple constraints?
      // I'm not sure what would cause something like this... The
      // versions with no preconditioning and no restart or
      // preconditioning found a few converged solutions in the d=7
      // case.  Seems to run reasonably fast.
      // nlopt_algorithm local_alg = NLOPT_LD_TNEWTON_PRECOND_RESTART;
      // nlopt_algorithm local_alg = NLOPT_LD_TNEWTON_PRECOND; // no restarting
      // nlopt_algorithm local_alg = NLOPT_LD_TNEWTON_RESTART; // no preconditioning
      // nlopt_algorithm local_alg = NLOPT_LD_TNEWTON; // no restart or preconditioning

      // Shifted limited-memory variable metric, rank-2 and rank-1 versions.
      // Seems to run quite fast and does find a few converged solutions (30/263
      // and 29/253, respectively).
      // nlopt_algorithm local_alg = NLOPT_LD_VAR2;
      // nlopt_algorithm local_alg = NLOPT_LD_VAR1;

      // So far this is the best local minimizer I've found to use in
      // conjunction with AUGLAG... For the degree=7 problem, it
      // converged for an impressive 79/115 initial guesses, but
      // perhaps it does not scale as well as some of the other
      // methods for larger problems?
      // nlopt_algorithm local_alg = NLOPT_LD_SLSQP;

      // The objective function, bounds, and nonlinear-constraint
      // parameters of local_opt are ignored. A copy is made when
      // calling set_local_optimizer, so it is safe to let the local
      // copy be destroyed.
      nlopt_opt local_opt = nlopt_create(local_alg, dim);
      // if (local_alg == NLOPT_GN_ISRES || local_alg == NLOPT_GN_CRS2_LM)
      //   nlopt_set_population(local_opt, 40*(r.dim() + 1));
      nlopt_set_xtol_rel(local_opt, xtol_rel);
      nlopt_set_ftol_rel(local_opt, ftol_rel);
      nlopt_set_maxeval(local_opt, maxeval);
      nlopt_set_local_optimizer(opt, local_opt);
    }

  // lower bounds: all parameters must be positive, so we just use the
  // default-constructed values of the std vector for lb. All weights
  // for 3-point orbits must be less than 1/6, which is (1/3)*(ref. element area).
  // Everything else must be less than 1.
  std::vector<double> lb;
  std::vector<double> ub;
  r.bounds(lb, ub);

  // Debugging: Print entries of ub
  // std::cout << "ub=" << std::endl;
  // for (const auto & val : ub)
  //   std::cout << val << std::endl;

  nlopt_set_lower_bounds(opt, lb.data());
  nlopt_set_upper_bounds(opt, ub.data());
  nlopt_set_min_objective(opt, myfunc, static_cast<void *>(&solver_data));

  // Add inequality constraint. We need call the "add" function once
  // per inequality constraint, and each time pass it the index of the
  // x DOF of the orbit which is to be constrained.
  std::vector<unsigned int> data;
  r.inequality_constraint_indices(data);

  // Debugging:
  // std::cout << "data=" << std::endl;
  // print(data);

  for (unsigned int c=0; c<data.size(); ++c)
    nlopt_add_inequality_constraint(opt, myconstraint, &data[c], 1e-8);

  // Set tolerance(s)
  nlopt_set_xtol_rel(opt, xtol_rel);
  nlopt_set_ftol_rel(opt, ftol_rel);

  // You can use this parameter to control how long some of the
  // gradient free optimization algorithms run, otherwise it seems
  // that they can just run forever.
  nlopt_set_maxeval(opt, maxeval);

  // Solution vector
  std::vector<double> x;

  while (true)
    {
      // Generate random initial guess.
      r.guess(x);

      // // Debugging: Print entries of x
      // std::cout << "initial guess=" << std::endl;
      // for (const auto & val : x)
      //   std::cout << val << std::endl;

      // A degree=3 rule with 7 QPs.
      // One can solve for this rule "by hand" since if you
      // assume the edge orbit is at 1/2 it becomes a linear system of
      // equations that has a solution. This rule is not optimal since
      // it has both points on the boundary and it has 7 QPs while the
      // best known degree 3 rule has only 4 QPs. However, the weights
      // are all positive and it was a useful test case for making
      // sure that the VERTEX residual and Jacobian contributions were
      // correct. The exact solution = [9/40, 1/40, 1/15, 1/2] is in Q.
      // if (r.has_orbits(1,1,1,0,0))
      //   x =
      //     {
      //       2.2500000000000000000000000000000e-1,
      //       2.5000000000000000000000000000000e-2,
      //       6.6666666666666666666666666666667e-2,
      //       5.0000000000000000000000000000000e-1
      //     };

      // A new (?) degree=3 rule with 6 QPs composed of two median orbits.
      // I'm not totally sure what is happening here... maybe there are
      // infinitely many solutions and we can show this algebraically?
      if (r.has_orbits(0,0,0,2,0))
        {
          // 1.) OK
          //x =
          //  {
          //     7.7933915797683249045490620091667e-2,
          //     1.2122610860169603123104488174293e-1,
          //     8.8732750868983417621176046575000e-2,
          //     4.4585335732240594543339113007394e-1
          //  };

          // 2.) This initial guess converges to itself.
          // x =
          //   {
          //     5.5430438592710409604505936979277e-2,
          //     9.2321219251293589274056663622764e-2,
          //     1.1123622807395625706216072968739e-1,
          //     4.4591389966593863706329438603431e-1
          //   };

          // 3.) This initial guess converges to itself.
          // x =
          //   {
          //     5.7117599987543327605737790321850e-2,
          //     4.5233569099773816153684027274282e-1,
          //     1.0954906667912333906092887634482e-1,
          //     1.4657918044063831712367103699393e-1
          //   };
        }

      // degree=3 rules with 7 QPs. Many solutions found for this
      // case, even though only a *single* solution found for
      // (1,1,1,0,0) case.
      if (r.has_orbits(1,1,0,1,0))
        {
          // 1.)
          // x =
          //   {
          //     1.9842138198894232087756355234257e-1,
          //     2.4695416925122307603625640854424e-2,
          //     7.5830789078563585437186508364720e-2,
          //     4.9102649835703920944141032123298e-1
          //   };

          // 2.)
          // x =
          //   {
          //     2.1821738936979620876016536190816e-1,
          //     2.4917905757005788534116782239244e-2,
          //     6.9009631119728808545828097124702e-2,
          //     4.9754924428627855803597412887029e-1
          //   };

          // 3.)
          // x =
          //   {
          //     5.4029040412956874973240568941493e-2,
          //     2.3556967318214072640266582461888e-2,
          //     1.2510001921080030236865322789095e-1,
          //     4.6015856878625571724228393577505e-1
          //   };

          // 4.)
          // x =
          //   {
          //     2.2461639812973749630910230345743e-1,
          //     2.4995269549536343748391966980495e-2,
          //     6.6799264407217824148573931867029e-2,
          //     4.9985812675355217882255942640768e-1
          //   };

          // 5.)
          // x =
          //   {
          //     2.1397421359078118356241649552057e-1,
          //     2.4868152433225093792956322275115e-2,
          //     7.0473776369847845019571512551361e-2,
          //     4.9607561820817369646618999468046e-1
          //   };

          // And many more...
        }

      // degree=4 with 6 QPs
      // This rule is different from the D3-invariant rule with 6 QP
      // in libmesh which is from Lyness and Jesperson and has all
      // points inside the reference element.  The rule also has 6
      // points (ne=1, ng=1) but it contains 3 points on the element
      // boundary, so the D3 rule should probably be considered superior,
      // but this one seems .
      // if (r.has_orbits(0,0,1,0,1))
      //   x =
      //     {
      //       3.4009490250701580549995389768470e-2,
      //       8.5705066495465124974799173176085e-1,
      //       1.3265717641596508611667127689820e-1,
      //       5.6577725671689651088265387031247e-1,
      //       3.1774517297791196146365253427748e-1
      //     };

      // degree=5 with 7 QPs
      // if (r.has_orbits(1,0,0,0,2))
      //   x =
      //     {
      //       1.1250000000000000000000000000000e-1,
      //       6.6197076394253090368824693916576e-2,
      //       4.7014206410511508977044120951345e-1,
      //       5.9715871789769820459117580973105e-2,
      //       6.2969590272413576297841972750091e-2,
      //       1.0128650732345633880098736191512e-1,
      //       7.9742698535308732239802527616975e-1
      //     };

      // degree=5 rule with 10 QPs
      // if (r.has_orbits(1,0,0,3,0))
      //   x =
      //     {
      //       1.2186423932585219984549820562014e-1,
      //       1.1700433911324622627645941889857e-2,
      //       3.8273347483138637226611226870869e-2,
      //       4.2316140269087041665556366610484e-2,
      //       4.9096810155229997362992320344920e-1,
      //       7.2028679377637602424964956292946e-2,
      //       1.4482809461230114506784241834028e-1
      //     };

      // degree=6 with 12 points, which matches the degree=6 rule
      // which is currently in libmesh. There are some points on the
      // boundary of the domain, so for that reason, the rule
      // currently in libmesh may be preferred.
      // if (r.has_orbits(0,0,2,0,2))
      //   x =
      //     {
      //       1.7521923484831786607547758325659e-2,
      //       4.2749900628191084402313725002908e-1,
      //       1.2850039868925906983627292049321e-2,
      //       9.3271642599256179913317342394557e-2,
      //       5.2253795110990941510714178014166e-2,
      //       6.8836320947862802695823181103657e-2,
      //       1.9226287054542767293536118369964e-1,
      //       8.4040908201918031564777438277521e-2,
      //       5.0269934990534993982516857981038e-1,
      //       3.3904897391721909524901496402551e-1
      //     };

      // Gatermann's degree=7 rule with 12 points.
      // if (r.has_orbits(0,0,0,0,4))
      //   x =
      //     {
      //       2.6517028157436251428754180460739e-2,
      //       6.2382265094402118173683000996350e-2,
      //       6.7517867073916085442557131050869e-2,
      //       4.3881408714446055036769903139288e-2,
      //       5.5225456656926611737479190275645e-2,
      //       3.2150249385198182266630784919920e-1,
      //       2.8775042784981585738445496900219e-2,
      //       3.4324302945097146469630642483938e-2,
      //       6.6094919618673565761198031019780e-1,
      //       6.7493187009802774462697086166421e-2,
      //       5.1584233435359177925746338682643e-1,
      //       2.7771616697639178256958187139372e-1
      //     };

      // A degree=7 rule with 13 QPs and some points on the boundary.
      // if (r.has_orbits(1,0,1,0,3))
      //   x =
      //     {
      //       5.5952307797500451168270351180298e-2,
      //       1.3472905782326158579409369641855e-2,
      //       2.9983830132009239524945782161234e-1,
      //       2.6082365542691963824789940662658e-2,
      //       5.9482416053645952551714087077197e-2,
      //       6.9397373107702483247241789818035e-2,
      //       6.1273891079961226273752784819268e-2,
      //       2.7907750657353493179879419573566e-1,
      //       1.4598341071586341030182362058933e-1,
      //       4.7186734995853834265957787816119e-2,
      //       3.2916411424668018651106653409879e-1,
      //       6.1215138856035843807465228900117e-1
      //     };

      // A degree=7 rule with 15 QPs and some points on the boundary.
      // if (r.has_orbits(0,0,3,0,2))
      //   x =
      //     {
      //       5.3055935561341320656442591841393e-3,
      //       3.4023228823368393399065798430885e-2,
      //       1.5002198243061481721648841389365e-2,
      //       7.7426490649533686580770326141160e-1,
      //       1.5301251644897155959135909968714e-2,
      //       4.0051406706170526286311293473985e-1,
      //       5.0624560660856324047249434886360e-2,
      //       1.5395946225900545445477142492632e-1,
      //       1.0329490144591944867168415690346e-1,
      //       8.0433062561717572872988221238088e-2,
      //       3.8824142455567101056261725619719e-1,
      //       4.6001634230332020766420547596794e-1
      //     };

      // A degree=7 rule with 16 QPs and points on boundary/at vertices.
      // if (r.has_orbits(1,1,2,0,2))
      //   x =
      //     {
      //       8.8156307800262736848734269834507e-2,
      //       2.0181692212805178377771076491231e-3,
      //       1.1842518154963341192378271866334e-2,
      //       8.3929917227292351714406960839597e-1,
      //       1.1842518154963341192378271866334e-2,
      //       1.6070082772707648285593039160403e-1,
      //       5.3251995907086993738381694116003e-2,
      //       6.1698507712375957239912339414237e-2,
      //       4.6915074614381202138004383029288e-1,
      //       5.8326029294951560422839897890704e-2,
      //       6.9012787955247933305015129005951e-1,
      //       1.5493606022376033347492435497024e-1
      //     };

      // A degree=8 rule with 16 QPs
      // if (r.has_orbits(1,0,1,0,4))
      //   x =
      //     {
      //       7.2521345610590588615251673993725e-2,
      //       1.1691671982623978980576722781465e-2,
      //       2.8464284047341892779748152552693e-1,
      //       1.6036995852050015675610346503080e-2,
      //       5.3867543687557559308134383474885e-2,
      //       4.6837555921867442934220084396571e-2,
      //       4.7920586242366783286185404347777e-2,
      //       4.3827550158929177203507813868026e-1,
      //       4.8025920827024682057103085428712e-1,
      //       1.4800720141846362085949314479152e-2,
      //       7.4347327164768248134793287452843e-1,
      //       1.4096246060080109128363068092795e-2,
      //       5.2042910577582663766594320557284e-2,
      //       1.8711966318476081140504944380348e-1,
      //       1.5372274637782641123563856802741e-1
      //     };

      // A degree=8 rule with 18 QPs
      // if (r.has_orbits(0,1,1,0,4))
      //   x =
      //     {
      //       1.6365619450731091079309677113646e-3,
      //       7.1563040840844654360044305010581e-3,
      //       1.4150003160281202221467222411207e-1,
      //       2.3032399676650253713746312132884e-2,
      //       2.6135943412979344118211349907110e-2,
      //       3.9683870899759689901048321440939e-1,
      //       3.0412879070999470272247013794143e-2,
      //       5.8705218537414718746555765593064e-2,
      //       1.2867350137855184316304358427429e-1,
      //       6.3077071453596363648364086313746e-2,
      //       5.0034691425209789094544109307179e-1,
      //       1.9824402124136825067504259045131e-1,
      //       4.1351450436263004488373856213471e-2,
      //       6.3013360188613761901960007277632e-1,
      //       2.9639453351548959482253455585386e-1
      //     };

      // A degree=8 rule with 19 QPs
      // if (r.has_orbits(1,1,2,0,3))
      //   x =
      //     {
      //       5.3249913011017943647121653650309e-2,
      //       1.5093537559558725922092008419750e-3,
      //       7.7922345016508128137115446990742e-3,
      //       8.5151321170633136497105403382591e-1,
      //       1.2773850763700670089813161350187e-2,
      //       4.1035393094358068585396592117871e-1,
      //       5.2019626337460031874380610292933e-2,
      //       3.0858419507892152432388479445908e-1,
      //       1.4167589821801641930724853282077e-1,
      //       2.8196727788590522799411409283777e-2,
      //       1.2245778857082282250180202616490e-1,
      //       5.6094505881456623061225804438646e-2,
      //       4.6624902515636108614766855648618e-2,
      //       6.1512488037728468907811464063495e-1,
      //       8.1799808561019703180235467829817e-2
      //     };

      // A degree=9 1-0-0-6 rule with 19 QPs, this ties the current minimum!
      // Actually this is the same as the previously-known D3-invariant rule
      // with 19 QPs.
      // if (r.has_orbits(1,0,0,0,6))
      // x =
      //   {
      //     4.8567898141399416909620991253644e-2,
      //     3.9823869463605126516445887132023e-2,
      //     1.8820353561903273024096128046734e-1,
      //     1.8820353561903273024096128046734e-1,
      //     1.5667350113569535268427415643605e-2,
      //     4.8968251919873762778370692483619e-1,
      //     2.0634961602524744432586150327614e-2,
      //     2.1641769688644688644688644688645e-2,
      //     7.4119859878449802069007987352342e-1,
      //     2.2196298916076569567510252769319e-1,
      //     1.2788837829349015630839399279500e-2,
      //     4.4729513394452709865106589966276e-2,
      //     4.4729513394452709865106589966276e-2,
      //     2.1641769688644688644688644688645e-2,
      //     7.4119859878449802069007987352342e-1,
      //     3.6838412054736283634817598783385e-2,
      //     3.8913770502387139658369678149702e-2,
      //     4.3708959149293663726993036443535e-1,
      //     1.2582081701412672546013927112929e-1
      //   };

      // A degree=9 rule with 21 QPs
      // if (r.has_orbits(0,1,0,0,6))
      // x =
      //   {
      //     4.6208136672890810190208010916353e-4,
      //     1.3937873649405077572538903674386e-2,
      //     4.4412127969532946609530810618352e-2,
      //     6.0455484530038181374983232537504e-2,
      //     8.7959612804017479710722128519961e-3,
      //     1.9013225698280544912811185627099e-1,
      //     7.1987546384528925078927829097248e-3,
      //     4.0344560584550655466935190288604e-2,
      //     2.1110808416324593232253596761464e-1,
      //     1.2631554612265104469150870016886e-1,
      //     2.5068384777481441427534757854270e-2,
      //     4.4513058891472071846505101597941e-1,
      //     3.5535480200215007101851885821557e-2,
      //     2.5175973969753207576301584782348e-2,
      //     4.0997052576549121963973333096097e-2,
      //     2.4960151776578923653305787716817e-1,
      //     5.2881831038345628550381937105899e-2,
      //     1.9128515995430339172415385144765e-1,
      //     3.5897434142891728391398528200056e-1
      //   };

      // A degree=9 rule with 22 QPs
      // if (r.has_orbits(1,1,1,0,5))
      //   x =
      //     {
      //       5.6681242229959782378455068998313e-2,
      //       1.7709563164108618654073199527745e-4,
      //       8.0056858564321250888415067519545e-3,
      //       5.0000000000000000000000000000000e-1,
      //       2.1641769688644688644688644688645e-2,
      //       3.6838412054736283634817598783385e-2,
      //       2.2196298916076569567510252769319e-1,
      //       1.3032090619145091983675537032338e-2,
      //       4.6947443199090612267241036420043e-2,
      //       4.6947443199090612267241036420043e-2,
      //       4.1203790999884165472995247646690e-2,
      //       1.9187191273744911850366265609742e-1,
      //       1.9187191273744911850366265609742e-1,
      //       2.1641769688644688644688644688645e-2,
      //       2.2196298916076569567510252769319e-1,
      //       3.6838412054736283634817598783385e-2,
      //       4.2070716772288226519084664197013e-2,
      //       1.0044132362596745752372601888727e-1,
      //       4.4977933818701627123813699055636e-1
      //     };

      // A degree=9 rule with 24 QPs
      // if (r.has_orbits(0,1,3,0,4))
      //   x =
      //     {
      //       9.9920681638586186356375443283050e-4,
      //       7.0690247201981272962308384289864e-3,
      //       6.3193040073722041295255627701824e-1,
      //       5.9266081839999399947176520476464e-3,
      //       8.8043673507833069274213166314801e-1,
      //       5.3475597462147498356454535340421e-3,
      //       1.0667687351399658302468746396423e-1,
      //       3.7760534137683438155770128636229e-2,
      //       8.9980550840639517183429370342513e-2,
      //       3.2711502822800339810214903476533e-1,
      //       2.9572113473164001765589944002256e-2,
      //       1.0777246753344816261464172868087e-1,
      //       1.0651976200065317286577734289541e-1,
      //       5.2240014546826852011186560781300e-2,
      //       2.1197643637448513468424883083615e-1,
      //       4.8594896708812698156314636765236e-1,
      //       2.7751605042193695743962334803376e-2,
      //       6.2512240599886847137981590902263e-1,
      //       3.3384503238421439815486250637193e-1
      //     };

      // A degree=10 rule with 24 QPs. New minimum!
      // if (r.has_orbits(0,0,2,0,6))
      //   x =
      //     {
      //       5.3681429742859643492357377689009e-3,
      //       1.8721976304020603799148510523827e-1,
      //       2.2815303655641663532038669323421e-3,
      //       3.4075334397348650193113634396732e-2,
      //       4.6508993462059758854967647559714e-2,
      //       3.4147086142038659869812874556596e-1,
      //       2.0360915320008500090421409688573e-1,
      //       1.0919981154490145628957387886100e-2,
      //       6.9817400501731927826050298390783e-1,
      //       1.4413404224312299237271155349378e-2,
      //       2.4240746515267806539218234118407e-2,
      //       4.2633514579462347877517092365945e-1,
      //       3.6779605207993416066366804528172e-2,
      //       1.4936119213925489681656723297916e-2,
      //       4.0093527438242910369389026226690e-2,
      //       9.1103646481204342896318426601411e-2,
      //       3.4963363159025821525120057092851e-2,
      //       1.1690104736983285255079342004697e-1,
      //       2.7867154732579039483251817674802e-1,
      //       2.7447789822047513734307012010436e-2,
      //       7.2000132708127838543972050558105e-1,
      //       1.9134246734700412201903713737393e-1
      //     };

      // A *second* degree=10 rule with 24 QPs and 0026 configuration.
      // Also a new minimum!
      // if (r.has_orbits(0,0,2,0,6))
      //   x =
      //     {
      //       6.8034549450389437911542695284968e-3,
      //       2.4033695082490661641610501527861e-1,
      //       6.4701145825025661268470771647561e-3,
      //       5.2540484162127156065220999286902e-1,
      //       7.6297245866351294827778042218373e-3,
      //       1.2970322053864413659668358505098e-2,
      //       1.4081107816864061630771237989169e-1,
      //       2.7968592259403195267388456279467e-2,
      //       7.5136392619669750371252757075214e-1,
      //       1.5004696166167503172923442728480e-1,
      //       3.6263197047755420520828693175881e-2,
      //       5.3735517659978976913561285260804e-1,
      //       3.7661987230500154446477586010322e-1,
      //       7.4674601101896288064569309205665e-3,
      //       4.4107239622531995584271061559211e-2,
      //       2.6087950887797941493159400743790e-2,
      //       2.6160683863563246187991128905119e-2,
      //       6.5311292067804412162841018713658e-1,
      //       5.9875144970959509735752771006843e-2,
      //       4.7903439271578536483222306470543e-2,
      //       2.8533958139225425086186103213101e-1,
      //       4.8509072998273930259505659173245e-1
      //     };

      // A *third* degree=10 rule with 24 QPs and 0026 configuration.
      // I wasn't even looking for this... it just popped up while I
      // was experimenting with different optimization algorithms! It
      // also almost didn't converge... it required exactly 50 Newton
      // iterations!
      // if (r.has_orbits(0,0,2,0,6))
      //   x =
      //     {
      //       5.4182602431180742419012975840987e-3,
      //       1.8436361439082912239287076199852e-1,
      //       5.9852576657074839510432326113733e-3,
      //       6.1050316703542844868812608863530e-1,
      //       8.0691227940699078197481579119651e-3,
      //       9.2893410705677578721193588750841e-1,
      //       3.5753059859079028015587068659359e-2,
      //       2.7929259575676991206594548597666e-2,
      //       1.0070799873421512282498660467458e-1,
      //       7.3715384361281192006447189345953e-1,
      //       1.3329322854499305969742173140315e-2,
      //       1.7425713782396221670836596810662e-1,
      //       7.9674412091269421568078805095696e-1,
      //       3.5352219432776408242842120396444e-2,
      //       3.2019851929290632085891892083740e-1,
      //       5.7670773292297347321026202955429e-1,
      //       2.4337851756964119586513300176318e-2,
      //       4.3963494980239674841811610522069e-2,
      //       5.7511639509054177348822161582207e-1,
      //       4.6245372343854375648281836248486e-2,
      //       3.1146210305934346282199719556151e-1,
      //       2.1624307751824776491821571751996e-1
      //     };

      // Wow, a *fourth* degree=10 rule with 24 QPs
      // if (r.has_orbits(0,0,2,0,6))
      //   x =
      //     {
      //       6.1935726910001758541078241539263e-3,
      //       3.3937615435249893352746365930789e-1,
      //       5.2925293385584474156952124159760e-3,
      //       8.1604706571429431480245669966563e-1,
      //       2.7830353344135349778500405243835e-2,
      //       9.4133792711222378238760106555204e-2,
      //       1.7631879015417056904139925982017e-1,
      //       4.6873870579045364081863575403270e-2,
      //       4.6859333350512851674998237955624e-1,
      //       2.1171166682586492544087434741492e-1,
      //       2.3970523489400444168780206716243e-2,
      //       3.9192270748963862949012192065698e-2,
      //       4.0281896555207495154945708709387e-1,
      //       1.3626423836434353260840172788744e-2,
      //       1.4757623947530962134060686657472e-1,
      //       3.5613227466971039550554430193004e-2,
      //       3.5870519028083430398177360897913e-2,
      //       1.0439471886636364342444603201521e-1,
      //       5.8802274042993997988484650010385e-1,
      //       7.0088743600091017087019090467601e-3,
      //       3.5987336367983871336129310176438e-2,
      //       9.3343603776198938768545803095892e-1
      //     };

      // Ridiculous, a *fifth* degree=10 rule with 24 QPs
      // if (r.has_orbits(0,0,2,0,6))
      //   x =
      //     {
      //       7.1633170697004196998871675324386e-3,
      //       7.5955588633458936806336676273387e-1,
      //       5.1049969920165429199237231873571e-3,
      //       1.5193663248046395463352611125183e-1,
      //       7.8568994386027624616574131921037e-3,
      //       4.2237680049868393035216758301495e-2,
      //       9.2862155781321497106698489852560e-1,
      //       2.4055550786638763160189013574995e-2,
      //       6.7373748617165142309703171777769e-1,
      //       2.6370177258295204404766466671693e-1,
      //       2.4598961140793922509347817346151e-2,
      //       1.3421000044845207799557776284413e-1,
      //       7.6764505900725274771488592099709e-1,
      //       1.5072036191824669938696300381058e-2,
      //       2.2585748215909707408243899147017e-2,
      //       5.1250467241310646796363694038491e-1,
      //       3.4678422346643332938784934500371e-2,
      //       5.6779938220150505301675570462183e-1,
      //       9.9751753925878964685649398982363e-2,
      //       4.8136482700446253038180296952192e-2,
      //       3.0785415280170697910264379441776e-1,
      //       2.1458187835052252754456597776010e-1
      //     };

      // A *sixth* degree=10 rule with 24 QPs. It just keeps going...
      // if (r.has_orbits(0,0,2,0,6))
      //   x =
      //     {
      //       5.3834706252454055829176909195913e-3,
      //       3.0855885949458465505208876769013e-1,
      //       2.2585614023837660045101392700689e-3,
      //       9.7067360043836614156360874897275e-1,
      //       2.3617073911209024419373702587775e-2,
      //       5.5596379430986266793841350014444e-1,
      //       3.6346785109772127345585658484591e-2,
      //       4.4759864182831451230220789451556e-2,
      //       4.0577405760968800735440681116372e-1,
      //       1.9765132573961310778448643201566e-1,
      //       1.4952769339638923689958680821648e-2,
      //       2.9419545386540148820586514239455e-2,
      //       1.7528426317643462839331536695820e-1,
      //       2.6220682345205438475350628365720e-2,
      //       2.9793346146766825012738112385943e-1,
      //       7.4101005002342904421619518548260e-2,
      //       1.5811554292507884609065109648664e-2,
      //       8.6050319153478836901733989334052e-1,
      //       1.0014651384993564882983705026052e-1,
      //       3.3662690567644772655269925601644e-2,
      //       6.4214981788197325641907116304210e-1,
      //       1.5645681947791997908258482531983e-1
      //     };

      // A degree=10 rule with 25 QPs
      // if (r.has_orbits(1,1,1,0,6))
      //   x =
      //     {
      //       4.7309440928094852885976950813341e-2,
      //       7.4443068205713019768945426025270e-4,
      //       3.1829359254805817751202442366943e-3,
      //       8.3061025308964767122329536525606e-2,
      //       3.8799021945743061001599802530126e-2,
      //       1.7544304433743764192225628495153e-1,
      //       2.0933705012838168942077572054598e-1,
      //       2.0030412127900619224360889008975e-2,
      //       3.2176407820770162393255046010956e-2,
      //       3.0726941261349314805052750203576e-1,
      //       1.5850177031042592646271274554700e-2,
      //       3.9026767455703841640294956578970e-2,
      //       9.9070711082074272239198884743969e-2,
      //       1.3553817908665389393448328866720e-2,
      //       1.9269035323026781372640492307404e-2,
      //       5.7154160680610615243446739606247e-1,
      //       2.0804506634455471119732622715725e-2,
      //       1.9770352662977282642262230494035e-1,
      //       4.9561536726352818108215112762139e-2,
      //       3.7931550768623537013118400222361e-2,
      //       4.1147963472022346886180725166791e-1,
      //       1.2777722809830669309582013618709e-1
      //     };

      // A degree=10 rule with 25 QPs
      // if (r.has_orbits(1,0,3,0,5))
      //   x =
      //     {
      //       4.0232634070269439301667832239894e-2,
      //       3.3422311840204575942302033783202e-3,
      //       5.2253403465953518825868105768782e-2,
      //       6.1424211433885610965819016355015e-3,
      //       2.7344809236135841598530720476098e-1,
      //       5.9764454459057452955269360703495e-3,
      //       5.2410938325006369651660186410547e-1,
      //       1.0057301050542603016660734438860e-2,
      //       9.0895615934244191891990428663041e-2,
      //       8.8452523398346606829848673809009e-1,
      //       3.8799010163202809566473695352386e-2,
      //       5.1227366986894646065998675430995e-1,
      //       3.9602679189978358269945323907258e-1,
      //       2.6174937160968349701196686726870e-2,
      //       7.6203914975198293014324875248846e-1,
      //       1.5901928319293625479920526665605e-1,
      //       3.9575280522419545364193216395460e-2,
      //       1.9753738413334949871287005642309e-1,
      //       2.5612186241563726204876368462414e-1,
      //       2.3188161972795448597914015255621e-2,
      //       6.8579014955455508229105658625928e-1,
      //       4.5690890104473341999976206420424e-2
      //     };

      // A *second* degree=10 rule with 25 QPs
      // if (r.has_orbits(1,0,3,0,5))
      //   x =
      //     {
      //       5.1345668536024914411167997684058e-2,
      //       3.4581948178630205438576188361475e-3,
      //       5.9390050736371604294803794039014e-2,
      //       7.0888025444097578448485241065769e-3,
      //       4.1492064251981583325843804152774e-1,
      //       2.1064303554996973040694501159955e-3,
      //       9.4092794891315452906467386615319e-1,
      //       3.9017241786683034241205005116621e-2,
      //       1.0284393318430089371584671935755e-1,
      //       4.8530147947070841080705363096520e-1,
      //       1.8985766577595162514393470423551e-2,
      //       3.0673465604941172960039426783390e-2,
      //       3.0799818976234575755778917061732e-1,
      //       2.2101552554910438933768447236166e-2,
      //       4.8385919093417191677167303459693e-2,
      //       7.5251385329642788292413917481847e-1,
      //       1.6405929030770321222066991770666e-2,
      //       5.2236414355455385930427554647279e-2,
      //       1.1351835493802583579205862674649e-1,
      //       4.0387526153593595925401159832924e-2,
      //       5.9280494105185473668624375618082e-1,
      //       1.8256468126521531326784284288188e-1
      //     };

      // A degree=10 rule with 27 QPs
      // if (r.has_orbits(0,1,3,0,5))
      //   x =
      //     {
      //       7.0349169304985602109942089161121e-4,
      //       4.2428107036736195296584198917347e-3,
      //       8.9246341808720256261073890380435e-1,
      //       7.0016114161227155365721513157167e-3,
      //       3.0494756762697497590213637925978e-1,
      //       6.5757487042465969693037379103481e-3,
      //       6.0349372376167113477272117355728e-1,
      //       4.5367047076144017880509249966141e-2,
      //       2.7702918680309928564466502274207e-1,
      //       4.8213298384380707906450631732514e-1,
      //       3.3434147221659577726825933401689e-2,
      //       4.5537573281654840359879597498098e-1,
      //       8.1628735765941535287742569977895e-2,
      //       2.7372947537570407982425816418687e-2,
      //       6.8321697249249291081054474783040e-2,
      //       2.1461041150208394379962799554115e-1,
      //       2.8008115730342111834738616377667e-2,
      //       9.5199220006813800448433295908254e-2,
      //       6.7380547129004477728227558622451e-1,
      //       1.3960746583857763185533320493072e-2,
      //       8.6097589069840721814036297772080e-2,
      //       3.8205053391734468157969505461413e-2
      //     };

      // A degree=11 rule with 27 QPs. This is a new minimum!
      // if (r.has_orbits(0,0,1,0,8))
      //   x =
      //     {
      //       5.2501900628803483265691402520251e-3,
      //       2.7263237565209247203086721087576e-1,
      //       9.6078311676320423628745787031043e-3,
      //       1.2566813660658779461518373277871e-2,
      //       3.6422495909528624483300674843998e-1,
      //       2.8397128511550294298445438150029e-3,
      //       9.4380452905028433573549593440564e-3,
      //       3.3307532214123278261894624459744e-2,
      //       9.1226212614723748339671998630591e-3,
      //       8.8159235080580116486130047655196e-1,
      //       9.0597208132597846103287274377368e-2,
      //       4.1325979108748536413659626103765e-2,
      //       2.1435335265095806153369337057104e-1,
      //       4.5275183764840794396905381350312e-1,
      //       1.5584140407676158910413646153654e-2,
      //       3.6673070947612617954340116533663e-2,
      //       1.5740290188190833241335042691653e-1,
      //       3.3066696149169029769073153196462e-2,
      //       5.8600348844560926069470015142148e-1,
      //       1.2352912853885500748429654291735e-1,
      //       2.3284027204112150759145400913219e-2,
      //       4.1377415154448198901737239785191e-1,
      //       5.1624512474990088671128327739257e-2,
      //       2.6585468453820995861119377666375e-2,
      //       9.8616918269915361670885340191111e-2,
      //       7.0688498243425849919240661534406e-1
      //     };

      // A second degree=11 rule with 27 QPs. This is a new minimum!
      // if (r.has_orbits(0,0,1,0,8))
      //   x =
      //     {
      //       4.1612804317314842767555766406199e-3,
      //       6.8002906083945000961186411060696e-1,
      //       3.1338137025347552815928236201535e-2,
      //       5.6433855195090291230615957448581e-1,
      //       3.1435732845762994573128512907851e-1,
      //       2.4652523258928223931542332898165e-2,
      //       1.2322869010646763885831827091858e-1,
      //       1.6927635734210890199455789431847e-1,
      //       1.0184449346229461815371958266942e-2,
      //       2.6486113476716200406483433449668e-2,
      //       1.4387747522813605056916070495475e-1,
      //       1.5012583397281533333668358360737e-2,
      //       1.6801782841091323672916814789423e-1,
      //       3.2980816870835011314453213507343e-2,
      //       3.9258394807048784825903710855647e-2,
      //       2.2209434646096813071330105288043e-1,
      //       3.2146756719314824820725269959936e-1,
      //       6.2502000020799461189244380099301e-3,
      //       3.3754120263242869089885850618322e-2,
      //       2.8770220336106529610117019606918e-2,
      //       2.3964763845196892791134629522367e-2,
      //       5.6757782226097864625226237614412e-1,
      //       6.1105709652421022972202840731351e-2,
      //       1.1844334552822786757437425910724e-2,
      //       3.8705519294784296506910844535742e-1,
      //       1.7297366791290889546544445963562e-2
      //     };

      // A degree=11 rule with 28 QPs.
      // if (r.has_orbits(1,0,2,0,7))
      //   x =
      //     {
      //       4.2122982350847949694727649296954e-2,
      //       3.9464280009823691929272033278939e-3,
      //       4.3851364648565862779521140096313e-1,
      //       3.8898791225765856909372783156346e-3,
      //       1.3995698827372101943554341473740e-1,
      //       3.5262305538711489842410819392275e-2,
      //       5.8566468487181347917878803731368e-1,
      //       1.7036392744105047524076839706264e-1,
      //       3.3361116413273105358961963180009e-2,
      //       4.8306262254575673934200276414849e-1,
      //       3.7827174154058150282594080624972e-1,
      //       1.9560808916875607711291904779182e-2,
      //       1.3137337444561324011566482502052e-1,
      //       8.9880877534325906643646059438573e-2,
      //       1.7473024328675662660544053438309e-2,
      //       6.6027376485714331282233861970713e-1,
      //       2.9644506798669171221562610931799e-1,
      //       1.9560067376719361931546552789315e-2,
      //       3.8903304298748738240372330102006e-2,
      //       3.7960586267452617988153688836539e-1,
      //       6.4204924585268918107940059108414e-3,
      //       3.7008141341327453433650952017530e-2,
      //       9.3611834286247349005855727560366e-1,
      //       1.3151550393376275902343669100888e-2,
      //       1.8410877771043024982846962887486e-1,
      //       7.8854437699959397338611667104577e-1
      //     };

      // A degree=11 rule with 30 QPs.
      // if (r.has_orbits(0,1,2,0,7))
      //   x =
      //     {
      //       7.5412493611041301804962310010194e-4,
      //       6.3161815270687152350168791503399e-3,
      //       5.6385014612129497529881148904435e-1,
      //       3.7494474079695397188491120716770e-3,
      //       8.3488305870667139228616877262442e-1,
      //       2.4433178455646771114724388698792e-2,
      //       2.5325397027450240080362940174495e-1,
      //       6.8174631574538163787988886551214e-1,
      //       4.2288098861187531062746735284753e-2,
      //       4.8177392799531548713787151950669e-1,
      //       2.5889350380839276410035974192840e-1,
      //       2.5394796909503515368692355944612e-2,
      //       1.8223607005333125199210988159907e-1,
      //       1.1333303081543736946090107077428e-1,
      //       1.4815826371661608754249307042139e-2,
      //       2.9093074037820991137978739688485e-2,
      //       6.5399402380789810595116446692129e-1,
      //       5.1517923138868893670696918220655e-3,
      //       1.3720341444785722286806905133316e-1,
      //       6.5733729231839964979308737167756e-3,
      //       3.1942784309260200650148324153126e-2,
      //       4.5314310583192857405867396170595e-1,
      //       4.4830066793135978961050156779897e-1,
      //       1.1820435574371482377120249399060e-2,
      //       8.8628106858323135974444732732874e-1,
      //       5.0240542970396197597182374396274e-2
      //     };

      // A degree=11 rule with 30 QPs.
      // if (r.has_orbits(0,0,4,0,6))
      //   x =
      //     {
      //       3.4042935191212521125417637423319e-3,
      //       1.0533831485611468951827070887515e-1,
      //       1.4136460443059679379104136378358e-3,
      //       9.8080568711710162525846816632421e-1,
      //       3.8685011609546712039597447831703e-3,
      //       8.2904397392277897613273802134403e-1,
      //       6.5265066749499600082919813841364e-3,
      //       4.9913242609386226378123612282765e-1,
      //       1.4631682182610466538653073674428e-2,
      //       6.3629896186521054170467908823881e-2,
      //       7.5696723176356166876951114964586e-2,
      //       3.8127981634651491545796781560602e-2,
      //       2.5158471141276727123152166420305e-1,
      //       2.7737293953432875042582608471750e-1,
      //       1.8283036918212341492710216612943e-2,
      //       2.6038385172770044846218234873667e-1,
      //       3.7096444360507273292978510853588e-2,
      //       2.2403768844131921054416252342201e-2,
      //       2.8694297686892572885172482330583e-1,
      //       6.6107239076911926661580007115788e-1,
      //       2.4990073832444596183619315819850e-2,
      //       6.7417466711525312625004511594622e-1,
      //       1.6946898544628711377722057702335e-1,
      //       3.3017175855283998588767123109168e-2,
      //       1.0112420559093728857058675508909e-1,
      //       4.7395330890231320202868021211856e-1
      //     };

      // A degree=11 rule with 31 QPs.
      // if (r.has_orbits(1,1,3,0,6))
      //   x =
      //     {
      //       2.8451667886196583694006718095887e-2,
      //       7.3055714467327739638952453597763e-4,
      //       3.6811918916502771321294475236903e-3,
      //       8.5887028128263670403917393805835e-1,
      //       6.5351044593471576161934645215936e-3,
      //       5.0000000000000000000000000000000e-1,
      //       3.6811918916502771321294475236903e-3,
      //       1.4112971871736329596082606194165e-1,
      //       3.5315570503778766188602955268914e-2,
      //       2.3833330710732888004699981720368e-1,
      //       2.3833330710732888004699981720368e-1,
      //       2.0528157714644283320826157453647e-2,
      //       6.7793765488259040154212614118875e-1,
      //       4.4841677589130443309052391468801e-2,
      //       3.3610672307980639973419663983883e-2,
      //       4.4733898105487205823982810218710e-1,
      //       1.0532203789025588352034379562580e-1,
      //       2.1108427194052090712444172658803e-2,
      //       7.2463396466285646833648738316732e-1,
      //       1.3768301766857176583175630841634e-1,
      //       1.1463746548846752642370103044192e-2,
      //       5.6406468667843899761177546691495e-2,
      //       5.6406468667843899761177546691495e-2,
      //       2.0528157714644283320826157453647e-2,
      //       4.4841677589130443309052391468801e-2,
      //       6.7793765488259040154212614118875e-1
      //     };

      // A _second_ degree=11 rule with 31 QPs.
      // if (r.has_orbits(1,1,3,0,6))
      //   x =
      //     {
      //       3.7542955251948995401672680293581e-2,
      //       3.6533451963881083443684362034211e-4,
      //       3.5762967678650843788591230666722e-3,
      //       9.1015474311166688608512214305766e-1,
      //       2.6073745012037032789340447798328e-3,
      //       6.2385191332680303902529447964571e-2,
      //       4.5551605943322813979348587630311e-3,
      //       6.9954864306968144192703052658983e-1,
      //       3.3338612898292408774050540627908e-2,
      //       1.5645214415692873582456121817917e-1,
      //       3.7135696041962229240381346558013e-1,
      //       3.6051275129777332126733353401397e-2,
      //       5.9709474708140337993704348861480e-1,
      //       2.4881761956721132779037905742563e-1,
      //       2.2879360289804412969026940836061e-2,
      //       6.8342015633946618852882997579606e-2,
      //       2.3706083495843791490270902111653e-1,
      //       1.4275934919031088293713358286828e-2,
      //       2.4122406608702034263703592419432e-1,
      //       2.6461088264622624769359649119523e-2,
      //       1.6271506250344408816976566105223e-2,
      //       8.4217751123168425530061679836483e-2,
      //       7.6290214863539934178354308794439e-2,
      //       2.0231492379060803995443477081512e-2,
      //       5.0229888907582665997864575023573e-1,
      //       4.6319966110752449250886263049113e-1
      //     };

      // A degree=14 rule with 43 QPs
      // if (r.has_orbits(1,0,3,0,11))
      //   x =
      //     {
      //       8.7602731941396900630735811947358e-3,
      //       2.9370048439533343099581854148131e-3,
      //       6.2569253700367656437599058972695e-1,
      //       1.6722513059433337416138535964424e-3,
      //       8.8680974748427433769925101370325e-2,
      //       5.7390658454334448532157896567854e-4,
      //       1.2799961839418472894118805041784e-2,
      //       9.5083874542120944322747876929124e-3,
      //       1.1710792836049738347970910075001e-1,
      //       4.8441813220397451000822101304212e-2,
      //       1.8325481630759103084118515572825e-2,
      //       4.4083626500893020361210042699666e-1,
      //       1.0532653095577596382329864081012e-1,
      //       2.4395501721248301734401840100893e-2,
      //       1.8192544578916314911386413213589e-1,
      //       2.6274708469190746410044887848486e-1,
      //       5.4733878981357532110672298946041e-3,
      //       2.3614034101103065689516149293159e-2,
      //       5.5365008543376860411186591718459e-2,
      //       7.1580877387819418984424546237794e-3,
      //       1.7805982025480138666710873201622e-2,
      //       1.7988024351778182678727751998295e-1,
      //       1.5382037527279956789320393864015e-2,
      //       7.4413372874920121574666100186901e-1,
      //       1.0096598858785286765560093338538e-1,
      //       2.0770322555420189374486180714216e-2,
      //       2.5755784039249370167366496740320e-1,
      //       1.0490301583281400901946442066470e-1,
      //       7.8294830776101924134496693934753e-3,
      //       7.3350913065469368221226676165934e-1,
      //       2.4892308231352947537903472219670e-1,
      //       1.1093259959636639241345639892050e-2,
      //       4.2927867298388176661985950025306e-1,
      //       2.3217124123051175846242571349823e-2,
      //       1.5272066925997130278476353288891e-2,
      //       5.2264287574525983396476919517189e-2,
      //       3.2561032239980910848464202797350e-1,
      //       2.3355396378432121651365456587160e-2,
      //       3.4651207179134267642714477679978e-1,
      //       2.2593051901338893848890937988649e-1
      //     };

      // A degree=14 rule with 45 QPs
      // if (r.has_orbits(0,1,3,0,11))
      //   x =
      //     {
      //       3.3341899592961390239330691261075e-4,
      //       3.0804695636428859775836998484786e-3,
      //       6.3949508844000101530322090060017e-1,
      //       1.3467792723446368940584261873350e-3,
      //       8.9865581126975221365647380324030e-1,
      //       1.5808612604308310282152033095180e-3,
      //       8.2584471866031358855731003566467e-2,
      //       1.5517179962432558030079423848826e-2,
      //       3.2737328694989254537979410999909e-1,
      //       6.1224845612522570738777194715772e-1,
      //       5.8662491339398932703718253350689e-3,
      //       3.4850864435982021385019048236532e-2,
      //       4.5197798839280629001552417799081e-2,
      //       1.3100196828241954768175386975552e-2,
      //       1.0211868224797030928841192894071e-1,
      //       1.1807074456720143580973186137872e-1,
      //       5.7384174774244469979829896128327e-3,
      //       1.2103222706495347272600039299542e-2,
      //       6.9871209216410311139155281266307e-1,
      //       8.8257700451123060137611409110076e-3,
      //       2.6608472404768591546336780338015e-2,
      //       1.9182589356278665756923992509375e-1,
      //       2.0328946435433698309360010474528e-2,
      //       2.7548450177012763711574099989473e-1,
      //       9.8017236746459420569439397915450e-2,
      //       1.2337108101872235990810963734127e-2,
      //       4.5072480205695292700062605730619e-1,
      //       2.7784943748007040744651153745007e-2,
      //       8.7367504162729152038456560338039e-3,
      //       1.5863306115852773909310196554799e-1,
      //       3.3841769725060766849846547666604e-2,
      //       2.6166451247689517218959272227204e-2,
      //       4.4733762244797522632545143088875e-1,
      //       2.8250580286951273987516423182381e-1,
      //       2.0031633502453665041956251163153e-2,
      //       1.6277387226616365612347189484367e-1,
      //       2.2247245619331818629709299504105e-1,
      //       2.3676434423445508019113110092622e-2,
      //       4.3244235515783959316476754659812e-1,
      //       1.3816549030011917242370306253053e-1
      //     };

      // std::cout << "x=" << std::endl;
      // print(x);

      // Will be set to the minimum objective value, upon return.
      double minf;
      nlopt_result res = nlopt_optimize(opt, x.data(), &minf);

      // Call the optimization routine
      if (res < 0) {
        std::cout << "nlopt failed!\n";
      }
      else
        {
          std::cout << "nlopt converged with reason: " << nlopt_result_to_string(res) << std::endl;
          std::cout << "found minimum " << minf << " at " << std::endl;
          print(x);

          // Set solver_data.u = x
          solver_data.u.resize(dim);
          for (unsigned int i=0; i<dim; ++i)
            solver_data.u[i] = x[i];

          // Call our Newton solver, using the nlopt solution as the initial guess.
          // One really frustrating thing that I found in the d==7 case was that,
          // even when started very close to the true solution by using an initial
          // guess converged by nlopt, Newton could still diverge.
          bool converged = false;
          converged = newton_min(solver_data);
          converged = newton(solver_data);
          // converged = gradient_descent(solver_data);
          // converged = nlcg(solver_data);

          if (!converged)
            std::cout << "Solution is *not* converged." << std::endl;
          else
            {
              std::cout << "Solution is converged!" << std::endl;
              print(solver_data.u);
            }
        }

      // Debugging: exit after one iteration.
      // break;
    } // end while(true)

  // Clean up
  nlopt_destroy(opt);

  return 0;
}


double myfunc(unsigned n, const double * x, double * grad, void * my_func_data)
{
  auto solver_data = static_cast<SolverData *>(my_func_data);

  // Get references to residual, Jacobian, and solution vector stored on the
  // SolverData object.
  Matrix<mpfr_class> & jac = solver_data->jac;
  std::vector<mpfr_class> & u = solver_data->u;
  std::vector<mpfr_class> & r = solver_data->r;

  // Copy data. It would be better if we didn't have to do this. One idea
  // would be to make the SolverData class templated on numerical type...
  u.assign(x, x+n);
  solver_data->ro3.residual_and_jacobian(&r, grad ? &jac : nullptr, u);

  // The gradient is grad_f = J^T * r
  if (grad)
    {
      std::vector<mpfr_class> grad_f =
        jac.matvec_transpose(r);
      for (unsigned int i=0; i<n; ++i)
        grad[i] = grad_f[i].get_d();
    }

  // The objective function value is 0.5 * dot(r,r)
  return 0.5 * dot(r,r).get_d();
}



double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
  // x-index of orbit to which are applying the constraint.
  auto kptr = static_cast<unsigned int *>(data);
  auto k = *kptr;

  if (grad)
    {
      grad[k] = 1.;
      grad[k+1] = 1.;
    }
  return x[k] + x[k+1] - 1;
}


void usage()
{
  std::cout << "Usage: ./drivers/nlopt_test -d 10 -c 0 -v 1 -e 3 - g 5" << std::endl;
  std::cout << "-d = polynomial degree of exactness" << std::endl;
  std::cout << "-c = number of centroid orbits" << std::endl;
  std::cout << "-v = number of vertex orbits" << std::endl;
  std::cout << "-e = number of edge orbits" << std::endl;
  std::cout << "-m = number of median orbits" << std::endl;
  std::cout << "-g = number of general orbits" << std::endl;
  std::cout << "-r nc,nv,ne,nm,ng = comma separated list specifying "
            << "quantity of each orbit" << std::endl;
}


std::string nlopt_result_to_string(nlopt_result res)
{
  switch (res)
    {
    case 1: // generic success code
      return std::string("NLOPT_SUCCESS");
    case 2:
      return std::string("NLOPT_STOPVAL_REACHED");
    case 3:
      return std::string("NLOPT_FTOL_REACHED");
    case 4:
      return std::string("NLOPT_XTOL_REACHED");
    case 5:
      return std::string("NLOPT_MAXEVAL_REACHED");
    case 6:
      return std::string("NLOPT_MAXTIME_REACHED");
    default:
      return std::string("Unrecognized nlopt_result " + std::to_string(res));
    }
}

#endif // HAVE_NLOPT
