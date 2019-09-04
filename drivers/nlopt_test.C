// C++ includes
#include <iostream>

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
#include <math.h>
#include <stdlib.h> // random()

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
  auto solver_data = static_cast<SolverData *>(my_func_data);

  // Copy data. It would be better if we didn't have to do this. One idea
  // would be to make the SolverData class templated on numerical type...
  std::vector<mpfr_class> r;
  Matrix<mpfr_class> jac;
  std::vector<mpfr_class> u(x, x+n);
  solver_data->residual_and_jacobian(&r, grad ? &jac : nullptr, u);

  // The gradient is grad_f = J^T * r
  if (grad)
    {
      std::vector<mpfr_class> grad_f = jac.matvec_transpose(r);
      for (unsigned int i=0; i<n; ++i)
        grad[i] = grad_f[i].get_d();
    }

  // The objective function value is 0.5 * dot(r,r)
  return 0.5 * dot(r,r).get_d();
}


// Constraints are of the form fc(x) <= 0, where fc is required to have
// the same parameters as the objective function. In our case, there is
// one inequality constraint associated with each ro3 orbit of the form:
// 1 - x - y >= 0,
// which we multiply by -1 to write as:
// x + y - 1 <= 0
// So in this case, the incoming "data" pointer just needs to tell us
// the index of the x dof for the orbit being constrained.
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

// Convert a numerical nlopt_result successful return code to a human readable string.
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

int main()
{
  // You can't trust all these digits from doubles, but you can with
  // mpfr_class objects.
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  // It's faster to comment this out while debugging.
  // mpfr_set_default_prec(256);

  // Use the current time since epoch as a seed.
  time_t seed = time(nullptr);
  // long seed = 1567251384;
  std::cout << "seed=" << seed << std::endl;
  srandom(seed);

  // Test making d==3 Ro3 objects
  // unsigned int d=3;
  // Ro3 r(d, /*nc*/1, /*nv*/0, /*ne*/0, /*ng*/1); // 4 QP
  // Ro3 r(d, /*nc*/0, /*nv*/0, /*ne*/2, /*ng*/0); // 6 QP
  // Ro3 r(d, /*nc*/1, /*nv*/1, /*ne*/1, /*ng*/0); // 7 QP

  // d==4, dim==5
  // unsigned int d=4;
  // Ro3 r(d, /*nc*/0, /*nv*/0, /*ne*/1, /*ng*/1); // 6 QP
  // Ro3 r(d, /*nc*/1, /*nv*/1, /*ne*/0, /*ng*/1); // 7 QP
  // Ro3 r(d, /*nc*/1, /*nv*/0, /*ne*/2, /*ng*/0); // 7 QP
  // Ro3 r(d, /*nc*/0, /*nv*/1, /*ne*/2, /*ng*/0); // 9 QP

  // d==7, dim=12
  unsigned int d=7;
  Ro3 r(d, /*nc*/0, /*nv*/0, /*ne*/0, /*ng*/4); // 12 QP

  // d==10, dim=22
  // unsigned int d=10;
  // Ro3 r(d, /*nc*/1, /*nv*/0, /*ne*/0, /*ng*/7); // 22 QP
  // Ro3 r(d, /*nc*/0, /*nv*/1, /*ne*/0, /*ng*/7); // 24 QP
  // Ro3 r(d, /*nc*/0, /*nv*/0, /*ne*/2, /*ng*/6); // 24 QP
  // Ro3 r(d, /*nc*/1, /*nv*/1, /*ne*/1, /*ng*/6); // 25 QP

  // Print information.
  // std::cout << "Rule has " << r.n_qp() << " quadrature points." << std::endl;
  //
  // std::cout << "CENTROID: ["
  //           << r.first_dof(Ro3::CENTROID)
  //           << ","
  //           << r.last_dof(Ro3::CENTROID)
  //           << ")" << std::endl;
  //
  // std::cout << "VERTEX: ["
  //           << r.first_dof(Ro3::VERTEX)
  //           << ","
  //           << r.last_dof(Ro3::VERTEX)
  //           << ")" << std::endl;
  //
  // std::cout << "EDGE: ["
  //           << r.first_dof(Ro3::EDGE)
  //           << ","
  //           << r.last_dof(Ro3::EDGE)
  //           << ")" << std::endl;
  //
  // std::cout << "GENERAL: ["
  //           << r.first_dof(Ro3::GENERAL)
  //           << ","
  //           << r.last_dof(Ro3::GENERAL)
  //           << ")" << std::endl;
  //
  // return 0;

  SolverData solver_data(r);
  solver_data.verbose = true;
  solver_data.maxits = 50;
  solver_data.residual_and_jacobian = ResidualAndJacobian(d);
  solver_data.check_feasibility = CheckFeasibility(d);
  solver_data.do_backtracking = true;
  solver_data.residual_reduction_required = true;

  // "GN" = global derivative-free optimization
  // nlopt_algorithm alg = NLOPT_GN_DIRECT_L;

  // "LD" = local gradient based optimization
  // nlopt_algorithm alg = NLOPT_LD_MMA; // method of moving asymptotes
  // nlopt_algorithm alg = NLOPT_LD_SLSQP; // sequential quadratic programming
  // nlopt_algorithm alg = NLOPT_LD_TNEWTON_PRECOND_RESTART; // Preconditioned truncated Newton
  // nlopt_algorithm alg = NLOPT_LD_LBFGS; // Low storage BFGS

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

  // The problem dimension depends only on "d".
  unsigned int dim = r.dim();

  // Set minimization algorithm and dimensionality
  nlopt_opt opt = nlopt_create(alg, dim);

  if (alg == NLOPT_AUGLAG ||
      alg == NLOPT_G_MLSL)
    {
      // The objective function, bounds, and nonlinear-constraint
      // parameters of local_opt are ignored.
      nlopt_opt local_opt = nlopt_create(NLOPT_LD_SLSQP, dim);
      nlopt_set_xtol_rel(local_opt, 1.e-8);
      nlopt_set_ftol_rel(local_opt, 1.e-16);
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
  nlopt_set_xtol_rel(opt, 1.e-16);
  nlopt_set_ftol_rel(opt, 1.e-16);

  // You can use this parameter to control how long some of the
  // gradient free optimization algorithms run, otherwise it seems
  // that they can just run forever.
  // nlopt_set_maxeval(opt, 9999);

  // Solution vector
  std::vector<double> x;

  while (true)
    {
      // Generate random initial guess.
      r.guess(x);

      // Debugging: Print entries of x
      // std::cout << "initial guess=" << std::endl;
      // for (const auto & val : x)
      //   std::cout << val << std::endl;

      // The d=7 solution with many digits of accuracy is:
      // 2.6517028157436251428754180460739e-2
      // 6.2382265094402118173683000996350e-2
      // 6.7517867073916085442557131050869e-2
      // 4.3881408714446055036769903139288e-2
      // 5.5225456656926611737479190275645e-2
      // 3.2150249385198182266630784919920e-1
      // 2.8775042784981585738445496900219e-2
      // 3.4324302945097146469630642483938e-2
      // 6.6094919618673565761198031019780e-1
      // 6.7493187009802774462697086166421e-2
      // 5.1584233435359177925746338682643e-1
      // 2.7771616697639178256958187139372e-1

      if (d==7)
        x =
          {
            2.65e-02, // w1
            6.23e-02, // x1
            6.75e-02, // y1

            4.38e-02, // w2
            5.52e-02, // x2
            3.21e-01, // y2

            2.87e-02, // w3
            3.43e-02, // x3
            6.60e-01, // y3

            6.74e-02, // w4
            5.15e-01, // x4
            2.77e-01, // y4
          };

      // For the d=10 case, the smallest objective function value
      // (2.4432860276909146e-18) found was the
      // following. Interestingly, it has all non-zero
      // weights. Unfortunately, standard Newton iterations were not
      // able to find a true solution nearby.  By doing some
      // multiprecision Hessian iterations, I was able to reduce the
      // objective function value ever so slightly to
      // 2.277695250723979e-18, but it does in fact seem that it's a
      // local minimum and *not* a zero of the original
      // equations. This is immensely frustrating!
      //
      // After increasing the mpfr default precision, it takes *much*
      // longer to run! But I was able to get slightly smaller
      // objective function values without too many attempts, as shown
      // below. Unfortunately, none of these appears to be an actual
      // root still.

      // objective = 2.16035791285266049414540599285916e-18
      // 4.10631353805066753870356421884935e-02
      // 1.88263162566837394251706427894533e-02
      // 2.66190499925494095112554759907653e-02
      // 6.19350346092532255681817332515493e-01
      // 3.64742436934004515824980785509979e-02
      // 4.07637940402920118110330349736614e-01
      // 4.45710580937403177959055255996645e-01
      // 1.50520631448374259814260156531418e-02
      // 3.44224206753352135468304595633526e-02
      // 8.36719130923662679499841487995582e-01
      // 2.05418424536244169364973544134045e-02
      // 6.04777485241435819318667199695483e-01
      // 3.77007248872530120742396775312955e-02
      // 1.69114072028801500080508191103945e-02
      // 3.88911463098898513290890832649893e-02
      // 1.62309672492393031184576557279797e-01
      // 4.09982936736961664503375857293577e-02
      // 2.05625153604004712315500569275173e-01
      // 1.65768691363852166409387223211525e-01
      // 4.17478844541306539500657990515720e-03
      // 1.60246543761499195268616091425429e-02
      // 3.20127671327287907643288633607881e-02

      // objective = 1.21009979268099643727942979965228e-18
      // 3.95649782304018093892494789542980e-02
      // 2.86047395103584652445438685219870e-03
      // 2.59335946313303082610968175458765e-02
      // 9.87804849985655410160578782097218e-03
      // 1.96949301321903275097113805713889e-02
      // 3.87351313220835824502330524410354e-01
      // 2.79783975644267643434481840358785e-02
      // 1.45133875849392003642979176447625e-02
      // 3.62201620629575421483892228025070e-02
      // 1.12322721643994552498391215067386e-01
      // 1.54199549310062179047209340865265e-02
      // 1.58000632193682960524583336336946e-01
      // 3.31591202109597985026034905331471e-02
      // 4.04065140568377403895716781789815e-02
      // 1.98851122520592171749598264796077e-01
      // 1.68900307402240323906283947508200e-01
      // 3.75996717288595364014902600047208e-02
      // 4.37187661901254132068572744174162e-01
      // 1.48440879464915048702167155170173e-01
      // 2.29834082040305878735786393463059e-02
      // 3.85311556909678915028294454714342e-02
      // 3.10766485177656093252807067983667e-01

      // objective = 2.63171244825735365771036449016650e-19
      // 3.73354919714630334448557391624490e-02 // w
      // 2.25537358367580570428501118840359e-02 // w
      // 3.52659023901809726142531076220621e-02
      // 3.21692820195059514531976674334146e-01
      // 1.59930863210363018167559090443319e-02 // w
      // 1.56149617888017178124471229239134e-01
      // 3.55907295085211433649519108257664e-02
      // 1.45358661539488447089762956920822e-02 // w
      // 3.35310579748630127605579787086754e-02
      // 1.07652804917763383896200934941589e-01
      // 2.19863013990911643907866235281290e-03 // w
      // 2.68687576262772939705847363711655e-02
      // 0.00000000000000000000000000000000e+00
      // 4.00674911522872836955322384255851e-02
      // 1.83468183141260843260766932871775e-01
      // 1.81963881954311457178619093610905e-01
      // 3.79932658791976157752756648733339e-02
      // 4.21281271967661685717843056409038e-01
      // 1.51993172760920947084173349139746e-01
      // 2.08794271919127832903839703249105e-02
      // 3.75003057621806634713834682770539e-01
      // 3.19731089240110205595968295710918e-02

      // if (d==10)
      //   x = {
      //     4.0982937046221364e-02,
      //     1.9798492190714671e-02,
      //     3.3095608903349873e-02,
      //     3.8161345981027001e-01,
      //     1.4994858159817931e-02,
      //     1.0723552400713073e-01,
      //     3.5017307536618969e-02,
      //     3.6839421470875953e-02,
      //     4.4070591601247633e-01,
      //     1.4583079686880440e-01,
      //     1.7605672332677621e-02,
      //     3.8460650957081052e-02,
      //     1.7059525415369317e-01,
      //     2.0591306611632748e-02,
      //     6.4450289076099920e-01,
      //     3.2415668641540085e-01,
      //     4.0452451678720927e-02,
      //     1.7496389690339609e-01,
      //     6.2927469367745104e-01,
      //     2.7234852018090342e-03,
      //     3.3649635977032091e-02,
      //     9.6446038980717408e-01,
      //   };

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

#endif // HAVE_NLOPT
