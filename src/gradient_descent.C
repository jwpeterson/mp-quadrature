#include "solver_data.h"
#include "gradient_descent.h"
#include "vect.h"

bool gradient_descent(SolverData & solver_data)
{
  // Gradient descent parameters.
  mpfr_class tol = solver_data.tol;
  mpfr_class divtol = solver_data.divtol;
  unsigned int maxits = solver_data.maxits;
  // const mpfr_class alphamin(1.e-10);
  bool do_backtracking = solver_data.do_backtracking;
  unsigned iter = 0;
  bool converged = false;

  ResidualAndJacobian & residual_and_jacobian =
    solver_data.residual_and_jacobian;

  CheckFeasibility & check_feasibility =
    solver_data.check_feasibility;

  std::vector<mpfr_class> & u = solver_data.u;

  // Problem size
  unsigned int n = u.size();

  // Storage needed for algorithm.
  std::vector<mpfr_class> r(n), du(n), grad_f(n), grad_f_old(n), u_old(n), trial_u(n);
  Matrix<mpfr_class> jac(n,n);

  while (true)
    {
      ++iter;
      if (iter > maxits)
        break;

      // Compute the residual and Jacobian at the current guess.
      residual_and_jacobian(&r, &jac, u);

      // Compute the size of the scalar function "f"
      // which we are trying to minimize, _not_ dot(r,r)^0.5.
      mpfr_class residual = 0.5 * dot(r,r);

      // Debugging:
      if (solver_data.verbose)
        std::cout << "Iteration " << iter << ", residual = " << residual << std::endl;

      // If the residual is small enough, return true.
      if (residual < tol)
        {
          converged = true;
          break;
        }

      // If the residual is too large, just give up.
      if (residual > divtol)
        break;

      // Compute grad(f) = jac^T * r at the current guess.
      grad_f = jac.matvec_transpose(r);

      // Debugging
      // std::cout << "grad_f=" << std::endl;
      // print(grad_f);

      // Compute scalar coefficient gamma that tells how far down the gradient direction to
      // travel.
      mpfr_class gamma_n(1);
      if (iter > 1)
        {
          std::vector<mpfr_class> z = u - u_old;
          std::vector<mpfr_class> y = grad_f - grad_f_old;

          gamma_n = abs(dot(z, y)) / dot(y, y);

          // Debugging:
          // std::cout << "gamma_n = " << gamma_n << std::endl;
        }

      // Store old information for next iteration. Note: do this before updating u.
      grad_f_old = grad_f;
      u_old = u;

      // If the solution becomes infeasible, reduce gamma_n by half
      // until it's feasible again.
      if (do_backtracking)
        {
          for (unsigned int back=0; back<50; ++back)
            {
              trial_u = u - gamma_n * grad_f;
              bool feasible = check_feasibility(trial_u);

              if (!feasible)
                {
                  // Debugging:
                  // std::cout << "gradient_descent: reducing gamma_n for feasibility" << std::endl;
                  gamma_n /= 2;
                }
              else
                break;
            }
        }

      // Compute the new iterate.
      subtract_scaled(u, gamma_n, grad_f);
    }

  return converged;
}
