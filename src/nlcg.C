#include "solver_data.h"
#include "nlcg.h"
#include "vect.h"

bool nlcg(SolverData & solver_data)
{
  // nlcg parameters.
  mpfr_class tol = solver_data.tol;
  mpfr_class divtol = solver_data.divtol;
  unsigned int maxits = solver_data.maxits;
  unsigned iter = 0;
  bool converged = false;

  Ro3 & ro3 = solver_data.ro3;

  std::vector<mpfr_class> & u = solver_data.u;
  std::vector<mpfr_class> & r = solver_data.r;
  Matrix<mpfr_class> & jac = solver_data.jac;

  // Problem size
  unsigned int n = u.size();

  // Storage needed for algorithm.
  std::vector<mpfr_class> du(n), /*u_old(n),*/ du_old(n), grad_f(n), s(n), s_old(n), trial_u(n);

  while (true)
    {
      ++iter;
      if (iter > maxits)
        break;

      // Compute the residual and Jacobian at the current guess.
      ro3.residual_and_jacobian(&r, &jac, u);

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

      // Compute the gradient of the objective function.
      grad_f = jac.matvec_transpose(r);

      // du = -1 * grad_f; // does not work for some reason?
      for (unsigned int i=0; i<n; ++i)
        du[i] = -grad_f[i];

      if (iter == 1)
        {
          // Initially we take a gradient descent step.
          s = du;
        }
      if (iter > 1)
        {
          // Debugging
          // std::cout << "du_old=" << std::endl;
          // print(du_old);
          // std::cout << "du=" << std::endl;
          // print(du);

          // Compute beta using "FR" formula
          mpfr_class beta = dot(du, du) / dot(du_old, du_old);

          // Compute beta using "PR" formula
          // mpfr_class beta = dot(du, du - du_old) / dot(du_old, du_old);

          // Compute beta using "HS" formula
          // mpfr_class beta = -dot(du, du - du_old) / dot(s_old, du - du_old);

          // Debugging:
          // std::cout << "beta = " << beta << std::endl;

          // If beta is negative then we have to set beta to zero and use
          // the s = -grad(f). The "FR" method will not give a negative beta...
          // if (beta < 0)
          //   beta = 0;

          // It does not make sense to do more than n iterations of nlcg,
          // you are supposed to take a gradient descent step (beta=0) before continuing.
          if (iter % n == 0)
            {
              // std::cout << "Taking gradient descent step to 'reset' nlcg." << std::endl;
              beta = 0;
            }

          s = du + beta * s_old;
        }

      // Compute next iterate.
      // TODO: we should actually take a step u += alpha * s, where
      // alpha is chosen by line search.
      mpfr_class alpha = mpfr_class(1);

      for (unsigned int back=0; back<10; ++back)
        {
          // std::cout << "Trying step with alpha=" << alpha << std::endl;
          trial_u = u + alpha * s;
          ro3.residual_and_jacobian(&r, nullptr, u);
          mpfr_class trial_residual = 0.5 * dot(r,r);
          if (trial_residual < residual)
            break;
          else
            alpha /= 2;
        }

      // std::cout << "Accepting step with alpha=" << alpha << std::endl;

      // In gradient descent, the following "alpha" is used:
      // if (iter > 1)
      //   {
      //     std::vector<mpfr_class> z = u - u_old;
      //     std::vector<mpfr_class> y = du_old - du; // grad_f - grad_f_old;
      //     alpha = abs(dot(z, y)) / dot(y, y);
      //   }
      u += alpha * s;

      du_old = du;
      s_old = s;
      // u_old = u;
    } // end while(true)

  return converged;
}

