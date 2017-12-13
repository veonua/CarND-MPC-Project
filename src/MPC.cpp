#include "MPC.h"
#include <cppad/cppad.hpp>

#ifndef HAVE_CSTDDEF
#define HAVE_CSTDDEF 1
#endif

#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"

using CppAD::AD;

///  Set the timestep length and duration
const size_t N = 5;
const double dt = 0.10;

// set reference speed in m/s
// set to 40.234 m/s = 90 mph
const double ref_v = 20.234;

#define K_CTE  50
#define K_EPSI 1000
#define K_V 100
#define K_DEL 100
#define K_A 1
#define K_DEL_DIFF 1000
#define K_A_DIFF 100

// Solver takes all state and actuator variables in a singular vector.
// Set up indeces for easier access later
const size_t x_start = 0;
const size_t y_start = x_start + N;
const size_t psi_start = y_start + N;
const size_t v_start = psi_start + N;
const size_t cte_start = v_start + N;
const size_t epsi_start = cte_start + N;
const size_t delta_start = epsi_start + N;
const size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    // `fg` is a vector containing the cost and constraints.
    // `vars` is a vector containing the variable values (state & actuators).
    void operator()(ADvector& fg, const ADvector& vars) {
      // The cost is stored is the first element of `fg`.
      // Any additions to the cost should be added to `fg[0]`.
      fg[0] = 0;

        // penalize for cross track error, orientation error and not keeping ref velocity
        for (size_t i = 0; i < N; i++) {
            fg[0] += K_CTE  * CppAD::pow(vars[cte_start + i], 2);
            fg[0] += K_EPSI * CppAD::pow(vars[epsi_start + i], 2);
            fg[0] += K_V    * CppAD::pow(vars[v_start + i] - ref_v, 2);
        }
        // penalize the use of actuators
        for (size_t i = 0; i < N-1; i++) {
            fg[0] += K_DEL * CppAD::pow(vars[delta_start + i], 2);
            fg[0] += K_A   * CppAD::pow(vars[a_start + i], 2);
        }
        // penalize big value gaps in sequential actuations
        for (size_t i = 0; i < N-2; i++) {
            fg[0] += K_DEL_DIFF * CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
            fg[0] += K_A_DIFF   * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
        }

      // Initial constraints
      //
      // We add 1 to each of the starting indices due to cost being located at
      // index 0 of `fg`.
      // This bumps up the position of all the other values.
      fg[1 + x_start] = vars[x_start];
      fg[1 + y_start] = vars[y_start];
      fg[1 + psi_start] = vars[psi_start];
      fg[1 + v_start] = vars[v_start];
      fg[1 + cte_start] = vars[cte_start];
      fg[1 + epsi_start] = vars[epsi_start];

      // The rest of the constraints
      for (size_t t = 1; t < N; t++) {
        // The state at time t+1 .
        AD<double> x1 = vars[x_start + t];
        AD<double> y1 = vars[y_start + t];
        AD<double> psi1 = vars[psi_start + t];
        AD<double> v1 = vars[v_start + t];
        AD<double> cte1 = vars[cte_start + t];
        AD<double> epsi1 = vars[epsi_start + t];

        // The state at time t.
        AD<double> x0 = vars[x_start + t - 1];
        AD<double> y0 = vars[y_start + t - 1];
        AD<double> psi0 = vars[psi_start + t - 1];
        AD<double> v0 = vars[v_start + t - 1];
        AD<double> cte0 = vars[cte_start + t - 1];
        AD<double> epsi0 = vars[epsi_start + t - 1];

        // Only consider the actuation at time t.
        AD<double> delta0 = vars[delta_start + t - 1];
        AD<double> a0 = vars[a_start + t - 1];

        AD<double> f0;
        AD<double> psides0;

          f0 = coeffs[0] +
               coeffs[1] * x0 +
               coeffs[2] * CppAD::pow(x0, 2);

          auto tsides = coeffs[1] + 2 * coeffs[2] * x0;
          // make it possible to switch between 2nd and 3rd degree polynomial fit
          if (coeffs.size() >= 4) {
              f0 += coeffs[3] * CppAD::pow(x0, 3);
              tsides+=3 * coeffs[3] * CppAD::pow(x0,2);
          }
          psides0 = CppAD::atan(tsides);

        // Here's `x` to get you started.
        // The idea here is to constraint this value to be 0.
        //
        // Recall the equations for the model:
        // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
        // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
        // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
        // v_[t+1] = v[t] + a[t] * dt
        // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
        // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
        fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
        fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
        fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / MPC::Lf * dt);
        fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
        fg[1 + cte_start + t] =
                cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
        fg[1 + epsi_start + t] =
                epsi1 - ((psi0 - psides0) + v0 * delta0 / MPC::Lf * dt);
      }
    }
};

//
// MPC class definition
//

MPC::MPC() = default;
MPC::~MPC() = default;

vector<double> MPC::Solve(const Eigen::VectorXd& x0, const Eigen::VectorXd& coeffs) {
    typedef CPPAD_TESTVECTOR(double) Dvector;

    const double x = x0[0];
    const double y = x0[1];
    const double psi = x0[2];
    const double v = x0[3];
    const double cte = x0[4];
    const double epsi = x0[5];

    // number of independent variables
    // N timesteps == N - 1 actuations
    const size_t n_vars = N * 6 + (N - 1) * 2;
    // Number of constraints
    const size_t n_constraints = N * 6;

    // Initial value of the independent variables.
    // Should be 0 except for the initial values.
    Dvector vars(n_vars);
    for (size_t i = 0; i < n_vars; i++) {
        vars[i] = 0.0;
    }
    // Set the initial variable values
    vars[x_start] = x;
    vars[y_start] = y;
    vars[psi_start] = psi;
    vars[v_start] = v;
    vars[cte_start] = cte;
    vars[epsi_start] = epsi;

    // Lower and upper limits for x
    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);

    // Set all non-actuators upper and lowerlimits
    // to the max negative and positive values.
    for (size_t i = 0; i < delta_start; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }

    // The upper and lower limits of delta are set to -25 and 25
    // degrees (values in radians).
    // NOTE: Feel free to change this to something else.
    for (size_t i = delta_start; i < a_start; i++) {
        vars_lowerbound[i] = -0.436332;
        vars_upperbound[i] = 0.436332;
    }

    // Acceleration/decceleration upper and lower limits.
    // NOTE: Feel free to change this to something else.
    for (size_t i = a_start; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
    }

    // Lower and upper limits for constraints
    // All of these should be 0 except the initial
    // state indices.
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (size_t i = 0; i < n_constraints; i++) {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }
    constraints_lowerbound[x_start] = x;
    constraints_lowerbound[y_start] = y;
    constraints_lowerbound[psi_start] = psi;
    constraints_lowerbound[v_start] = v;
    constraints_lowerbound[cte_start] = cte;
    constraints_lowerbound[epsi_start] = epsi;

    constraints_upperbound[x_start] = x;
    constraints_upperbound[y_start] = y;
    constraints_upperbound[psi_start] = psi;
    constraints_upperbound[v_start] = v;
    constraints_upperbound[cte_start] = cte;
    constraints_upperbound[epsi_start] = epsi;

    // Object that computes objective and constraints
    FG_eval fg_eval(coeffs);

    // options
    std::string options;
    options += "Integer print_level  0\n";
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";

    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval>(
            options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
            constraints_upperbound, fg_eval, solution);

    //
    // Check some of the solution values
    //
    bool ok = solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    auto cost = solution.obj_value;
    std::cout << "Cost " << cost << " ok:"<< ok << std::endl;

    vector<double> result;
    result.push_back(solution.x[delta_start]);
    result.push_back(solution.x[a_start]);

    for (size_t i = 0; i < N-1; i++) {
        result.push_back(solution.x[x_start + i]);
        result.push_back(solution.x[y_start + i]);
    }
    return result;
}