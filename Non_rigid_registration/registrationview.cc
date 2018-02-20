//=============================================================================
//
//   Code framework for the lecture
//
//   "Surface Representation and Geometric Modeling"
//
//   Mark Pauly, Mario Botsch, Balint Miklos, and Hao Li
//
//   Copyright (C) 2007 by  Applied Geometry Group and
//                          Computer Graphics Laboratory, ETH Zurich
//
//-----------------------------------------------------------------------------
//
//                                License
//
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License
//   as published by the Free Software Foundation; either version 2
//   of the License, or (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin Street, Fifth Floor,
//   Boston, MA  02110-1301, USA.
//
//=============================================================================

//#include <geodesic/geodesic_algorithm_exact.h>
#include <cstdio>
#include <fstream>
#include "RegistrationViewer.hh"
using namespace std;

static bool getFileNames(char *name, std::vector<std::string> &filenames) {
  ifstream fin(name);
  if (!fin) return false;
  // put all lines to filenames
  std::string line;
  while (getline(fin, line)) {
    filenames.push_back(line);
  }
  fin.close();
  return true;
}
int main(int argc, char **argv) {
  glutInit(&argc, argv);
  RegistrationViewer window("Registration Viewer", 512, 512);

  if (argc == 2) {
    std::vector<std::string> filenames;
    if (!getFileNames(argv[1], filenames)) {
      printf("Cannot open file %s.\n", argv[1]);
      return 1;
    }

    if (window.init(filenames)) {
      glutMainLoop();
    } else {
      printf("Could not init\n");
    }
  } else {
    printf("Usage: %s <meshes.txt>\n", argv[0]);
  }
}

// #include <ceres/ceres.h>
// //#include <glog/logging.h>
// using ceres::AutoDiffCostFunction;
// using ceres::CostFunction;
// using ceres::Problem;
// using ceres::Solve;
// using ceres::Solver;
// // A templated cost functor that implements the residual r = 10 -
// // x. The method operator() is templated so that we can then use an
// // automatic differentiation wrapper around it to generate its
// // derivatives.
// struct F1 {
//   template <typename T>
//   bool operator()(const T* const x1,const T* const x2, T* residual) const {
//     residual[0] = x1[0] + (T)10 * x2[0];
//     return true;
//   }
// };
// struct F2 {
//   template <typename T>
//   bool operator()(const T* const x3,const T* const x4, T* residual) const {
//     residual[0] = T(sqrt(5.0)) * (x3[0]-x4[0]);
//     return true;
//   }
// };
// struct F3 {
//   template <typename T>
//   bool operator()(const T* const x2,const T* const x3, T* residual) const {
//     residual[0] = (x2[0]- (T)2 * x3[0]) * (x2[0] - (T)2 * x3[0]);
//     return true;
//   }
// };
// struct F4 {
//   template <typename T>
//   bool operator()(const T* const x1,const T* const x4, T* residual) const {
//     residual[0] = T(sqrt(10.0)) * (x1[0]-x4[0]) * (x1[0]-x4[0]);
//     return true;
//   }
// };
// int main(int argc, char** argv) {
//   //google::InitGoogleLogging(argv[0]);
//   // The variable to solve for with its initial value. It will be
//   // mutated in place by the solver.
//   double x1 =  3.0; double x2 = -1.0; double x3 =  0.0; double x4 = 1.0;
//   // Build the problem.
//   Problem problem;
//   // Set up the only cost function (also known as residual). This uses
//   // auto-differentiation to obtain the derivative (jacobian).
//   problem.AddResidualBlock(
//   new AutoDiffCostFunction<F1, 1, 1, 1>(new F1), NULL, &x1, &x2);
//   problem.AddResidualBlock(
//   new AutoDiffCostFunction<F2, 1, 1, 1>(new F2), NULL, &x3, &x4);
//   problem.AddResidualBlock(
//   new AutoDiffCostFunction<F3, 1, 1, 1>(new F3), NULL, &x2, &x3);
//   problem.AddResidualBlock(
//   new AutoDiffCostFunction<F4, 1, 1, 1>(new F4), NULL, &x1, &x4);
//   // Run the solver!
//   Solver::Options options;
//   options.minimizer_progress_to_stdout = true;
//   Solver::Summary summary;
//   Solve(options, &problem, &summary);
//   std::cout << summary.BriefReport() << "\n";
//   std::cout << "x1 : " << x1 << " x2 : " << x2<<" x3 : " << x3 <<"  x4 : " <<
//   x4<< "\n"; return 0;
// }

// #include <ceres/ceres.h>
// #include "glog/logging.h"
// // Data generated using the following octave code.
// //   randn('seed', 23497);
// //   m = 0.3;
// //   c = 0.1;
// //   x=[0:0.075:5];
// //   y = exp(m * x + c);
// //   noise = randn(size(x)) * 0.2;
// //   outlier_noise = rand(size(x)) < 0.05;
// //   y_observed = y + noise + outlier_noise;
// //   data = [x', y_observed'];
// const int kNumObservations = 67;
// const double data[] = {
//     0.000000e+00, 1.133898e+00, 7.500000e-02, 1.334902e+00,
//     1.500000e-01, 1.213546e+00, 2.250000e-01, 1.252016e+00,
//     3.000000e-01, 1.392265e+00, 3.750000e-01, 1.314458e+00,
//     4.500000e-01, 1.472541e+00, 5.250000e-01, 1.536218e+00,
//     6.000000e-01, 1.355679e+00, 6.750000e-01, 1.463566e+00,
//     7.500000e-01, 1.490201e+00, 8.250000e-01, 1.658699e+00,
//     9.000000e-01, 1.067574e+00, 9.750000e-01, 1.464629e+00,
//     1.050000e+00, 1.402653e+00, 1.125000e+00, 1.713141e+00,
//     1.200000e+00, 1.527021e+00, 1.275000e+00, 1.702632e+00,
//     1.350000e+00, 1.423899e+00, 1.425000e+00, 5.543078e+00,  // Outlier point
//     1.500000e+00, 5.664015e+00,                              // Outlier point
//     1.575000e+00, 1.732484e+00, 1.650000e+00, 1.543296e+00,
//     1.725000e+00, 1.959523e+00, 1.800000e+00, 1.685132e+00,
//     1.875000e+00, 1.951791e+00, 1.950000e+00, 2.095346e+00,
//     2.025000e+00, 2.361460e+00, 2.100000e+00, 2.169119e+00,
//     2.175000e+00, 2.061745e+00, 2.250000e+00, 2.178641e+00,
//     2.325000e+00, 2.104346e+00, 2.400000e+00, 2.584470e+00,
//     2.475000e+00, 1.914158e+00, 2.550000e+00, 2.368375e+00,
//     2.625000e+00, 2.686125e+00, 2.700000e+00, 2.712395e+00,
//     2.775000e+00, 2.499511e+00, 2.850000e+00, 2.558897e+00,
//     2.925000e+00, 2.309154e+00, 3.000000e+00, 2.869503e+00,
//     3.075000e+00, 3.116645e+00, 3.150000e+00, 3.094907e+00,
//     3.225000e+00, 2.471759e+00, 3.300000e+00, 3.017131e+00,
//     3.375000e+00, 3.232381e+00, 3.450000e+00, 2.944596e+00,
//     3.525000e+00, 3.385343e+00, 3.600000e+00, 3.199826e+00,
//     3.675000e+00, 3.423039e+00, 3.750000e+00, 3.621552e+00,
//     3.825000e+00, 3.559255e+00, 3.900000e+00, 3.530713e+00,
//     3.975000e+00, 3.561766e+00, 4.050000e+00, 3.544574e+00,
//     4.125000e+00, 3.867945e+00, 4.200000e+00, 4.049776e+00,
//     4.275000e+00, 3.885601e+00, 4.350000e+00, 4.110505e+00,
//     4.425000e+00, 4.345320e+00, 4.500000e+00, 4.161241e+00,
//     4.575000e+00, 4.363407e+00, 4.650000e+00, 4.161576e+00,
//     4.725000e+00, 4.619728e+00, 4.800000e+00, 4.737410e+00,
//     4.875000e+00, 4.727863e+00, 4.950000e+00, 4.669206e+00};
// using ceres::AutoDiffCostFunction;
// using ceres::CauchyLoss;
// using ceres::CostFunction;
// using ceres::Problem;
// using ceres::Solve;
// using ceres::Solver;
// struct ExponentialResidual {
//   ExponentialResidual(double x, double y) : x_(x), y_(y) {}
//   template <typename T>
//   bool operator()(const T* const m, const T* const c, T* residual) const {
//     residual[0] = y_ - exp(m[0] * x_ + c[0]);
//     return true;
//   }

//  private:
//   const double x_;
//   const double y_;
// };
// int main(int argc, char** argv) {
//   //google::InitGoogleLogging(argv[0]);
//   double m = 0.0;
//   double c = 0.0;
//   Problem problem;
//   for (int i = 0; i < kNumObservations; ++i) {
//     CostFunction* cost_function =
//         new AutoDiffCostFunction<ExponentialResidual, 1, 1, 1>(
//             new ExponentialResidual(data[2 * i], data[2 * i + 1]));
//     problem.AddResidualBlock(cost_function, new CauchyLoss(0.5), &m, &c);
//   }
//   Solver::Options options;
//   options.linear_solver_type = ceres::DENSE_QR;
//   options.minimizer_progress_to_stdout = true;
//   Solver::Summary summary;
//   Solve(options, &problem, &summary);
//   std::cout << summary.BriefReport() << "\n";
//   std::cout << "Initial m: " << 0.0 << " c: " << 0.0 << "\n";
//   std::cout << "Final   m: " << m << " c: " << c << "\n";
//   return 0;
// }