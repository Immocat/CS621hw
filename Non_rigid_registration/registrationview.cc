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

#include <cstdio>
#include "RegistrationViewer.hh"
#include <fstream>
using namespace std;

static bool getFileNames(char* name,std::vector<std::string> &filenames){
  ifstream fin(name);
  if(!fin) return false;
  //put all lines to filenames
  std::string line;
  while(getline(fin, line)){
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
    if(!getFileNames(argv[1], filenames)){
      printf("Cannot open file %s.\n",argv[1]);
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
// #include <glog/logging.h>
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
//   std::cout << "x1 : " << x1 << " x2 : " << x2<<" x3 : " << x3 <<"  x4 : " << x4<< "\n";
//   return 0;
// }