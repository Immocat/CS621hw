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
//=============================================================================
//
//  CLASS Registration - IMPLEMENTATION
//
//=============================================================================

#include <cstdio>
#include <cstring>

#include <math.h>
#include "Registration.hh"

using namespace std;

//=============================================================================

// point-2-point registration
Transformation Registration::register_point2point(
    const std::vector<Vector3d>& _src,     // moving points (source)
    const std::vector<Vector3d>& _target)  // target points  (target)
{
  // we minimise the following distance measure:
  // min e = sum(i=1.._n)(|| R src[i] + T - target[i] ||^2)

  int numRows = _src.size() * 3;
  double* A = new double[numRows * 6];
  double* b = new double[numRows];
  memset(A, 0, sizeof(double) * numRows * 6);
  memset(b, 0, sizeof(double) * numRows);

  for (int i = 0; i < (int)_src.size(); i++) {
    // EXERCISE 2.4
    // /////////////////////////////////////////////////////////////
    // point-2-point constraints
    // set up matrix A and b with the linear constraints

    ////////////////////////////////////////////////////////////////////////////
    A[i * 3 * 6 + 1] = _src[i].v[2];     // Piz
    A[i * 3 * 6 + 2] = -(_src[i].v[1]);  //-Piy
    A[i * 3 * 6 + 3] = (double)1;

    A[(i * 3 + 1) * 6] = -(_src[i].v[2]);   //-Piz
    A[(i * 3 + 1) * 6 + 2] = _src[i].v[0];  // Pix
    A[(i * 3 + 1) * 6 + 4] = (double)1;

    A[(i * 3 + 2) * 6] = _src[i].v[1];         // Piy
    A[(i * 3 + 2) * 6 + 1] = -(_src[i].v[0]);  //-Pix
    A[(i * 3 + 2) * 6 + 5] = (double)1;        //-Pix

    b[i * 3] = _target[i].v[0] - _src[i].v[0];
    b[i * 3 + 1] = _target[i].v[1] - _src[i].v[1];
    b[i * 3 + 2] = _target[i].v[2] - _src[i].v[2];

    ////////////////////////////////////////////////////////////////////////////
  }

  // solve overdetermined Ax=b
  double x[6];

  Transformation tr;

  if (Solve(A, b, x, numRows)) {
    // get the transformation from the rotation angles and translation vector
    tr.rotation_ = GetRotation(x[0], x[1], x[2]);
    tr.translation_[0] = x[3];
    tr.translation_[1] = x[4];
    tr.translation_[2] = x[5];
  } else {
    printf("Registration::ComputeTransformation() => Cholesky failed\n");
  }

  if (A) delete[] A;
  if (b) delete[] b;

  return tr;
}

//=============================================================================

// point-2-surface registration
Transformation Registration::register_point2surface(
    const std::vector<Vector3d>& _src,  // moving points (source)
    const std::vector<Vector3d>&
        _target,  // points on the tangent plane (target)
    const std::vector<Vector3d>&
        _target_normals  // the normal of the tangent plane
) {
  // we minimise the following distance measure:
  // min e = sum(i=1.._n)(|| nTarget[i] . ( R src[i] + T - target[i] ) ) ||^2)

  int numRows = _src.size();
  double* A = new double[numRows * 6];
  double* b = new double[numRows];
  memset(A, 0, sizeof(double) * numRows * 6);
  memset(b, 0, sizeof(double) * numRows);

  for (int i = 0; i < (int)_src.size(); i++) {
    // EXERCISE 2.5
    // /////////////////////////////////////////////////////////////
    // point-2-surface constraints
    // set up matrix A and b with the linear constraints

    ////////////////////////////////////////////////////////////////////////////
    Vector3d pi_x_ni = cross_product(_src[i], _target_normals[i]);
    A[i * 6] = pi_x_ni.v[0];
    A[i * 6 + 1] = pi_x_ni.v[1];
    A[i * 6 + 2] = pi_x_ni.v[2];
    A[i * 6 + 3] = _target_normals[i].v[0];
    A[i * 6 + 4] = _target_normals[i].v[1];
    A[i * 6 + 5] = _target_normals[i].v[2];

    b[i] = -(dot_product(_src[i] - _target[i], _target_normals[i]));
    ////////////////////////////////////////////////////////////////////////////
  }

  // solve overdetermined Ax=b
  double x[6];

  Transformation tr;

  if (Solve(A, b, x, numRows)) {
    // get the transformation from the rotation angles and translation vector
    tr.rotation_ = GetRotation(x[0], x[1], x[2]);
    tr.translation_[0] = x[3];
    tr.translation_[1] = x[4];
    tr.translation_[2] = x[5];
  } else {
    printf("Registration::ComputeTransformation() => Cholesky failed\n");
  }

  if (A) delete[] A;
  if (b) delete[] b;

  return tr;
}

//=============================================================================
// solve the linear system Ax = b with 6 unknowns
bool Registration::Solve(double* A, double* b, double x[6], int numRows) {
  double AtA[6][6];
  double Atb[6];
  for (int i = 0; i < 6; i++) {
    Atb[i] = 0;
    for (int j = 0; j < 6; j++) {
      AtA[i][j] = 0;
    }
  }

  for (int r = 0; r < numRows; r++) {
    for (int i = 0; i < 6; i++) {
      Atb[i] += A[r * 6 + i] * b[r];
      for (int j = 0; j < 6; j++) {
        AtA[i][j] += A[r * 6 + i] * A[r * 6 + j];
      }
    }
  }

  return CholeskySolve(AtA, Atb, x);
}

//=============================================================================

// Solve x from AtAx=b using Cholesky decomposition.
bool Registration::CholeskySolve(double AtA[6][6], double Atb[6], double x[6]) {
  int i, j, k;
  double sum;

  for (i = 0; i < 6; i++) {
    for (sum = AtA[i][i], k = 0; k < i; k++) {
      sum -= AtA[i][k] * AtA[i][k];
    }

    if (sum < 0.0) {
      // std::cerr << "CRegistration::CholeskySolve() => matrix not
      // pos.semidef." << std::endl;
      return false;
    } else if (sum < 1.0e-7) {
      // std::cerr << "CRegistration::CholeskySolve() => matrix not pos.def." <<
      // std::endl;
      return false;
    } else {
      AtA[i][i] = sqrt(sum);
    }

    for (j = i + 1; j < 6; j++) {
      for (sum = AtA[i][j], k = 0; k < i; k++) {
        sum -= AtA[i][k] * AtA[j][k];
      }

      AtA[j][i] = sum / AtA[i][i];
    }
  }

  for (i = 0; i < 6; i++)  // Forward elimination;
  {
    for (sum = Atb[i], j = 0; j < i; j++)  // solve Ly=b, store y in x
    {
      sum -= AtA[i][j] * x[j];
    }
    x[i] = sum / AtA[i][i];
  }

  for (i = 5; i >= 0; i--)  // Backward elimination;
  {
    for (sum = x[i], j = i + 1; j < 6; j++)  // solve L'x = y
    {
      sum -= AtA[j][i] * x[j];
    }
    x[i] = sum / AtA[i][i];
  }

  return true;
}

//=============================================================================

// returns the rotation matrix for 3 rotation angles
Matrix3x3d Registration::GetRotation(double alpha, double beta, double gamma) {
  Matrix3x3d R;

  double sa = sin(alpha);
  double ca = sqrt(1 - sa * sa);
  double sb = sin(beta);
  double cb = sqrt(1 - sb * sb);
  double sr = sin(gamma);
  double cr = sqrt(1 - sr * sr);

  R[0][0] = cb * cr;
  R[0][1] = -cb * sr;
  R[0][2] = sb;

  R[1][0] = sa * sb * cr + ca * sr;
  R[1][1] = -sa * sb * sr + ca * cr;
  R[1][2] = -sa * cb;

  R[2][0] = -ca * sb * cr + sa * sr;
  R[2][1] = ca * sb * sr + sa * cr;
  R[2][2] = ca * cb;

  return R;
}

//=============================================================================
