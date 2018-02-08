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
//  CLASS Registration
//
//=============================================================================

#ifndef _REGISTRATION_HH_
#define _REGISTRATION_HH_

#include <vector>
#include "Transformation.hh"


// point-2-point and point-2-surface registration using linearized rotation matrices
class Registration
{
public:

    // point-2-point registration
    Transformation register_point2point(
        const std::vector< Vector3d > & _src,
        const std::vector< Vector3d > & _target );

    // point-2-surface registration
    Transformation register_point2surface(
        const std::vector< Vector3d > & _src,
        const std::vector< Vector3d > & _target,
        const std::vector< Vector3d > & _target_normals );

private:

    // solve the linear system Ax = b with 6 unknowns
    bool Solve( double * A, double * b, double x[6], int numRows );

    // solves the linear equation AtA x = Atb using cholesky decomposition
    bool CholeskySolve(double AtA[6][6], double Atb[6], double x[6]);

    // returns the rotation matrix for 3 rotation angles
    Matrix3x3d GetRotation(double alpha, double beta, double gamma);

};

#endif

