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
//  CLASS Transformation
//
//=============================================================================

#ifndef TRANSFORMATION_HPP_
#define TRANSFORMATION_HPP_

#include <vector>
#include "Matrix.hh"
#include "Vector.hh"

/**
 * Transformation class
 *
 * contains a rotation and translation defining a rigid transformation
 */
class Transformation {
 public:
  /// constructor: identity transformation
  Transformation();

  /// constructor: translation
  Transformation(float tx, float ty, float tz);

  /// constructor: rotation around axis
  Transformation(float angle, Vector3f axis);

  /// set identity transformation
  void set_identity();

  /// apply transformation to current OpenGL Matrix
  void apply_gl();

  /// retrieve curren OpenGL transformation
  static Transformation retrieve_gl();

  /// concatenate two transformations
  Transformation operator*(const Transformation& o);

  /// return inverse transformation
  Transformation inverse();

  /// Transform point
  Vector3d transformPoint(const Vector3d& p);

  /// Transform vector
  Vector3d transformVector(const Vector3d& v);

  /// Transform points
  std::vector<Vector3d> transformPoints(const std::vector<Vector3d>& ps);

  /// Transform vectors
  std::vector<Vector3d> transformVectors(const std::vector<Vector3d>& vs);

  Matrix3x3d rotation_;
  Vector3d translation_;
};

#endif /* TRANSFORMATION_HPP_ */
