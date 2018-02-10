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
//  CLASS ClosestPoint
//
//=============================================================================

#ifndef CLOSESTPOINT_HPP_
#define CLOSESTPOINT_HPP_

#include <vector>
#include "ANN/ANN.h"
#include "Vector.hh"

/**
 * ClosestPoint class
 *
 * class that allows efficient closest point lookup using KD-tree (ANN library)
 */
class ClosestPoint {
 public:
  /// constructor
  ClosestPoint();

  /// destructor
  ~ClosestPoint();

  /// build ANN KD-tree for _pts
  void init(const std::vector<Vector3d>& _pts);

  /// release data
  void release();

  /// retrieve closest point of query
  int getClosestPoint(const Vector3d& _queryVertex);

 private:
  /// data points only used when ANN search is performed
  ANNpointArray* dataPoints_;

  /// kd tree search structure
  ANNkd_tree* kDTree_;
};

#endif /* CLOSESTPOINT_HPP_ */
