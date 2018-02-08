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
//  CLASS Registration - ClosestPoint
//
//=============================================================================

#include "ClosestPoint.hh"
#include "ANN/pr_queue_k.h"

ClosestPoint::
ClosestPoint()
{
    dataPoints_ = NULL;
    kDTree_ = NULL;
}


ClosestPoint::
~ClosestPoint()
{
    release();
}


void
ClosestPoint::
release()
{
    if( kDTree_ ) {
        delete kDTree_;
        kDTree_ = NULL;
    }
    if( dataPoints_ ) {
        annDeallocPts(*dataPoints_);
        delete dataPoints_;
        dataPoints_ = NULL;
    }
}

void
ClosestPoint::
init(
    const std::vector< Vector3d > & _pts
)
{
    release();

    dataPoints_ = new ANNpointArray;
    *dataPoints_ = annAllocPts(_pts.size(),3); // allocate data points

    for(int i = 0; i < (int) _pts.size(); i++) {
        (*dataPoints_)[i][0] = _pts[i][0];
        (*dataPoints_)[i][1] = _pts[i][1];
        (*dataPoints_)[i][2] = _pts[i][2];
    }

    kDTree_ = new ANNkd_tree(       // build search structure
            *dataPoints_,                   // the data points
            _pts.size(),    // number of points
            3               // dimension of space
    );

}

int     // returns index
ClosestPoint::
getClosestPoint(
        const Vector3d & _queryVertex
)
{

    // initialize ANN types for wrapping
    ANNpoint queryPt;                           // query point
    ANNidx nnIdx[1];                            // near neighbor indices
    ANNdist dists[1];                           // near neighbor distances

    queryPt = annAllocPt(3);                    // allocate query point

    // assign values to ANN query point

    queryPt[0] = _queryVertex[0];
    queryPt[1] = _queryVertex[1];
    queryPt[2] = _queryVertex[2];

    ANNmin_k mink;
    int ptsVisited=0;

    kDTree_->ann1Search(            // search
        &queryPt,               // query point
        nnIdx,                  // nearest neighbors (returned)
        dists,                  // distance (returned)
        0,
        &mink,
        ptsVisited);                        // epsilon error bound


    annDeallocPt(queryPt);

    return nnIdx[0];
}
