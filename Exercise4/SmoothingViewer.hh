//=============================================================================
//                                                                            
//   Example code for the full-day course
//
//   M. Botsch, M. Pauly, C. Roessl, S. Bischoff, L. Kobbelt,
//   "Geometric Modeling Based on Triangle Meshes"
//   held at SIGGRAPH 2006, Boston, and Eurographics 2006, Vienna.
//
//   Copyright (C) 2006 by  Computer Graphics Laboratory, ETH Zurich, 
//                      and Computer Graphics Group,      RWTH Aachen
//
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
//  CLASS SmoothingViewer
//
//=============================================================================

#ifndef SMOOTHING_VIEWER_HH
#define SMOOTHING_VIEWER_HH

//== INCLUDES =================================================================

#include "QualityViewer.hh"

//== CLASS DEFINITION =========================================================

class SmoothingViewer : public QualityViewer
{
public:
   
    /// default constructor
    SmoothingViewer(const char* _title, int _width, int _height);

    /// iterative Laplacian smoothing
    void smooth(unsigned int _iters);

    void uniform_smooth(unsigned int _iters);

private:

    virtual void keyboard(int key, int x, int y);


    // easier access to new vertex positions
    Mesh::Point& new_pos(Mesh::VertexHandle _vh) 
    { return mesh_.property(vpos_, _vh); }

private:

    OpenMesh::VPropHandleT<Mesh::Point>   vpos_;
};

//=============================================================================
#endif // SMOOTHING_VIEWER_HH defined
//=============================================================================

