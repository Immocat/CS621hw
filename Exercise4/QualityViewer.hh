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
//  CLASS QualityViewer
//
//=============================================================================

#ifndef QUALVIEWERWIDGET_HH
#define QUALVIEWERWIDGET_HH

//== INCLUDES =================================================================

#include "MeshViewer.hh"

//== CLASS DEFINITION =========================================================

class QualityViewer : public MeshViewer
{
public:
   
    /// default constructor
    QualityViewer(const char* _title, int _width, int _height);

    // destructor
    ~QualityViewer();

    /// open mesh
    virtual bool open_mesh(const char* _filename);

protected:

    typedef OpenMesh::VPropHandleT<Mesh::Scalar> Vertex_property;
    std::vector<float>  face_colors_;

    virtual void init();
    virtual void draw(const std::string& _draw_mode);


    /// calculate vertex and edge weights
    void calc_weights();

    /// calculate mean curvature per vertex
    void calc_mean_curvature();
    void calc_uniform_mean_curvature();

    void calc_gauss_curvature();

    /// calculate triangle shape indices
    void calc_triangle_quality();
    void face_color_coding();

    void find_min_max(Vertex_property prop, Mesh::Scalar& min, Mesh::Scalar& max);
    Mesh::Color value_to_color(float value, float min, float max);

    /// set vertex color from vertex curvature
    void color_coding(Vertex_property prop);


    OpenMesh::VPropHandleT<Mesh::Scalar>  vweight_, vunicurvature_, vcurvature_, vgausscurvature_;
    OpenMesh::VPropHandleT<Mesh::Point> vLu_, vL_B;
    OpenMesh::EPropHandleT<Mesh::Scalar>  eweight_;
    OpenMesh::FPropHandleT<Mesh::Scalar>  tshape_;

    GLuint  textureID_;
};


//=============================================================================
#endif // QUALVIEWERWIDGET_HH defined
//=============================================================================

