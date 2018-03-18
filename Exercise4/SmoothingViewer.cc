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
//  CLASS SmoothingViewer - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

#include "SmoothingViewer.hh"

//== IMPLEMENTATION ==========================================================

SmoothingViewer::SmoothingViewer(const char* _title, int _width, int _height)
    : QualityViewer(_title, _width, _height) {
  mesh_.add_property(vpos_);
}

//-----------------------------------------------------------------------------

void SmoothingViewer::keyboard(int key, int x, int y) {
  switch (toupper(key)) {
    case 'N': {
      std::cout << "10 Laplace-Beltrami smoothing iterations: " << std::flush;
      smooth(10);
      calc_weights();
      calc_mean_curvature();
      calc_uniform_mean_curvature();
      calc_gauss_curvature();
      calc_triangle_quality();
      face_color_coding();

      glutPostRedisplay();
      std::cout << "done\n";
      break;
    }
    case 'U': {
      std::cout << "10 uniform smoothing iterations: " << std::flush;
      uniform_smooth(10);
      calc_weights();
      calc_mean_curvature();
      calc_uniform_mean_curvature();
      calc_gauss_curvature();
      calc_triangle_quality();
      face_color_coding();

      glutPostRedisplay();
      std::cout << "done\n";
      break;
    }

    default: {
      QualityViewer::keyboard(key, x, y);
      break;
    }
  }
}

//-----------------------------------------------------------------------------

void SmoothingViewer::smooth(unsigned int _iters) {
  Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
  Mesh::HalfedgeHandle h;
  Mesh::EdgeHandle e;
  Mesh::VertexVertexIter vv_it;
  Mesh::Point laplace(0.0, 0.0, 0.0);
  Mesh::Scalar w, ww;

  // ------------- IMPLEMENT HERE ---------
  // TASK 4.3.b Smoothing using the Laplace-Beltrami.
  // Use eweight_ properties for the individual edge weights
  // and their sum for the normalization term.
  // ------------- IMPLEMENT HERE ---------
  for(unsigned int i = 0; i < _iters; ++i){
      calc_mean_curvature();
      // update position
      for(v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it){
          mesh_.point(*v_it) += 0.5 * mesh_.property(vL_B, *v_it);
      }
      //Hint: do not forget to update normals after 
      //vertex coordinates change.
      mesh_.update_normals();
  }
}

//-----------------------------------------------------------------------------

void SmoothingViewer::uniform_smooth(unsigned int _iters) {
  Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
  //Mesh::VertexVertexIter vv_it;
  Mesh::Point Lu(0.0, 0.0, 0.0);
  //unsigned w;

  // ------------- IMPLEMENT HERE ---------
  // TASK 4.1.b Smoothing using the uniform Laplacian approximation
  // ------------- IMPLEMENT HERE ---------
  for (unsigned int i = 0; i < _iters; ++i) {
      // Firstly calculate uniform mean curvature for each vertice
      calc_uniform_mean_curvature();
      // update position
      for(v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it){
          mesh_.point(*v_it) += 0.5 * mesh_.property(vLu_, *v_it);
      }
      //Hint: do not forget to update normals after 
      //vertex coordinates change.
      mesh_.update_normals();
  }
}

//=============================================================================
