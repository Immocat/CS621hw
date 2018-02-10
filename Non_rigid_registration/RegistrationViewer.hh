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
//  CLASS RegistrationViewer
//
//=============================================================================

#ifndef REGISTRATIONVIEWERWIDGET_HH
#define REGISTRATIONVIEWERWIDGET_HH

//== INCLUDES =================================================================

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "GlutExaminer.hh"
#include "Transformation.hh"

//== CLASS DEFINITION =========================================================

class RegistrationViewer : public GlutExaminer {
  typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

 public:
  /// default constructor
  RegistrationViewer(const char* _title, int _width, int _height);

  // destructor
  ~RegistrationViewer();

  /// set output filename
  void set_output(const std::string& filename);

  /// open meshes
  bool open_meshes(const std::vector<std::string>& _filenames);

 protected:
  virtual void draw(const std::string& _draw_mode);
  virtual void keyboard(int key, int x, int y);
  virtual void motion(int x, int y);
  virtual void mouse(int button, int state, int x, int y);

 private:
  /// update buffer with face indices
  void update_face_indices();

  /// draw the mesh in scene
  virtual void draw(int index, const OpenMesh::Vec3f& color);

  /// save current points
  void save_points();

  /// clean mesh by removing "bad" triangles
  void clean_mesh(Mesh& mesh);

  /// perform registration
  void perform_registration(bool tangential_motion);

  /// subsample points
  std::vector<int> subsample(const std::vector<Vector3d>& pts);

  /// calculate correspondences
  void calculate_correspondences(std::vector<Vector3d>& src,
                                 std::vector<Vector3d>& target,
                                 std::vector<Vector3d>& target_normals);

  /// get points of mesh
  std::vector<Vector3d> get_points(const Mesh& mesh);

  /// get normals of mesh
  std::vector<Vector3d> get_normals(const Mesh& mesh);

  /// get border vertices of mesh
  std::vector<bool> get_borders(const Mesh& mesh);

  /// get average vertex distance
  float get_average_vertex_distance(const Mesh& mesh);

 protected:
  enum Mode { VIEW, MOVE } mode_;

 protected:
  std::string outputFilename_;

  float averageVertexDistance_;
  int currIndex_;
  int numProcessed_;
  std::vector<Mesh> meshes_;
  std::vector<std::vector<unsigned int> > indices_;
  std::vector<Transformation> transformations_;

  std::vector<int> sampledPoints_;
};

//=============================================================================
#endif  // REGISTRATIONVIEWERWIDGET_HH defined
//=============================================================================
