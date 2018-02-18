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
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "DeformationGraph.hh"
#include "EventList.hh"
#include "GlutExaminer.hh"
#include "Transformation.hh"
class EventList;
//== CLASS DEFINITION =========================================================

class RegistrationViewer : public GlutExaminer {
  typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

 public:
  /// default constructor
  RegistrationViewer(const char* _title, int _width, int _height);

  // destructor
  ~RegistrationViewer();

  /// copy all filenames, open Mesh S0, set M->S0, improve M
  bool init(const std::vector<std::string>& _filenames);

 protected:
  virtual void draw(const std::string& _draw_mode);
  virtual void keyboard(int key, int x, int y);
  virtual void motion(int x, int y);
  virtual void mouse(int button, int state, int x, int y);

 private:
  /// draw the mesh in scene
  void draw(const Mesh& mesh, const std::vector<unsigned int>& indices,
            const OpenMesh::Vec3f& color, const std::string& _draw_mode);

  /// clean mesh by removing "bad" triangles
  void clean_mesh(Mesh& mesh);

  /// get points of mesh
  void get_points(const Mesh& mesh, std::vector<Vector3d>& points,
                  std::vector<OpenMesh::VertexHandle>& vhandles) const;

  /// get normals of mesh
  std::vector<Vector3d> get_normals(const Mesh& mesh);

  /// get border vertices of mesh
  std::vector<bool> get_borders(const Mesh& mesh);

  /// get average vertex distance
  float get_average_vertex_distance(const Mesh& mesh) const;

  // void init_improve_weights();
  // void new_mid_point(const OpenMesh::TriMesh_ArrayKernelT<>
  // &M,OpenMesh::EdgeHandle heh,OpenMesh::Vec3f &newPoint);

  // improve mesh triangle mesh M with split and collapse, save operation to
  // eventList
  void improveMesh(Mesh& mesh, EventList* eventList,
                   std::vector<unsigned int>& indices);

  // automaticly start the whole process with or without fine alignment
  void start_morph(bool bFine_align = true);

  // Coarse alignment M to S
  void coarseNonRigidAlignment(Mesh* src_mesh,
                               const std::vector<unsigned int>& src_indices,
                               const Mesh& target_mesh);

  // fineLiearAlignment M to S
  void fineLiearAlignment(Mesh& src_mesh, const Mesh& target_mesh);

  // TODO: need more parameters
  // calculate signed distance
  // void calculateSignedDistance(const Mesh& mesh);

  // TODO: need more parameters
  // void constranTopology(Mesh* mesh);

  // save M to Disk
  bool saveMeshToDisk(const Mesh& mesh) const;

  bool loadTargetMesh(Mesh& mesh, const std::string& filename,
                      std::vector<unsigned int>& indices);

  void update_face_indices(const Mesh& mesh,
                           std::vector<unsigned int>& indices);

  void subsample(const std::vector<Vector3d>& pts, std::vector<int>& sample_ids,
                 const Mesh& mesh);

 protected:
  enum Mode { VIEW, MOVE } mode_;

 protected:
  float averageVertexDistance_;
  std::vector<std::string> filenames;
  std::vector<unsigned int> M_indices;
  std::vector<unsigned int> S_indices;
  DeformationGraph M_DG;
  int S_id;
  Mesh M;
  Mesh S;
  EventList m_eventList;
  bool draw_DG;
};

//=============================================================================
#endif  // REGISTRATIONVIEWERWIDGET_HH defined
//=============================================================================
