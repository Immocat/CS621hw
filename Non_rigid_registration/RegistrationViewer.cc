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
//  CLASS RegistrationViewer - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================
#include "RegistrationViewer.hh"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
#include "ClosestPoint.hh"
//#include "Registration.hh"
#include "gl.hh"
//== IMPLEMENTATION ==========================================================
// template <typename Elem>
// class randomaccesstable {
//  public:
//   randomaccesstable(std::size_t initial_size)
//       : data_(initial_size), count_(0) {}

//   inline randomaccesstable &push_back(const Elem &elem) {
//     if (count_ < data_.size())
//       data_[count_++] = elem;
//     else {
//       data_.push_back(elem);
//       ++count_;
//     }
//     return *this;
//   }

//   inline randomaccesstable &remove(const std::size_t index) {
//     if (index < count_) {
//       std::swap(data_[index], data_[count_ - 1]);
//       --count_;
//     }
//     return *this;
//   }

//   inline const Elem &operator[](const std::size_t index) const {
//     return data_[index];
//   }

//   inline Elem &operator[](const std::size_t index) { return data_[index]; }

//   inline std::size_t size() const { return count_; }

//  public:
//   std::vector<Elem> data_;
//   std::size_t count_;
// };
//=============================================================================
RegistrationViewer::RegistrationViewer(const char *_title, int _width,
                                       int _height)
    : GlutExaminer(_title, _width, _height) {
  // clear_draw_modes();

  S_id = 0;

  mode_ = VIEW;

  draw_DG = true;

  draw_fineAlign_intermediate = false;

  draw_S = true;

  draw_M = true;
}

//-----------------------------------------------------------------------------

RegistrationViewer::~RegistrationViewer() {
  // wait for all thread
  // for(std::thread & thread:thread_pool){
  //   if(thread.joinable()){
  //     thread.join();
  //   }
  // }
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool RegistrationViewer::build_target_bvh() {
  // try to build bvh for target

  S_BVH = nanort::BVHAccel<float>();
  nanort::TriangleMesh<float> triangle_mesh(
      (const float *)S.points(), S_indices.data(), sizeof(float) * 3);
  nanort::TriangleSAHPred<float> triangle_pred(
      (const float *)S.points(), S_indices.data(), sizeof(float) * 3);
  nanort::BVHBuildOptions<float> build_options;  // Use default option
  if (!S_BVH.Build(S.n_faces(), triangle_mesh, triangle_pred, build_options)) {
    printf("\t[Build Target BVH]: build BVH failed\n");
    return false;
  }
  // print bvh info
  nanort::BVHBuildStatistics stats = S_BVH.GetStatistics();
  printf(
      "\t[BVH statistics]: %d leaf nodes, %d branch nodes, max tree depth %d\n",
      stats.num_leaf_nodes, stats.num_branch_nodes, stats.max_tree_depth);
}

bool RegistrationViewer::init(const std::vector<std::string> &_filenames) {
  if (_filenames.size() < 2) {
    printf("Need at least 2 meshes for registration.\n");
    return false;
  }
  S_id = 0;
  // copy all filenames
  filenames = _filenames;
  // init EventList
  m_eventList.init(filenames.size() - 1);
  // open Mesh M,S0
  if (!loadTargetMesh(M, filenames[S_id++], M_indices)) {
    printf("Cannot load mesh %s.\n", filenames[S_id - 1].c_str());
    return false;
  }
  if (!loadTargetMesh(S, filenames[S_id++], S_indices)) {
    printf("Cannot load mesh %s.\n", filenames[S_id - 1].c_str());
    return false;
  }
  // improve Mesh
  improveMesh(M, &m_eventList, M_indices);
  improveMesh(S, nullptr, S_indices);
  // build BVH for S
  // build_target_bvh();

  // get average vertex radius for rednering sampled point
  averageVertexDistance_ = get_average_vertex_distance(M);

  // init secene
  // calculate bounding box of all points
  Mesh::Point bbMin(1e9, 1e9, 1e9), bbMax(-1e9, -1e9, -1e9);
  Mesh::ConstVertexIter v_it(M.vertices_begin()), v_end(M.vertices_end());

  Mesh::Point gravity(0, 0, 0);
  for (; v_it != v_end; ++v_it) {
    bbMin.minimize(M.point(*v_it));
    bbMax.maximize(M.point(*v_it));
    gravity += M.point(*v_it);
  }
  m_bbDiagnol = (bbMax - bbMin).length();
  if (M.n_vertices() != 0) gravity /= float(M.n_vertices());
  set_scene(gravity, 0.3 * (bbMin - bbMax).norm());

  // M is at frame 0, M0, M_i+1 = F(frame_events_i,M_i)
  glutPostRedisplay();
  //
  // for (int i = 0; i < (int)_filenames.size(); i++) {
  //   Mesh mesh;
  //   mesh.request_vertex_status();
  //   mesh.request_edge_status();
  //   mesh.request_face_status();
  //   // mesh.request_face_normals();
  //   // mesh.request_vertex_normals();

  //   // load mesh
  //   if (OpenMesh::IO::read_mesh(mesh, _filenames[i].c_str())) {
  //     // clean mes
  //     clean_mesh(mesh);
  //     // mesh.update_normals();
  //     EventList event_list;
  //     improveMesh(mesh, event_list);
  //     // mesh.request_vertex_status();
  //     // mesh.request_edge_status();
  //     // mesh.request_face_status();
  //     // mesh.request_face_normals();
  //     // mesh.request_vertex_normals();
  //     Mesh::ConstVertexIter v_begin(mesh.vertices_begin()),
  //         v_end(mesh.vertices_end()), v_it;
  //     Mesh::Point gravity(0, 0, 0);

  //     // calculate center of gravity
  //     for (v_it = v_begin; v_it != v_end; ++v_it) {
  //       gravity += mesh.point(*v_it);
  //     }

  //     gravity /= float(mesh.n_vertices());

  //     // move to center of gravity
  //     for (v_it = v_begin; v_it != v_end; ++v_it) {
  //       mesh.point(*v_it) -= gravity;
  //     }

  //     mesh.request_face_normals();
  //     mesh.request_vertex_normals();
  //     // compute face & vertex normals
  //     mesh.update_normals();

  //     clean_mesh(mesh);
  //     // compute average vertex distance
  //     averageVertexDistance_ = get_average_vertex_distance(mesh);

  //     // info
  //     std::cerr << _filenames[i] << ": " << mesh.n_vertices() << " vertices,
  //     "
  //               << mesh.n_faces() << " faces\n";
  //   } else {
  //     success = false;
  //   }

  //   meshes_.push_back(mesh);
  //   transformations_.push_back(Transformation());
  // }

  // // calculate bounding box of all points
  // Mesh::Point gravity(0, 0, 0);
  // float count = 0;
  // Mesh::Point bbMin(1e9, 1e9, 1e9), bbMax(-1e9, -1e9, -1e9);
  // for (int i = 0; i < (int)meshes_.size(); i++) {
  //   Mesh::ConstVertexIter v_it(meshes_[i].vertices_begin()),
  //       v_end(meshes_[i].vertices_end());

  //   for (; v_it != v_end; ++v_it) {
  //     bbMin.minimize(meshes_[i].point(v_it));
  //     bbMax.maximize(meshes_[i].point(v_it));
  //     gravity[0] += meshes_[i].point(v_it)[0];
  //     gravity[1] += meshes_[i].point(v_it)[1];
  //     gravity[2] += meshes_[i].point(v_it)[2];
  //     count++;
  //   }
  // }

  // set_scene(Vec3f(0, 0, 0), 0.3 * (bbMin - bbMax).norm());

  // if (success) {
  //   // update face indices for faster rendering
  //   update_face_indices();

  //   glutPostRedisplay();
  // }

  // M_id = 0;
  // S_id = 1;

  return true;
}

bool RegistrationViewer::loadTargetMesh(Mesh &mesh, const std::string &filename,
                                        std::vector<unsigned int> &indices) {
  mesh.clear();
  mesh.request_vertex_status();
  mesh.request_edge_status();
  mesh.request_face_status();
  mesh.request_face_normals();
  mesh.request_vertex_normals();

  // load mesh
  if (!OpenMesh::IO::read_mesh(mesh, filename.c_str())) {
    return false;
  }

  mesh.update_normals();
  update_face_indices(mesh, indices);
  // load info
  std::cerr << "[Load Mesh]: " << filename << ": " << mesh.n_vertices()
            << " vertices, " << mesh.n_faces() << " faces\n";

  // transformations_.push_back(Transformation());
  return true;
}

//-----------------------------------------------------------------------------

void RegistrationViewer::update_face_indices(
    const Mesh &mesh, std::vector<unsigned int> &indices) {
  indices.clear();

  Mesh::ConstFaceIter f_it(mesh.faces_sbegin()), f_end(mesh.faces_end());
  Mesh::ConstFaceVertexIter fv_it;

  indices.reserve(mesh.n_faces() * 3);

  for (; f_it != f_end; ++f_it) {
    indices.push_back((*(fv_it = mesh.cfv_iter(*f_it))).idx());
    indices.push_back((*(++fv_it)).idx());
    indices.push_back((*(++fv_it)).idx());
  }
}

//-----------------------------------------------------------------------------
void RegistrationViewer::start_morph(bool bFine_align) {
  //
}

// TODO: need more parameters
// calculate signed distance
// void calculateSignedDistance(const Mesh& mesh);

// TODO: need more parameters
// void constranTopology(Mesh* mesh);

// save M to Disk
bool RegistrationViewer::saveMeshToDisk(const Mesh &mesh) const {
  //
}

void RegistrationViewer::draw(const std::string &_draw_mode) {
  if (M_indices.size() == 0 || S_indices.size() == 0) {
    GlutExaminer::draw(_draw_mode);
    return;
  }
  if (draw_M) drawM(M, M_indices, OpenMesh::Vec3f(0.1, 0.5, 0.1), _draw_mode);
  if (draw_S) drawS(S, S_indices, OpenMesh::Vec3f(0.5, 0.5, 0.5), _draw_mode);

  // display deformation graph
  if (draw_DG) {
    // std::lock_guard<std::mutex> guard(M_DG.mutex);

    // draw sumplepoints
    glEnable(GL_COLOR_MATERIAL);
    // orange
    glColor3f(1, 0.64453125, 0);
    std::vector<Vector3d> X_transed;
    X_transed.reserve(M_DG.X.size());
    for (int i = 0; i < M_DG.X.size(); ++i) {
      X_transed.emplace_back(M_DG.X_T[i].transformPoint(M_DG.X[i]));
    }
    for (int i = 0; i < (int)M_DG.X.size(); ++i) {
      glPushMatrix();
      const Vector3d &pt(X_transed[i]);
      // Vector3d pt(M_DG.X_T[i].transformPoint(M_DG.X[i]));
      glTranslatef(pt[0], pt[1], pt[2]);
      glutSolidSphere(averageVertexDistance_ * 0.5, 10, 10);
      glPopMatrix();
    }
    // TODO: draw DG edges
    glLineWidth(averageVertexDistance_ * 0.5);
    glBegin(GL_LINES);
    for (int i = 0; i < (int)M_DG.X.size(); ++i) {
      const Vector3d &pt(X_transed[i]);
      for (const int &j : M_DG.X_edges[i]) {
        glVertex3f(pt[0], pt[1], pt[2]);
        glVertex3f(X_transed[j][0], X_transed[j][1], X_transed[j][2]);
      }
    }
    glEnd();
    glDisable(GL_COLOR_MATERIAL);
  }
}

//-----------------------------------------------------------------------------

void RegistrationViewer::drawM(const Mesh &mesh,
                               const std::vector<unsigned int> &indices,
                               const OpenMesh::Vec3f &color,
                               const std::string &_draw_mode) {
  glPushMatrix();
  // apply transformation matrix of scan
  // transformations_[index].apply_gl();
  if (_draw_mode == "Wireframe") {
    //
    std::lock_guard<std::mutex> lock(draw_finealign_mutex);

    glDisable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glColor3f(color[0], color[1], color[2]);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    if (draw_fineAlign_intermediate) {
      std::vector<Vector3d> transedPos;
      std::vector<Vector3d> transedNormal;
      transedPos.reserve(mesh.n_vertices());
      transedNormal.reserve(mesh.n_vertices());

      Mesh::ConstVertexIter v_it(mesh.vertices_begin()),
          v_end(mesh.vertices_end());
      // transform
      for (; v_it != v_end; ++v_it) {
        Vec3f p = mesh.point(*v_it);
        Vec3f n = mesh.normal(*v_it);
        const Transformation &trans(mesh.property(M_fineAlignTrans, *v_it));
        transedPos.emplace_back(
            trans.transformPoint(Vector3d(p[0], p[1], p[2])));
        transedNormal.emplace_back(
            trans.transformVector(Vector3d(n[0], n[1], n[2])));
      }
      // draw

      ::glVertexPointer(3, GL_DOUBLE, 0, transedPos.data());
      ::glNormalPointer(GL_DOUBLE, 0, transedNormal.data());
      glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT,
                     &(indices[0]));

    } else {
      GL::glVertexPointer(mesh.points());
      GL::glNormalPointer(mesh.vertex_normals());
      glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT,
                     &(indices[0]));
    }
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_COLOR_MATERIAL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  } else if (_draw_mode == "Solid Smooth") {
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);
    glShadeModel(GL_SMOOTH);
    glColor3f(color[0], color[1], color[2]);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    std::lock_guard<std::mutex> lock(draw_finealign_mutex);

    if (draw_fineAlign_intermediate) {
      //
      std::vector<Vector3d> transedPos;
      std::vector<Vector3d> transedNormal;
      transedPos.reserve(mesh.n_vertices());
      transedNormal.reserve(mesh.n_vertices());

      Mesh::ConstVertexIter v_it(mesh.vertices_begin()),
          v_end(mesh.vertices_end());
      // transform
      for (; v_it != v_end; ++v_it) {
        Vec3f p = mesh.point(*v_it);
        Vec3f n = mesh.normal(*v_it);
        const Transformation &trans(mesh.property(M_fineAlignTrans, *v_it));
        transedPos.emplace_back(
            trans.transformPoint(Vector3d(p[0], p[1], p[2])));
        transedNormal.emplace_back(
            trans.transformVector(Vector3d(n[0], n[1], n[2])));
      }
      // draw

      ::glVertexPointer(3, GL_DOUBLE, 0, transedPos.data());
      ::glNormalPointer(GL_DOUBLE, 0, transedNormal.data());
      glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT,
                     &(indices[0]));
    } else {
      GL::glVertexPointer(mesh.points());
      GL::glNormalPointer(mesh.vertex_normals());
      glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT,
                     &(indices[0]));
    }
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_COLOR_MATERIAL);
  }
  glPopMatrix();

  // draw pairs
  if (draw_fineAlign_intermediate) {
    // std::lock_guard<std::mutex> guard(M_DG.mutex);

    // draw sumplepoints
    glEnable(GL_COLOR_MATERIAL);
    // orange
    glColor3f(1, 0.64453125, 0);

    // TODO: draw DG edges
    glLineWidth(averageVertexDistance_ * 0.5);
    glBegin(GL_LINES);
    // for (int i = 0; i < (int)M_DG.X.size(); ++i) {
    //   const Vector3d &pt(X_transed[i]);
    //   for (const int &j : M_DG.X_edges[i]) {
    //     glVertex3f(pt[0], pt[1], pt[2]);
    //     glVertex3f(X_transed[j][0], X_transed[j][1], X_transed[j][2]);
    //   }
    // }
    Mesh::ConstVertexIter v_it(mesh.vertices_begin()),
        v_end(mesh.vertices_end());
    for (; v_it != v_end; ++v_it) {
      // pruned
      if (!mesh.property(hasTarget, *v_it)) continue;

      const OpenMesh::Vec3f &p = mesh.point(*v_it);
      const Transformation &trans = mesh.property(M_fineAlignTrans, *v_it);
      const Vector3d &ci = mesh.property(targetPoints, *v_it);
      Vector3d pos(p[0], p[1], p[2]);
      pos = trans.transformPoint(pos);
      glVertex3f(pos[0], pos[1], pos[2]);
      glVertex3f(ci[0], ci[1], ci[2]);
    }

    glEnd();
    glDisable(GL_COLOR_MATERIAL);
  }
}
void RegistrationViewer::drawS(const Mesh &mesh,
                               const std::vector<unsigned int> &indices,
                               const OpenMesh::Vec3f &color,
                               const std::string &_draw_mode) {
  glPushMatrix();
  // apply transformation matrix of scan
  // transformations_[index].apply_gl();
  if (_draw_mode == "Wireframe") {
    //
    glDisable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glColor3f(color[0], color[1], color[2]);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    {
      GL::glVertexPointer(mesh.points());
      GL::glNormalPointer(mesh.vertex_normals());
      glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT,
                     &(indices[0]));
    }
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_COLOR_MATERIAL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  } else if (_draw_mode == "Solid Smooth") {
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);
    glShadeModel(GL_SMOOTH);
    glColor3f(color[0], color[1], color[2]);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    {
      GL::glVertexPointer(mesh.points());
      GL::glNormalPointer(mesh.vertex_normals());
      glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT,
                     &(indices[0]));
    }
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_COLOR_MATERIAL);
  }

  glPopMatrix();
}
//-----------------------------------------------------------------------------

void RegistrationViewer::keyboard(int key, int x, int y) {
  switch (key) {
    case ' ': {
      std::cout << "Start Morphing with fine alignment..." << std::endl;
      start_morph(true);
      glutPostRedisplay();
      break;
    }
    case 'i': {
      std::cout << "Improve Mesh M..." << std::endl;
      improveMesh(M, &m_eventList, M_indices);
      glutPostRedisplay();
      break;
    }
    case 'n': {
      // sampledPoints_.clear();
      if (S_id >= filenames.size()) {
        printf("All target processed\n");
        break;
      }
      if (!loadTargetMesh(S, filenames[S_id++], S_indices)) {
        printf("Cannot load mesh %s.\n", filenames[S_id - 1].c_str());
        break;
      }
      improveMesh(S, nullptr, S_indices);
      glutPostRedisplay();
      break;
    }
    case 's': {
      if (!saveMeshToDisk(M)) {
        printf("Cannot save M to disk.\n");
        break;
      }
    }
    case ';': {
      draw_S = !draw_S;
      glutPostRedisplay();
      break;
    }
    case '.': {
      draw_M = !draw_M;
      glutPostRedisplay();
      break;
    }
    case 'c': {
      thread_pool.emplace_back(
          [&]() { coarseNonRigidAlignment(&M, M_indices, S); });
      // coarseNonRigidAlignment(&M, M_indices, S);
      glutPostRedisplay();
      break;
    }
    case 'f': {
       thread_pool.emplace_back([&]() { fineLinearAlignment(M, S, S_indices);
       });
      //fineLinearAlignment(M, S, S_indices);
      glutPostRedisplay();
      break;
    }
    case 't': {
      // TODO: do topologyChange
      glutPostRedisplay();
      break;
    }
    case 'd': {
      draw_DG = !draw_DG;
      printf("[Draw Deformation Graph]: %s\n", draw_DG ? "True" : "False");
      glutPostRedisplay();
      break;
    }
    default: {
      GlutExaminer::keyboard(key, x, y);
      break;
    }
  }
}

//-----------------------------------------------------------------------------

void RegistrationViewer::mouse(int button, int state, int x, int y) {
  // manual object transformation when pressing SHIFT
  if (glutGetModifiers() & GLUT_ACTIVE_SHIFT)
    mode_ = VIEW;
  else
    mode_ = VIEW;

  GlutExaminer::mouse(button, state, x, y);
}

//-----------------------------------------------------------------------------

/// called during mouse motion while button is pressed
void RegistrationViewer::motion(int x, int y) {
  switch (mode_) {
    default:
    case VIEW: {
      GlutExaminer::motion(x, y);
      break;
    }
      // case MOVE: {
      //   // manual object transformation when pressing SHIFT

      //   // zoom
      //   if (button_down_[0] && button_down_[1]) {
      //     float dy = y - last_point_2D_[1];
      //     float h = height_;
      //     Transformation mv_tr = Transformation::retrieve_gl();
      //     Transformation tr(0.0, 0.0, radius_ * dy * 3.0 / h);
      //     transformations_[currIndex_] =
      //         mv_tr.inverse() * tr * mv_tr * transformations_[currIndex_];
      //   }
      //   // rotation
      //   else if (button_down_[0]) {
      //     if (last_point_ok_) {
      //       Vec2i new_point_2D;
      //       Vec3f new_point_3D;
      //       bool new_point_ok;

      //       new_point_2D = Vec2i(x, y);
      //       new_point_ok = map_to_sphere(new_point_2D, new_point_3D);

      //       if (new_point_ok) {
      //         Vec3f axis = (last_point_3D_ % new_point_3D);
      //         float cos_angle = (last_point_3D_ | new_point_3D);

      //         if (fabs(cos_angle) < 1.0) {
      //           float angle = 2.0 * acos(cos_angle);
      //           // rotate(axis, angle);
      //           Transformation mv_tr = Transformation::retrieve_gl();
      //           mv_tr.translation_.fill(0);
      //           Transformation tr(angle, Vector3f(axis[0], axis[1],
      //           axis[2]));

      //           transformations_[currIndex_] =
      //               mv_tr.inverse() * tr * mv_tr *
      //               transformations_[currIndex_];
      //         }
      //       }
      //     }
      //   }
      //   // translation
      //   else if (button_down_[1]) {
      //     float dx = x - last_point_2D_[0];
      //     float dy = y - last_point_2D_[1];

      //     float z =
      //         -((modelview_matrix_[2] * center_[0] +
      //            modelview_matrix_[6] * center_[1] +
      //            modelview_matrix_[10] * center_[2] + modelview_matrix_[14])
      //            /
      //           (modelview_matrix_[3] * center_[0] +
      //            modelview_matrix_[7] * center_[1] +
      //            modelview_matrix_[11] * center_[2] +
      //            modelview_matrix_[15]));

      //     float aspect = (float)width_ / (float)height_;
      //     float up = tan(fovy_ / 2.0f * M_PI / 180.f) * near_;
      //     float right = aspect * up;

      //     Transformation mv_tr = Transformation::retrieve_gl();
      //     Transformation tr(2.0 * dx / width_ * right / near_ * z,
      //                       -2.0 * dy / height_ * up / near_ * z, 0.0f);
      //     transformations_[currIndex_] =
      //         mv_tr.inverse() * tr * mv_tr * transformations_[currIndex_];
      //   }

      //   // remeber points
      //   last_point_2D_ = Vec2i(x, y);
      //   last_point_ok_ = map_to_sphere(last_point_2D_, last_point_3D_);

      //   glutPostRedisplay();
      //   break;
      // }
  }
}

//=============================================================================

// /// save current points
// void RegistrationViewer::save_points() {
//   std::ofstream out(outputFilename_.c_str());
//   for (int i = 0; i < numProcessed_; i++) {
//     // get points on target meshes
//     std::vector<Vector3d> pts = get_points(meshes_[i]);
//     std::vector<Vector3d> normals = get_normals(meshes_[i]);

//     // transform using current scan transformations
//     pts = transformations_[i].transformPoints(pts);
//     normals = transformations_[i].transformVectors(normals);

//     for (int j = 0; j < (int)pts.size(); j++) {
//       out << "v " << pts[j][0] << " " << pts[j][1] << " " << pts[j][2] << "
//       "; out << "vn " << normals[j][0] << " " << normals[j][1] << " "
//           << normals[j][2] << std::endl;
//     }
//   }
//   out.close();
//   std::cout << "merged points saved to: " << outputFilename_ << std::endl;
// }

//=============================================================================

/// clean mesh by removing "bad" triangles
// void RegistrationViewer::clean_mesh(Mesh &_mesh) {
//   Mesh::FaceIter f_it = _mesh.faces_begin();
//   for (; f_it != _mesh.faces_end(); ++f_it) {
//     float maxEdge = 0;
//     float minEdge = 1e9;
//     Mesh::ConstFaceHalfedgeIter fh_it = _mesh.fh_iter(f_it.handle());
//     while (fh_it) {
//       OpenMesh::Vec3f p =
//       _mesh.point(_mesh.from_vertex_handle(fh_it.handle())); OpenMesh::Vec3f
//       q = _mesh.point(_mesh.to_vertex_handle(fh_it.handle())); float
//       edgeLength = sqrt((p - q) | (p - q)); maxEdge = std::max(maxEdge,
//       edgeLength); minEdge = std::min(minEdge, edgeLength);

//       ++fh_it;
//     }
//     if (minEdge / maxEdge < 0.2) _mesh.delete_face(f_it.handle(), true);
//   }
//   _mesh.garbage_collection();
// }

//=============================================================================

/// get average vertex distance
float RegistrationViewer::get_average_vertex_distance(const Mesh &_mesh) const {
  float accDist = 0;
  int accCount = 0;

  Mesh::ConstHalfedgeIter he_it = _mesh.halfedges_begin();
  for (; he_it != _mesh.halfedges_end(); ++he_it) {
    OpenMesh::Vec3f p = _mesh.point(_mesh.from_vertex_handle(*he_it));
    OpenMesh::Vec3f q = _mesh.point(_mesh.to_vertex_handle(*he_it));
    float edgeLength = sqrt((p - q) | (p - q));
    accDist += edgeLength;
    accCount++;
  }

  if (accCount > 0)
    return accDist / float(accCount);
  else
    return 0;
}

//=============================================================================

// get points of mesh
void RegistrationViewer::get_points(
    const Mesh &_mesh, std::vector<Vector3d> &pts,
    std::vector<OpenMesh::VertexHandle> &vHandles) const {
  pts.clear();
  vHandles.clear();
  Mesh::ConstVertexIter v_it(_mesh.vertices_begin()),
      v_end(_mesh.vertices_end());
  for (; v_it != v_end; ++v_it) {
    Vec3f p = _mesh.point(*v_it);
    pts.push_back(Vector3d(p[0], p[1], p[2]));
    vHandles.push_back(*v_it);
  }
}

//=============================================================================

// get normals of mesh
void RegistrationViewer::get_normals(const Mesh &_mesh,
                                     std::vector<Vector3d> &normals) {
  normals.clear();
  Mesh::ConstVertexIter v_it(_mesh.vertices_begin()),
      v_end(_mesh.vertices_end());
  for (; v_it != v_end; ++v_it) {
    Vec3f n = _mesh.normal(*v_it);
    normals.push_back(Vector3d(n[0], n[1], n[2]));
  }
}

//=============================================================================

// get border vertices of mesh
void RegistrationViewer::get_borders(const Mesh &_mesh,
                                     std::vector<bool> &borders) {
  borders.clear();
  Mesh::ConstVertexIter v_it(_mesh.vertices_begin()),
      v_end(_mesh.vertices_end());
  for (; v_it != v_end; ++v_it) {
    borders.push_back(_mesh.is_boundary(*v_it));
  }
}

//=============================================================================

// /// perform registration
// void RegistrationViewer::start_morph(bool bFine_align) {
//   std::vector<Vector3d> src;
//   std::vector<Vector3d> target;
//   std::vector<Vector3d> target_normals;

//   // calculate correspondences
//   calculate_correspondences(src, target, target_normals);

//   Registration reg;
//   printf("Num correspondences: %d\n", int(src.size()));

//   // calculate optimal transformation
//   Transformation opt_tr;
//   if (_tangential_motion) {
//     opt_tr = reg.register_point2surface(src, target, target_normals);
//   } else {
//     opt_tr = reg.register_point2point(src, target);
//   }

//   // set transformation
//   transformations_[currIndex_] = opt_tr * transformations_[currIndex_];
// }
// static bool sample_valid(const Vector3d &p, const std::vector<int> &ids,
//                          const std::vector<Vector3d> &_pts, float len2) {
//   for (int j = 0; j < ids.size(); ++j) {
//     if (length2(p - _pts[ids[j]]) < len2) return false;
//   }
//   return true;
// }
// //=============================================================================
// /// subsample points
// std::vector<int> RegistrationViewer::subsample(
//     const std::vector<Vector3d> &_pts) {
//   std::vector<int> indeces;

//   float subsampleRadius = 5 * averageVertexDistance_;

//   // EXERCISE 2.2
//   /////////////////////////////////////////////////////////////
//   // subsampling:
//   // perform uniform subsampling by storing the indeces of the subsampled
//   points
//   // in 'indeces' (Hint: take advantage of the fact that consecutive points
//   in
//   // the
//   //  vector _pts are also often close in the scan )
//   ////////////////////////////////////////////////////////////////////////////
//   // always add the first point...
//   // similar to poisson disk sampling, except that the surounding is choosed
//   // from nearby [-15,15] points
//   if (_pts.size() == 0) {
//     printf("[ERROR] subsample : no source point to choose!\n");
//   } else {
//     // suppose we need 1000 samples
//     float subsampleRadius2 = subsampleRadius * subsampleRadius;
//     int di = _pts.size() / 2000;
//     for (int i = 0; i < _pts.size(); i += di) {
//       while (i < _pts.size() &&
//              !sample_valid(_pts[i], indeces, _pts, subsampleRadius2))
//         ++i;
//       indeces.push_back(i);
//     }

//     // float subsampleRadius2 = subsampleRadius * subsampleRadius;
//     // std::unordered_set<int> indeces_set;
//     // while(indeces.size() < 1000){
//     //     int i = rand() % _pts.size();
//     //     if(indeces_set.find(i)!= indeces_set.end()) continue;
//     //     bool is_valid = true;
//     //     for(int j = 0; j < indeces.size(); ++j){
//     //         if(length2(_pts[i] - _pts[j]) < subsampleRadius2){
//     //             is_valid = false;
//     //             break;
//     //         }
//     //     }
//     //     if(!is_valid) continue;
//     //     // put in
//     //     indeces_set.insert(i);
//     //     indeces.push_back(i);
//     // }

//     // randomaccesstable<int> candidates(_pts.size());
//     // std::unordered_set<int> indeces_set;
//     // //1. random choose one, push it to ouput indeces and candidates
//     // int start_point_id = rand() % _pts.size();
//     // candidates.push_back(start_point_id);
//     // indeces.push_back(start_point_id);
//     // indeces_set.insert(start_point_id);
//     // while(candidates.size() > 0){
//     //     // randomly choose a point in candidate
//     //     int candidate = candidates[0];
//     //     candidates.remove(0);
//     //     int c_max = std::min(candidate + 15,(int)_pts.size());
//     //     // to keep it simple, generate
//     //     for(int i = std::max(0,candidate-15); i < c_max;++i){
//     //         // new points should not be chosen before
//     //         if(indeces_set.find(i)!= indeces_set.end()) continue;
//     //         // new points should be farther than the threshold
//     //         if(length2(_pts[i] - _pts[candidate]) < subsampleRadius2)
//     //         continue;
//     //         //
//     //         bool is_valid = true;
//     //         for(int j = 0; j < indeces.size(); ++j){
//     //             if(length2(_pts[i] - _pts[j]) < subsampleRadius2){
//     //                 is_valid = false;
//     //                 break;
//     //             }
//     //         }
//     //         if(!is_valid) continue;
//     //         // new point is valid, put it in
//     //         indeces_set.insert(i);
//     //         indeces.push_back(i);
//     //         candidates.push_back(i);
//     //     }
//     // }
//   }

//   ////////////////////////////////////////////////////////////////////////////

//   // keep indeces/samples for display
//   sampledPoints_ = indeces;
//   printf("subsample: choose %d samples\n", indeces.size());
//   return indeces;
// }

// //=============================================================================

// /// calculate correspondences
// void RegistrationViewer::calculate_correspondences(
//     std::vector<Vector3d> &_src, std::vector<Vector3d> &_target,
//     std::vector<Vector3d> &_target_normals) {
//   _src.clear();
//   _target.clear();
//   _target_normals.clear();

//   std::vector<Vector3d> srcCandidatePts;
//   std::vector<Vector3d> srcCandidateNormals;
//   std::vector<Vector3d> targetCandidatePts;
//   std::vector<Vector3d> targetCandidateNormals;
//   std::vector<double> src_target_dis2;

//   // get points on src mesh
//   std::vector<Vector3d> srcPts = get_points(meshes_[currIndex_]);
//   std::vector<Vector3d> srcNormals = get_normals(meshes_[currIndex_]);

//   // transform using current scan transformation
//   srcPts = transformations_[currIndex_].transformPoints(srcPts);
//   srcNormals = transformations_[currIndex_].transformVectors(srcNormals);

//   // subsample the points
//   std::vector<int> indeces = subsample(srcPts);

//   // iterate over all previously processed scans and find correspondences
//   // note that we perform registration to all other scans simultaneously, not
//   // only pair-wise
//   for (int i = 0; i < numProcessed_; i++) {
//     if (i == currIndex_) continue;

//     // get points on target meshes
//     std::vector<Vector3d> targetPts = get_points(meshes_[i]);
//     std::vector<Vector3d> targetNormals = get_normals(meshes_[i]);
//     std::vector<bool> targetBorders = get_borders(meshes_[i]);

//     // transform using current scan transformation
//     targetPts = transformations_[i].transformPoints(targetPts);
//     targetNormals = transformations_[i].transformVectors(targetNormals);

//     // set up closest point structure
//     ClosestPoint cp;
//     cp.init(targetPts);

//     // find closest points for each src vertex
//     for (int j = 0; j < (int)indeces.size(); j++) {
//       int index = indeces[j];

//       int bestIndex = cp.getClosestPoint(srcPts[index]);

//       // do not keep border correspondences
//       if (!targetBorders[bestIndex]) {
//         srcCandidatePts.push_back(srcPts[index]);
//         srcCandidateNormals.push_back(srcNormals[index]);
//         targetCandidatePts.push_back(targetPts[bestIndex]);
//         targetCandidateNormals.push_back(targetNormals[bestIndex]);
//         src_target_dis2.push_back(
//             length2(srcPts[index] - targetPts[bestIndex]));
//       }
//     }
//   }

//   printf("calculate_correspondences: candidate num: %d\n",
//          srcCandidatePts.size());

//   // EXERCISE 2.3
//   /////////////////////////////////////////////////////////////
//   // correspondence pruning:
//   // prune correspondence based on
//   // - distance threshold
//   // - normal compatability
//   //
//   // fill _src, _target, and _target_normals from the candidate pairs

//   // normals of correspondences do not deviate more than 60 degrees
//   float normalCompatabilityThresh = 60;
//   // distance threshold is 3 times the median distance
//   float distMedianThresh = 3;

//   ////////////////////////////////////////////////////////////////////////////
//   std::vector<double> dis2_copy = src_target_dis2;
//   std::nth_element(src_target_dis2.begin(),
//                    src_target_dis2.begin() + src_target_dis2.size() / 2,
//                    src_target_dis2.end());
//   double distThreshold =
//       distMedianThresh *
//       std::sqrt(
//           std::max((double)0, src_target_dis2[src_target_dis2.size() / 2]));
//   double distThreshold2 = distThreshold * distThreshold;
//   double cosineThresh = std::cos(normalCompatabilityThresh * M_PI / 180.0);
//   for (int i = 0; i < srcCandidatePts.size(); ++i) {
//     if (dis2_copy[i] > distThreshold2) continue;
//     // suppose distMedian is always less than 90 degree, which means that
//     // cosineThresh<=cos(theta)<=1
//     if (dot_product(srcCandidateNormals[i].normalize(),
//                     targetCandidateNormals[i].normalize()) < cosineThresh)
//       continue;
//     // put into output
//     _src.push_back(srcCandidatePts[i]);
//     _target.push_back(targetCandidatePts[i]);
//     _target_normals.push_back(targetCandidateNormals[i]);
//   }
//   printf("calculate_correspondences: final num: %d\n", (int)_src.size());

//   ////////////////////////////////////////////////////////////////////////////
// }

//=============================================================================
