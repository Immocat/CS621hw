#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <unordered_set>
#include <vector>
#include "Transformation.hh"
#include "Vector.hh"
class DeformationGraph {
  typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

 public:
  DeformationGraph() {}
  virtual ~DeformationGraph() {
    m_mesh->remove_property(m_weights);
    m_mesh->remove_property(m_wXids);
  }
  DeformationGraph(const std::vector<Vector3d> &xi_pos, Mesh *mesh,
                   const std::vector<OpenMesh::VertexHandle> &vHandles,
                   const std::vector<int> &sample_ids,
                   const std::vector<unsigned int> &src_indices, int numOfComps,
                   OpenMesh::VPropHandleT<int> comp_id);
  void init(const std::vector<Vector3d> &xi_pos, Mesh *mesh,
            const std::vector<OpenMesh::VertexHandle> &vHandles,
            const std::vector<int> &sample_ids,
            const std::vector<unsigned int> &src_indices, int numOfComps,
            OpenMesh::VPropHandleT<int> comp_id);
  void updateTransforms(const std::vector<Transformation> &newTransforms);
  // void clear(){
  //   X.clear();
  //   X_edges.clear();
  // }
  // public data
  std::vector<Vector3d> X;
  std::vector<Vector3d> X_normal;
  std::vector<std::unordered_set<int>> X_edges;
  std::vector<Transformation> X_T;  // 3x3 transformation + 3x1 translation
  // set by coarseAlignment, only for visualization
  std::vector<Vector3d> targetPoints;
  std::vector<OpenMesh::VertexHandle> targetVHandles;
  std::vector<Vector3d> targetNormals;
  std::vector<bool> targetBorders;
  std::vector<std::pair<int, int>> stPairs;
  Mesh *m_mesh;
  OpenMesh::VPropHandleT<OpenMesh::Vec4d> m_weights;
  OpenMesh::VPropHandleT<OpenMesh::Vec4i> m_wXids;

 private:
  // void flag_components(OpenMesh::VPropHandleT<int> comp_id, Mesh *mesh,
  //                      int &numOfComp) const;
};