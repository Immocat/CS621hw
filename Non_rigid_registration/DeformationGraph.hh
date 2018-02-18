#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <vector>
#include <unordered_set>
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
                   const std::vector<unsigned int> &src_indices);
  void init(const std::vector<Vector3d> &xi_pos, Mesh *mesh,
            const std::vector<OpenMesh::VertexHandle> &vHandles,
            const std::vector<int> &sample_ids,
            const std::vector<unsigned int> &src_indices);
  // public data
  std::vector<Vector3d> X;
  std::vector<std::unordered_set<int>> X_edges;
  std::vector<std::vector<int>> Xin;
  Mesh *m_mesh;
  OpenMesh::VPropHandleT<OpenMesh::Vec4d> m_weights;
  OpenMesh::VPropHandleT<OpenMesh::Vec4i> m_wXids;

 private:
  void flag_components(OpenMesh::VPropHandleT<int> comp_id, Mesh *mesh,
                       int &numOfComp) const;
};