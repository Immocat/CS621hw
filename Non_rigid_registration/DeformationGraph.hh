#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <vector>
#include "Vector.hh"
class DeformationGraph {
  typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

 public:
  DeformationGraph() {}
  DeformationGraph(const std::vector<Vector3d> &xi_pos, const Mesh &mesh,
                   const std::vector<OpenMesh::VertexHandle> &vHandles,
                   const std::vector<int> &sample_ids);
  void init(const std::vector<Vector3d> &xi_pos, const Mesh &mesh,
            const std::vector<OpenMesh::VertexHandle> &vHandles,
            const std::vector<int> &sample_ids);
  // public data
  std::vector<Vector3d> X;
  std::vector<std::vector<int>> Xin;
};