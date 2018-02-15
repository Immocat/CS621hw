#include "DeformationGraph.hh"
DeformationGraph::DeformationGraph(
    const std::vector<Vector3d> &xi_pos, const Mesh &mesh,
    const std::vector<OpenMesh::VertexHandle> &vHandles,
    const std::vector<int> &sample_ids) {
  init(xi_pos, mesh, vHandles, sample_ids);
}
void DeformationGraph::init(const std::vector<Vector3d> &xi_pos,
                            const Mesh &mesh,
                            const std::vector<OpenMesh::VertexHandle> &vHandles,
                            const std::vector<int> &sample_ids) {
  // 1. get all graph nodes
  X.clear();
  for(int i = 0;i < sample_ids.size(); ++i){
    X.push_back(xi_pos[sample_ids[i]]);
  }
}