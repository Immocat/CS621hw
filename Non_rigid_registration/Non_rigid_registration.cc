#include "ClosestPoint.hh"
#include "DeformationGraph.hh"
#include "RegistrationViewer.hh"
static bool sample_valid(const Vector3d &p, const std::vector<int> &ids,
                         const std::vector<Vector3d> &_pts, float len2) {
  for (int j = 0; j < ids.size(); ++j) {
    if (length2(p - _pts[ids[j]]) < len2) return false;
  }
  return true;
}
void RegistrationViewer::subsample(const std::vector<Vector3d> &pts,
                                   std::vector<int> &sample_ids,
                                   const Mesh &mesh) {
  //
  if (pts.size() == 0) return;
  averageVertexDistance_ = get_average_vertex_distance(mesh);
  double r_threshold = 4 * averageVertexDistance_;
  double r2 = r_threshold * r_threshold;
  for (int i = 1; i < pts.size(); ++i) {
    if (sample_valid(pts[i], sample_ids, pts, r2)) {
      sample_ids.push_back(i);
    }
  }
}
// Coarse alignment M to S
void RegistrationViewer::coarseNonRigidAlignment(
    Mesh *src_mesh, const std::vector<unsigned int> &src_indices,
    const Mesh &target_mesh) {
  // given a source mesh, transform it to tartget_mesh's pos (coarsely)
  // no topology changes happen here.
  // 1. uniform subsample source points and build deformation graph
  printf("[CoarseAlignment]: Start\n");
  std::vector<Vector3d> srcPoints;
  std::vector<int> srcSampleIds;
  std::vector<OpenMesh::VertexHandle> srcVHandles;
  get_points(*src_mesh, srcPoints, srcVHandles);
  //
  subsample(srcPoints, srcSampleIds, *src_mesh);
  printf("\t[Subsample]: sampled %d vertices.\n", (int)srcSampleIds.size());
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  M_DG.init(srcPoints, src_mesh, srcVHandles, srcSampleIds, src_indices); 

  // we only operate on deformation graph's nodes

  // based on graph nodes' transform and weights, get new vertex transform
}

// fineLiearAlignment M to S
void RegistrationViewer::fineLiearAlignment(Mesh &src_mesh,
                                            const Mesh &target_mesh) {
  // given a source mesh, transform it to tartget_mesh's pos (finely)
  // no topology changes happen here.
}