#include <ceres/ceres.h>
#include <algorithm>
#include <utility>
#include "ClosestPoint.hh"
#include "DeformationGraph.hh"
#include "ParallelFor.hh"
#include "RegistrationViewer.hh"
#include "Timer.hh"
static bool sample_valid(const Vector3d &p, const std::vector<int> &ids,
                         const std::vector<Vector3d> &_pts, float len2) {
  for (int j = 0; j < ids.size(); ++j) {
    if (length2(p - _pts[ids[j]]) < len2) return false;
  }
  return true;
}
void RegistrationViewer::subsample(const std::vector<Vector3d> &pts,
                                   std::vector<int> &sample_ids,
                                   const double avgDis) {
  //
  if (pts.size() == 0) return;
  double r_threshold = 4 * avgDis;
  double r2 = r_threshold * r_threshold;
  for (int i = 0; i < pts.size(); ++i) {
    if (sample_valid(pts[i], sample_ids, pts, r2)) {
      sample_ids.push_back(i);
    }
  }
}
void RegistrationViewer::flag_components(OpenMesh::VPropHandleT<int> comp_id,
                                         Mesh *mesh, int &numOfComps) {
  numOfComps = 0;
  // iterate over all vertices
  Mesh::VertexIter v_it = mesh->vertices_begin();
  Mesh::VertexIter v_end = mesh->vertices_end();
  for (Mesh::VertexIter v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
    mesh->property(comp_id, *v_it) =
        -1;  // init, all vertex belongs to -1(invalid) component

  Mesh::VertexHandle vh;
  Mesh::VertexIter current_pos = mesh->vertices_begin();

  while (true) {
    // find an unvisited vertex
    bool found = false;
    for (v_it = current_pos; v_it != v_end; ++v_it)
      if (mesh->property(comp_id, *v_it) == -1) {
        found = true;
        vh = *v_it;
        mesh->property(comp_id, *v_it) = numOfComps;
        current_pos = v_it;
        break;
      }

    // if none was found -> finished
    if (!found) break;

    numOfComps++;

    std::vector<Mesh::VertexHandle> handles;
    handles.push_back(vh);

    // grow from found vertex
    while (handles.size() > 0) {
      Mesh::VertexHandle current = handles.back();
      handles.pop_back();

      Mesh::VertexVertexIter vv_it;

      for (vv_it = mesh->vv_iter(current); vv_it.is_valid(); ++vv_it)
        if (mesh->property(comp_id, *vv_it) == -1) {
          mesh->property(comp_id, *vv_it) = numOfComps - 1;
          handles.push_back(*vv_it);
        }
    }
  }
}
// given deformationgraph, scale: al_reg, add Erigid cost function into problem.
static void addRigidE(std::vector<Transformation> &newTransforms, double al_reg,
                      ceres::Problem &coarse_align_problem) {
  // TODO
}
static void addSmoothE(const std::vector<Vector3d> &X_newPos,
                       const std::vector<std::unordered_set<int>> &X_edges,
                       double al_smooth,
                       std::vector<Transformation> &newTransforms,
                       ceres::Problem &coarse_align_problem) {
  // TODO
}
static void addFitE(const std::vector<Vector3d> &X_newPos,
                    const std::vector<Vector3d> &all_targetPoints,
                    const std::vector<Vector3d> &all_targetNormals,
                    const std::vector<std::pair<int, int>> &stPairs,
                    double al_fit, ceres::Problem &coarse_align_problem) {
  // TODO
}
inline static bool converge_coarseE(double Ek, double Ekm1) {
  static const double sigma = 0.005;
  return (std::abs(Ek - Ekm1) < sigma * Ek) ? true : false;
}
static void calculate_correspondences(
    const ClosestPoint &targetCP, DeformationGraph &DG,
    std::vector<Vector3d> &X_newPos, std::vector<Vector3d> &X_newNormal) {
  // 1. calculate new pos and normal for each DG node
  X_newPos.clear();
  X_newNormal.clear();
  DG.stPairs.clear();
  std::vector<std::pair<int, int>> candidate_pairs;
  std::vector<double> src_target_dis2;

  // at most each DG node has a candidate
  X_newPos.reserve(DG.X.size());
  X_newNormal.reserve(DG.X.size());
  DG.stPairs.reserve(DG.X.size());
  candidate_pairs.reserve(DG.X.size());
  for (int i = 0; i < DG.X.size(); ++i) {
    X_newPos.emplace_back(DG.X_T[i].transformPoint(DG.X[i]));
    X_newNormal.emplace_back(DG.X_T[i].transformVector(DG.X_normal[i]));
  }
  // 2. for each new pos, find nearest point from all_target_pts, get all
  // candidates pairs
  for (int i = 0; i < X_newPos.size(); ++i) {
    int bestIndex = targetCP.getClosestPoint(X_newPos[i]);
    // do not keep border correspondences
    if (!DG.targetBorders[bestIndex]) {
      candidate_pairs.emplace_back(i, bestIndex);
      src_target_dis2.emplace_back(
          length2(X_newPos[i] - DG.targetPoints[bestIndex]));
    }
  }
  // 3. prune
  // normals of correspondences do not deviate more than 60 degrees
  static const double normalCompatabilityThresh = 60;
  static const double cosineThresh =
      std::cos(normalCompatabilityThresh * M_PI / 180.0);
  // distance threshold is 3 times the median distance
  static const double distMedianThresh = 3;
  std::vector<double> dis2_copy = src_target_dis2;
  std::nth_element(src_target_dis2.begin(),
                   src_target_dis2.begin() + src_target_dis2.size() / 2,
                   src_target_dis2.end());
  double distThreshold =
      distMedianThresh *
      std::sqrt(
          std::max((double)0, src_target_dis2[src_target_dis2.size() / 2]));
  double distThreshold2 = distThreshold * distThreshold;
  // double cosineThresh = std::cos(normalCompatabilityThresh * M_PI / 180.0);
  for (const std::pair<int, int> &stPair : candidate_pairs) {
    int newPosid = stPair.first;
    int tarId = stPair.second;
    if (dis2_copy[newPosid] > distThreshold2) continue;
    // suppose distMedian is always less than 90 degree, which means that
    // cosineThresh<=cos(theta)<=1
    //!!!!!!!!!!!!!!!!!! Here I suppose that target(openmesh calculates) normals
    //! are all normalized!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (dot_product(X_newNormal[newPosid].normalize(),
                    DG.targetNormals[tarId]) < cosineThresh)
      continue;
    // put into output
    DG.stPairs.emplace_back(stPair);
  }
}
// Coarse alignment M to S
void RegistrationViewer::coarseNonRigidAlignment(
    Mesh *src_mesh, const std::vector<unsigned int> &src_indices,
    const Mesh &target_mesh) {
  // given a source mesh, transform it to tartget_mesh's pos (coarsely)
  // no topology changes happen here.
  // 1. Find all connected components
  printf("[CoarseAlignment]: Start\n");
  int numOfComps = 0;
  OpenMesh::VPropHandleT<int> comp_id;
  src_mesh->add_property(comp_id);
  flag_components(comp_id, src_mesh, numOfComps);
  printf("\t[DeformationGraph]: found %d connected component(s)\n", numOfComps);
  // 2. Uniform subsample source points in each connected component

  std::vector<Vector3d> srcPoints;
  std::vector<int> srcSampleIds;
  std::vector<OpenMesh::VertexHandle> srcVHandles;
  // std::vector<Vector3d> srcNormals;
  get_points(*src_mesh, srcPoints, srcVHandles);
  // get_normals(*src_mesh, srcNormals);
  //
  std::vector<std::vector<Vector3d>> srcPoints_cc(numOfComps);
  std::vector<std::vector<int>> srcSampleIds_cc(numOfComps);
  std::vector<std::vector<int>> srcSampleIdscc2srcPointsid(numOfComps);
  for (int i = 0; i < srcPoints.size(); ++i) {
    //
    OpenMesh::VertexHandle vh = srcVHandles[i];
    int cid = src_mesh->property(comp_id, vh);
    srcPoints_cc[cid].emplace_back(srcPoints[i]);
    srcSampleIdscc2srcPointsid[cid].emplace_back(i);
  }
  // get averagevertexdistance as subsample threshold

  averageVertexDistance_ = get_average_vertex_distance(*src_mesh);
  Timer subsample_timer;
  // ParallelFor(0,numOfComps,[&] (int i) {
  // 	subsample(srcPoints_cc[i], srcSampleIds_cc[i], averageVertexDistance_);
  // });
  for (int i = 0; i < numOfComps; ++i) {
    subsample(srcPoints_cc[i], srcSampleIds_cc[i], averageVertexDistance_);
  }
  std::string subsample_elapseString = subsample_timer.elapsedString();
  // gather all samples into srcSampleIds
  for (int i = 0; i < srcSampleIds_cc.size(); ++i) {
    for (int j = 0; j < srcSampleIds_cc[i].size(); ++j) {
      int id = srcSampleIdscc2srcPointsid[i][srcSampleIds_cc[i][j]];
      srcSampleIds.push_back(id);
    }
  }

  printf("\t[Subsample]: sampled %d vertices, took %s\n",
         (int)srcSampleIds.size(), subsample_elapseString.c_str());
  // 3. Build deformation graph
  M_DG.init(srcPoints, src_mesh, srcVHandles, srcSampleIds, src_indices,
            numOfComps, comp_id);
  // 4. build target mesh kd-tree
  // clear and set target's pos, normal, border infos
  get_points(target_mesh, M_DG.targetPoints, M_DG.targetVHandles);
  get_normals(target_mesh, M_DG.targetNormals);
  get_borders(target_mesh, M_DG.targetBorders);
  ClosestPoint targetCP;
  targetCP.init(M_DG.targetPoints);

  // 5. main iteration, core of this project!
  // !!!!!!!!!!!!if done too musch times outer, try use al_reg > 0.1
  // or al_reg *= 0.5
  double al_reg = 1000.0;  // init as 1000, relax it if converge
  static const double al_fit = 0.1;
  for (int n_out = 0; n_out < 100 && al_reg >= 0.1; ++n_out) {
    double E_km1 = 0;
    double E_k = 0;
    for (int n_in = 0;; ++n_in) {
      // ICP & prune
      
      std::vector<Vector3d> X_newPos;
      std::vector<Vector3d> X_newNormal;
      calculate_correspondences(targetCP, M_DG, X_newPos, X_newNormal);

      // global optimization
      ceres::Problem coarse_align_problem;
      // given deformationgraph, scale: al_reg, add Erigid cost function into
      // problem.
      std::vector<Transformation> newTransforms(M_DG.stPairs.size());
      addRigidE(newTransforms, al_reg, coarse_align_problem);
      addSmoothE(X_newPos, M_DG.X_edges, al_reg * 0.1, newTransforms,
                 coarse_align_problem);
      addFitE(X_newPos, M_DG.targetPoints, M_DG.targetNormals, M_DG.stPairs, al_fit,
              coarse_align_problem);

      // Solve it  based on Cholesky decomposition
      // TODO: solve it
      // update E_k
      // TODO :E_k = ceres thing
      // try if converge, and not quit at first inner iter
      if (converge_coarseE(E_k, E_km1) && n_in > 0) break;
      // not converge update E_k, E_km1
      E_km1 = E_k;
      // X_T[i] = newMatrix[i]  * X_T[i] for all pairs
      M_DG.updateTransforms(newTransforms);
    }
    // Great! we can relax al_reg by 10
    al_reg *= 0.1;
  }

  // 6. Based on graph nodes' transform and weights, get new vertex transform

  // 7. transform each vertex to new pos
}

// fineLiearAlignment M to S
void RegistrationViewer::fineLiearAlignment(Mesh &src_mesh,
                                            const Mesh &target_mesh) {
  // given a source mesh, transform it to tartget_mesh's pos (finely)
  // no topology changes happen here.
}