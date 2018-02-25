#include <ceres/ceres.h>
#include <Eigen/Core>
#include <algorithm>
#include <limits>
#include <utility>
#include "ClosestPoint.hh"
#include "CostFunctions.hh"
#include "DeformationGraph.hh"
#include "ParallelFor.hh"
#include "Point_Hash_Grid_Searcher3.hh"
#include "RegistrationViewer.hh"
#include "Timer.hh"
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::ScaledLoss;
using ceres::Solve;
using ceres::Solver;
// static bool sample_valid(const Vector3d &p, const std::vector<int> &ids,
//                          const std::vector<Vector3d> &_pts, float len2) {
//   for (int j = 0; j < ids.size(); ++j) {
//     if (length2(p - _pts[ids[j]]) < len2) return false;
//   }
//   return true;
// }
void RegistrationViewer::subsample(const std::vector<Vector3d> &pts,
                                   std::vector<int> &sample_ids,
                                   const double avgDis) {
  //
  if (pts.size() == 0) return;
  double r_threshold = 3 * avgDis;
  // double r2 = r_threshold * r_threshold;

  PointHashGridSearcher3 searcher;
  for (int i = 0; i < pts.size(); ++i) {
    // if (sample_valid(pts[i], sample_ids, pts, r2)) {
    //   sample_ids.push_back(i);
    // }
    if (!searcher.hasNearbyPoint(pts[i], r_threshold)) {
      sample_ids.push_back(i);
      searcher.add(pts[i]);
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
static void addRigidE(std::vector<Transformation> &X_T, double al_reg,
                      ceres::Problem &coarse_align_problem) {
  // for each X[i]'s
  for (int xid = 0; xid < X_T.size(); ++xid) {
    // 1 - ai.dot(ai)
    for (int i = 0; i < 3; ++i) {
      CostFunction *cost_function =
          new AutoDiffCostFunction<OneMinusSelfDotResidual, 1, 3>(
              new OneMinusSelfDotResidual());
      coarse_align_problem.AddResidualBlock(
          cost_function, new ScaledLoss(NULL, al_reg, ceres::TAKE_OWNERSHIP),
          X_T[xid].rotation_[i]);
    }
    CostFunction *cost_function01 =
        new AutoDiffCostFunction<Dot3dResidual, 1, 3, 3>(new Dot3dResidual());
    CostFunction *cost_function02 =
        new AutoDiffCostFunction<Dot3dResidual, 1, 3, 3>(new Dot3dResidual());
    CostFunction *cost_function12 =
        new AutoDiffCostFunction<Dot3dResidual, 1, 3, 3>(new Dot3dResidual());
    coarse_align_problem.AddResidualBlock(
        cost_function01, new ScaledLoss(NULL, al_reg, ceres::TAKE_OWNERSHIP),
        X_T[xid].rotation_[0], X_T[xid].rotation_[1]);
    coarse_align_problem.AddResidualBlock(
        cost_function02, new ScaledLoss(NULL, al_reg, ceres::TAKE_OWNERSHIP),
        X_T[xid].rotation_[0], X_T[xid].rotation_[2]);
    coarse_align_problem.AddResidualBlock(
        cost_function12, new ScaledLoss(NULL, al_reg, ceres::TAKE_OWNERSHIP),
        X_T[xid].rotation_[1], X_T[xid].rotation_[2]);
  }
}
static void addSmoothE(const std::vector<Vector3d> &X_Pos,
                       const std::vector<std::unordered_set<int>> &X_edges,
                       double al_smooth, std::vector<Transformation> &X_T,
                       ceres::Problem &coarse_align_problem) {
  // for each X[i]
  for (int xid = 0; xid < X_Pos.size(); ++xid) {
    const std::unordered_set<int> &X_ie(X_edges[xid]);
    for (const int &j : X_ie) {
      //
      CostFunction *cost_function =
          new AutoDiffCostFunction<SmoothResidual, 1, 3, 3, 3, 3, 3>(
              new SmoothResidual(X_Pos[xid], X_Pos[j]));
      coarse_align_problem.AddResidualBlock(
          cost_function, new ScaledLoss(NULL, al_smooth, ceres::TAKE_OWNERSHIP),
          X_T[xid].rotation_[0], X_T[xid].rotation_[1], X_T[xid].rotation_[2],
          &(X_T[xid].translation_[0]), &(X_T[j].translation_[0]));
    }
  }
}
static void addFitE(const std::vector<Vector3d> &X_pos,
                    std::vector<Transformation> &X_T,
                    const std::vector<Vector3d> &all_targetPoints,
                    const std::vector<Vector3d> &all_targetNormals,
                    const std::vector<std::pair<int, int>> &stPairs,
                    double al_fit, ceres::Problem &coarse_align_problem) {
  // TODO
  double al_point = al_fit * 0.1;
  double al_plane = al_fit;
  for (const auto &stPair : stPairs) {
    int xid = stPair.first;
    int tid = stPair.second;
    const Vector3d &xi(X_pos[xid]);
    const Vector3d &ci(all_targetPoints[tid]);
    const Vector3d &ni(all_targetNormals[tid]);
    CostFunction *cost_function_point =
        new AutoDiffCostFunction<PointResidual, 1, 3, 3, 3, 3>(
            new PointResidual(xi, ci));
    CostFunction *cost_function_plane =
        new AutoDiffCostFunction<PlaneResidual, 1, 3, 3, 3, 3>(
            new PlaneResidual(xi, ci, ni));
    // add point energy
    coarse_align_problem.AddResidualBlock(
        cost_function_point,
        new ScaledLoss(NULL, al_point, ceres::TAKE_OWNERSHIP),
        X_T[xid].rotation_[0], X_T[xid].rotation_[1], X_T[xid].rotation_[2],
        &(X_T[xid].translation_[0]));
    coarse_align_problem.AddResidualBlock(
        cost_function_plane,
        new ScaledLoss(NULL, al_plane, ceres::TAKE_OWNERSHIP),
        X_T[xid].rotation_[0], X_T[xid].rotation_[1], X_T[xid].rotation_[2],
        &(X_T[xid].translation_[0]));
  }
}
static void addFineFitE(OpenMesh::TriMesh_ArrayKernelT<> &mesh,
                        const OpenMesh::VPropHandleT<bool> &hasTarget,
                        const OpenMesh::VPropHandleT<Vector3d> &targetPoints,
                        OpenMesh::VPropHandleT<Transformation> &M_trans,
                        ceres::Problem &fine_align_problem) {
  //
  typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;
  Mesh::ConstVertexIter v_it(mesh.vertices_begin()), v_end(mesh.vertices_end());
  for (; v_it != v_end; ++v_it) {
    // pruned
    if (!mesh.property(hasTarget, *v_it)) continue;

    const OpenMesh::Vec3f &p = mesh.point(*v_it);
    Transformation &trans = mesh.property(M_trans, *v_it);
    const Vector3d &ci = mesh.property(targetPoints, *v_it);
    CostFunction *cost_function_point =
        new AutoDiffCostFunction<PointResidualFine, 1, 9, 3>(
            new PointResidualFine(p, ci));
    fine_align_problem.AddResidualBlock(cost_function_point, NULL,
                                        trans.rotation_[0],
                                        &(trans.translation_[0]));
  }
}
static void addFineReg(OpenMesh::TriMesh_ArrayKernelT<> &mesh,
                       const OpenMesh::VPropHandleT<bool> &hasTarget,
                       const OpenMesh::VPropHandleT<Vector3d> &targetPoints,
                       OpenMesh::VPropHandleT<Transformation> &M_trans,
                       ceres::Problem &fine_align_problem) {
  //
  typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;
  Mesh::ConstEdgeIter e_it(mesh.edges_begin()), e_end(mesh.edges_end());
  for (; e_it != e_end; ++e_it) {
    HalfedgeHandle he = mesh.halfedge_handle(*e_it, 0);
    OpenMesh::VertexHandle v0 = mesh.to_vertex_handle(he);
    OpenMesh::VertexHandle v1 = mesh.from_vertex_handle(he);
    bool hit0 = mesh.property(hasTarget, v0);
    bool hit1 = mesh.property(hasTarget, v1);
    // cannot add rig energy if any one of vertex's correspond is pruned
    if ((!hit0) && (!hit1)) continue;
    //
    const OpenMesh::Vec3f &p0 = mesh.point(v0);
    const OpenMesh::Vec3f &p1 = mesh.point(v1);
    Transformation &trans0 = mesh.property(M_trans, v0);
    Transformation &trans1 = mesh.property(M_trans, v1);
    const Vector3d &c0 = mesh.property(targetPoints, v0);
    const Vector3d &c1 = mesh.property(targetPoints, v1);
    CostFunction *cost_function_point =
        new AutoDiffCostFunction<RigidResidualFine, 1, 9, 3, 9, 3>(
            new RigidResidualFine(p0, p1, c0, c1));
    fine_align_problem.AddResidualBlock(
        cost_function_point, NULL, trans0.rotation_[0],
        &(trans0.translation_[0]), trans1.rotation_[0],
        &(trans1.translation_[0]));
  }
}
inline static bool converge_coarseE(double Ek, double Ekm1) {
  static const double sigma = 0.005;
  return (std::abs(Ek - Ekm1) < sigma * Ek) ||
                 (Ek < std::numeric_limits<double>::epsilon())
             ? true
             : false;
}
static void calculate_correspondences(const ClosestPoint &targetCP,
                                      DeformationGraph &DG) {
  // 1. calculate new pos and normal for each DG node
  std::vector<Vector3d> X_newPos;
  std::vector<Vector3d> X_newNormal;
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
  static const double normalCompatabilityThresh = 90;
  static const double cosineThresh =
      std::cos(normalCompatabilityThresh * M_PI / 180.0);
  // distance threshold is 3 times the median distance
  static const double distMedianThresh = 10;
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
                    DG.targetNormals[tarId]) <= cosineThresh)
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
  glutPostRedisplay();
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
  size_t n_lastPairs = 0;
  static const double al_fit = 0.1;
  for (int n_out = 0; n_out < 100 && al_reg >= 0.1; ++n_out) {
    double E_km1 = 0;
    double E_k = 0;
    for (int n_in = 0; n_in < 20; ++n_in) {
      // ICP & prune

      // std::vector<Vector3d> X_newPos;
      // std::vector<Vector3d> X_newNormal;
      calculate_correspondences(targetCP, M_DG);
      printf("\t[ICP]: find %d pairs\n", (int)M_DG.stPairs.size());

      // global optimizatio
      ceres::Problem coarse_align_problem;
      // given deformationgraph, scale: al_reg, add Erigid cost function into
      // problem.
      // std::vector<Transformation> newTransforms(M_DG.stPairs.size());
      addRigidE(M_DG.X_T, al_reg, coarse_align_problem);
      addSmoothE(M_DG.X, M_DG.X_edges, al_reg * 0.1, M_DG.X_T,
                 coarse_align_problem);
      addFitE(M_DG.X, M_DG.X_T, M_DG.targetPoints, M_DG.targetNormals,
              M_DG.stPairs, al_fit, coarse_align_problem);

      // Solve it  based on Cholesky decomposition
      Solver::Options options;
      options.max_num_iterations = 25;
      //options.max_solver_time_in_seconds = 10.0;
      options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
      options.minimizer_progress_to_stdout = false;
      Solver::Summary summary;
      /////////////////////////////!!!!!!!!!!!!!!!!!!!!!!!!!!
      {
        // when optimization, do not draw
        // std::lock_guard<std::mutex> guard(M_DG.mutex);
        Solve(options, &coarse_align_problem, &summary);
        glutPostRedisplay();
      }
      //////////////////////////////!!!!!!!!!!!!!!!!!!!!!!!
      // update E_k
      E_k = summary.final_cost;
      printf(
          "\t[Global Optimize]: out %d iters, inner %d iters, before E: %f, "
          "after E: %f, solver_iter %d steps, takes %f second\n",
          n_out, n_in, summary.initial_cost, E_k, summary.num_successful_steps,
          summary.total_time_in_seconds);
      // try if converge, and not quit at first inner iter
      if (converge_coarseE(E_k, E_km1) && n_in > 0 &&
          M_DG.stPairs.size() == n_lastPairs)
        break;
      // not converge update E_k, E_km1
      n_lastPairs = M_DG.stPairs.size();
      E_km1 = E_k;
      // X_T[i] = newMatrix[i]  * X_T[i] for all pairs
      // M_DG.updateTransforms(newTransforms);
    }
    // Great! we can relax al_reg by 10
    al_reg *= 0.1;
  }
  // 6. transform each vertex to new pos

  Timer transform_vn_timer;
  M_DG.transformVandN();
  printf("\t[transformVandN]: took %s\n",
         transform_vn_timer.elapsedString().c_str());

  // 7. clean
  src_mesh->remove_property(comp_id);
}

// fineLiearAlignment M to S
void RegistrationViewer::fineLinearAlignment(
    Mesh &src_mesh, const Mesh &target_mesh,
    const std::vector<unsigned int> &target_indices) {
  // given a source mesh, transform it to tartget_mesh's pos (finely)
  // no topology changes happen here.
  printf("[Fine Linear Alignment]: Start\n");
  build_target_bvh();
  // get source mesh
  std::vector<Vector3d> srcPoints;
  std::vector<OpenMesh::VertexHandle> srcVHandles;
  std::vector<Vector3d> srcNormals;
  // std::vector<Vector3d> srcNormals;
  get_points(src_mesh, srcPoints, srcVHandles);
  get_normals(src_mesh, srcNormals);
  // std::vector<Vector3d> targetPoints(srcPoints.size());

  src_mesh.add_property(targetPoints);
  src_mesh.add_property(M_fineAlignTrans);
  src_mesh.add_property(hasTarget);
  src_mesh.add_property(displacement);
  for (int i = 0; i < srcPoints.size(); ++i) {
    // init fine alignment trans
    src_mesh.property(M_fineAlignTrans, srcVHandles[i]) = Transformation();
    src_mesh.property(hasTarget, srcVHandles[i]) = false;
    src_mesh.property(displacement, srcVHandles[i]) = 0;
    src_mesh.property(targetPoints, srcVHandles[i]) = srcPoints[i];

  }

  ///
  ClosestPoint targetCP;
  std::vector<Vector3d> targetP;
  std::vector<OpenMesh::VertexHandle> targetVHandles;
  get_points(target_mesh, targetP, targetVHandles);
  targetCP.init(targetP);
  // draw_finealign_mutex.lock();
  draw_fineAlign_intermediate = true;

  ////////////////////////////////////////////
  // TODO: real alignment happens here
  // for each vertex on src_mesh, trace a undirected ray to get closet point on
  // target mesh, also init M_fineAlignTrans to I and (0,0,0)
  // nanort::BVHTraceOptions trace_options;  // default
  // for (int n_solve = 0; n_solve < 1; ++n_solve) {
  nanort::TriangleIntersector<float, nanort::TriangleIntersection<float>>
      triangle_intersecter((const float *)target_mesh.points(),
                           target_indices.data(), sizeof(float) * 3);
  Timer trace_timer;
  const float tFar = 1e30;  // m_bbDiagnol * 0.01;
  std::vector<float> ts;
  ts.reserve(srcPoints.size());
  std::mutex ts_mutex;
  // ParallelFor(0, (int)srcPoints.size(), [&](int i) {
  //   ////
  //   //const Vector3d& closetTargetP =
  //   targetP[targetCP.getClosestPoint(srcPoints[i])];
  //   //if(length2(closetTargetP - srcPoints[i]) < 0.01);
  //   ////
  //   nanort::Ray<float> ray0;
  //   srcNormals[i] = srcNormals[i].normalize();
  //   ray0.min_t = 0.001f;
  //   ray0.max_t = tFar;
  //   ray0.org[0] = srcPoints[i][0];
  //   ray0.org[1] = srcPoints[i][1];
  //   ray0.org[2] = srcPoints[i][2];

  //   ray0.dir[0] = srcNormals[i][0];
  //   ray0.dir[1] = srcNormals[i][1];
  //   ray0.dir[2] = srcNormals[i][2];
  //   nanort::TriangleIntersection<float> isect0;
  //   bool hit0 = S_BVH.Traverse(ray0, triangle_intersecter, &isect0);

  //   nanort::Ray<float> ray1;

  //   ray1.min_t = 0.001f;
  //   ray1.max_t = tFar;
  //   ray1.org[0] = srcPoints[i][0];
  //   ray1.org[1] = srcPoints[i][1];
  //   ray1.org[2] = srcPoints[i][2];
  //   ray1.dir[0] = -srcNormals[i][0];
  //   ray1.dir[1] = -srcNormals[i][1];
  //   ray1.dir[2] = -srcNormals[i][2];
  //   nanort::TriangleIntersection<float> isect1;
  //   bool hit1 = S_BVH.Traverse(ray1, triangle_intersecter, &isect1);
  //   if (hit0 || hit1) {
  //     // found an intersection, set target point for this vertex
  //     OpenMesh::VertexHandle vh(srcVHandles[i]);

  //     if ((hit0 && !hit1) || (hit0 && hit1 && isect0.t < isect1.t)) {
  //       Vector3d tar(srcPoints[i] + (double)isect0.t * srcNormals[i]);
  //       if (length2(tar - srcPoints[i]) < 0.01) {
  //         src_mesh.property(targetPoints, vh) = tar;
  //         src_mesh.property(displacement, vh) = isect0.t;
  //         src_mesh.property(hasTarget, vh) = true;

  //         {
  //           std::lock_guard<std::mutex> guard(ts_mutex);
  //           ts.emplace_back(isect0.t);
  //         }
  //       }
  //     } else if ((hit1 && !hit0) || (hit0 && hit1 && isect0.t >= isect1.t)) {
  //       Vector3d tar(srcPoints[i] + (double)(-isect1.t) * (srcNormals[i]));
  //       if (length2(tar - srcPoints[i]) < 0.01) {
  //         src_mesh.property(targetPoints, vh) = tar;
  //         src_mesh.property(displacement, vh) = isect1.t;
  //         src_mesh.property(hasTarget, vh) = true;

  //         {
  //           std::lock_guard<std::mutex> guard(ts_mutex);
  //           ts.emplace_back(isect1.t);
  //         }
  //       }
  //     }
  //   }
  // });
  for (int i = 0; i < srcPoints.size(); ++i) {
    ////
    const Vector3d &closetTargetP =
        targetP[targetCP.getClosestPoint(srcPoints[i])];
    if (length2(closetTargetP - srcPoints[i]) < 0.0001) continue;
    ////
    nanort::Ray<float> ray0;
    srcNormals[i] = srcNormals[i].normalize();
    ray0.min_t = 0.0f;
    ray0.max_t = tFar;
    ray0.org[0] = srcPoints[i][0];
    ray0.org[1] = srcPoints[i][1];
    ray0.org[2] = srcPoints[i][2];

    ray0.dir[0] = srcNormals[i][0];
    ray0.dir[1] = srcNormals[i][1];
    ray0.dir[2] = srcNormals[i][2];
    nanort::TriangleIntersection<float> isect0;
    bool hit0 = S_BVH.Traverse(ray0, triangle_intersecter, &isect0);

    nanort::Ray<float> ray1;

    ray1.min_t = 0.0f;
    ray1.max_t = tFar;
    ray1.org[0] = srcPoints[i][0];
    ray1.org[1] = srcPoints[i][1];
    ray1.org[2] = srcPoints[i][2];
    ray1.dir[0] = -srcNormals[i][0];
    ray1.dir[1] = -srcNormals[i][1];
    ray1.dir[2] = -srcNormals[i][2];
    nanort::TriangleIntersection<float> isect1;
    bool hit1 = S_BVH.Traverse(ray1, triangle_intersecter, &isect1);
    if (hit0 || hit1) {
      // found an intersection, set target point for this vertex
      OpenMesh::VertexHandle vh(srcVHandles[i]);

      if ((hit0 && !hit1) || (hit0 && hit1 && isect0.t < isect1.t)) {
        Vector3d tar(srcPoints[i] + (double)isect0.t * srcNormals[i]);
        if (length2(tar - srcPoints[i]) < 0.01) {
          src_mesh.property(targetPoints, vh) = tar;
          src_mesh.property(displacement, vh) = isect0.t;
          src_mesh.property(hasTarget, vh) = true;

          {
            std::lock_guard<std::mutex> guard(ts_mutex);
            ts.emplace_back(isect0.t);
          }
        }
      } else if ((hit1 && !hit0) || (hit0 && hit1 && isect0.t >= isect1.t)) {
        Vector3d tar(srcPoints[i] + (double)(-isect1.t) * (srcNormals[i]));
        if (length2(tar - srcPoints[i]) < 0.01) {
          src_mesh.property(targetPoints, vh) = tar;
          src_mesh.property(displacement, vh) = isect1.t;
          src_mesh.property(hasTarget, vh) = true;

          {
            std::lock_guard<std::mutex> guard(ts_mutex);
            ts.emplace_back(isect1.t);
          }
        }
      }
    }
  }
  printf("\t[trace undirected ray]: %s\n", trace_timer.elapsedString().c_str());

  // prune pairs with long correspond
  // std::nth_element(ts.begin(), ts.begin() + ts.size() / 2, ts.end());
  // static const double distMedianThresh = 10;

  double distThreshold =
      0.01 * m_bbDiagnol;  // distMedianThresh * ts[ts.size() / 2];
  for (int i = 0; i < srcVHandles.size(); ++i) {
    if (src_mesh.property(hasTarget, srcVHandles[i])) {
      //
      if (src_mesh.property(displacement, srcVHandles[i]) > distThreshold) {
        src_mesh.property(hasTarget, srcVHandles[i]) = false;
      }
    }
  }

  // // global optimization
  // // global optimizatio
  ceres::Problem fine_align_problem;
  // given deformationgraph, scale: al_reg, add Erigid cost function into
  // problem.
  // std::vector<Transformation> newTransforms(M_DG.stPairs.size());
  addFineFitE(src_mesh, hasTarget, targetPoints, M_fineAlignTrans,
              fine_align_problem);
  addFineReg(src_mesh, hasTarget, targetPoints, M_fineAlignTrans,
             fine_align_problem);

  // Solve it  based on Cholesky decomposition
  Solver::Options options;
  options.max_num_iterations = 100;
  // options.max_solver_time_in_seconds = 10.0;
  options.linear_solver_type = ceres::CGNR;
  options.minimizer_progress_to_stdout = true;
  Solver::Summary summary;
  /////////////////////////////!!!!!!!!!!!!!!!!!!!!!!!!!!
  {
    // when optimization, do not draw
    // std::lock_guard<std::mutex> guard(M_DG.mutex);
    Solve(options, &fine_align_problem, &summary);
    glutPostRedisplay();
  }

  //}
  // sleep(300);
  ////////////////////////////////////////////
  // finally transform mesh to new place
  // for (int i = 0; i < srcVHandles.size(); ++i) {
  //   OpenMesh::VertexHandle vh(srcVHandles[i]);
  //   const Transformation &trans(src_mesh.property(M_fineAlignTrans, vh));
  //   const OpenMesh::Vec3f &p(src_mesh.point(vh));
  //   Vector3d pos(p[0], p[1], p[2]);
  //   pos = trans.transformPoint(pos);
  //   src_mesh.point(vh) = OpenMesh::Vec3f(pos[0], pos[1], pos[2]);
  // }
  // src_mesh.update_normals();

  //
  // draw_fineAlign_intermediate = false;
  // draw_finealign_mutex.unlock();

  // clean up
  // src_mesh.remove_property(hasTarget);
  // src_mesh.remove_property(targetPoints);
  // src_mesh.remove_property(M_fineAlignTrans);
  // src_mesh.remove_property(displacement);
}