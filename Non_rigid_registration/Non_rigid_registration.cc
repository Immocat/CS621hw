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
  double r_threshold = 4 * avgDis;
  //double r2 = r_threshold * r_threshold;

  PointHashGridSearcher3 searcher;
  for (int i = 0; i < pts.size(); ++i) {
    // if (sample_valid(pts[i], sample_ids, pts, r2)) {
    //   sample_ids.push_back(i);
    // }
    if (!searcher.hasNearbyPoint(pts[i],r_threshold)) {
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
  static const double normalCompatabilityThresh = 70;
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
    // if (dis2_copy[newPosid] > distThreshold2) continue;
    // suppose distMedian is always less than 90 degree, which means that
    // cosineThresh<=cos(theta)<=1
    //!!!!!!!!!!!!!!!!!! Here I suppose that target(openmesh calculates) normals
    //! are all normalized!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // if (dot_product(X_newNormal[newPosid].normalize(),
    //                DG.targetNormals[tarId]) <= 0)
    //  continue;
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
  static const double al_fit = 0.1;
  for (int n_out = 0; n_out < 100 && al_reg >= 0.1; ++n_out) {
    double E_km1 = 0;
    double E_k = 0;
    for (int n_in = 0; n_in < 100; ++n_in) {
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
      options.max_num_iterations = 5;
      // options.max_solver_time_in_seconds = 10.0;
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
      if (converge_coarseE(E_k, E_km1) && n_in > 0) break;
      // not converge update E_k, E_km1
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
  // get_normals(src_mesh, srcNormals);
  // std::vector<Vector3d> targetPoints(srcPoints.size());

  OpenMesh::VPropHandleT<Vector3d> targetPoints;
  src_mesh.add_property(targetPoints);
  src_mesh.add_property(M_fineAlignTrans);

  draw_fineAlign_intermediate = true;
  for (int i = 0; i < srcPoints.size(); ++i) {
    // init fine alignment trans
    src_mesh.property(M_fineAlignTrans, srcVHandles[i]) = Transformation();
  }
  ////////////////////////////////////////////
  // TODO: real alignment happens here
  // for each vertex on src_mesh, trace a undirected ray to get closet point on
  // target mesh, also init M_fineAlignTrans to I and (0,0,0)
  // nanort::BVHTraceOptions trace_options;  // default
  nanort::TriangleIntersector<float, nanort::TriangleIntersection<float>>
      triangle_intersecter((const float *)target_mesh.points(),
                           target_indices.data(), sizeof(float) * 3);
  Timer trace_timer;
  static const float tFar = 1.0e+30f;
  ParallelFor(0, (int)srcPoints.size(), [&](int i) {
    nanort::Ray<float> ray;
    ray.min_t = 0.0f;
    ray.max_t = tFar;
    ray.org[0] = srcPoints[i][0];
    ray.org[1] = srcPoints[i][1];
    ray.org[2] = srcPoints[i][2];

    const OpenMesh::Vec3f &cur_normal(src_mesh.normal(srcVHandles[i]));
    ray.dir[0] = cur_normal[0];
    ray.dir[1] = cur_normal[1];
    ray.dir[2] = cur_normal[2];
    nanort::TriangleIntersection<float> isect;
    bool hit = S_BVH.Traverse(ray, triangle_intersecter, &isect);
    ray.dir[0] = -cur_normal[0];
    ray.dir[1] = -cur_normal[1];
    ray.dir[2] = -cur_normal[2];
    hit = S_BVH.Traverse(ray, triangle_intersecter, &isect);
  });
  // for (int i = 0; i < srcPoints.size(); ++i) {
  //   nanort::Ray<float> ray;
  //   ray.min_t = 0.0f;
  //   ray.max_t = tFar;
  //   ray.org[0] = srcPoints[i][0];
  //   ray.org[1] = srcPoints[i][1];
  //   ray.org[2] = srcPoints[i][2];

  //   const OpenMesh::Vec3f &cur_normal(src_mesh.normal(srcVHandles[i]));
  //   ray.dir[0] = cur_normal[0];
  //   ray.dir[1] = cur_normal[1];
  //   ray.dir[2] = cur_normal[2];
  //   nanort::TriangleIntersection<float> isect;
  //   bool hit = S_BVH.Traverse(ray, triangle_intersecter, &isect);
  // }
  printf("\t[trace ray]: %s\n", trace_timer.elapsedString().c_str());
  ////////////////////////////////////////////
  draw_fineAlign_intermediate = false;

  // clean up
  src_mesh.remove_property(targetPoints);
  src_mesh.remove_property(M_fineAlignTrans);
}