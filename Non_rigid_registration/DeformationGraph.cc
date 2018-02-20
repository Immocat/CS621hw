#include "DeformationGraph.hh"
//#include <geodesic/geodesic_algorithm_exact.h>
#include <unordered_map>
#include "ClosestPoint.hh"
#include "Timer.hh"
DeformationGraph::DeformationGraph(
    const std::vector<Vector3d> &xi_pos, Mesh *mesh,
    const std::vector<OpenMesh::VertexHandle> &vHandles,
    const std::vector<int> &sample_ids,
    const std::vector<unsigned int> &src_indices, int numOfComps,
    OpenMesh::VPropHandleT<int> comp_id) {
  init(xi_pos, mesh, vHandles, sample_ids, src_indices, numOfComps, comp_id);
}
// void DeformationGraph::flag_components(OpenMesh::VPropHandleT<int> comp_id,
//                                        Mesh *mesh, int &numOfComps) const {
//   numOfComps = 0;
//   // iterate over all vertices
//   Mesh::VertexIter v_it = mesh->vertices_begin();
//   Mesh::VertexIter v_end = mesh->vertices_end();
//   for (Mesh::VertexIter v_it = mesh->vertices_begin(); v_it != v_end; ++v_it)
//     mesh->property(comp_id, *v_it) =
//         -1;  // init, all vertex belongs to -1(invalid) component

//   Mesh::VertexHandle vh;
//   Mesh::VertexIter current_pos = mesh->vertices_begin();

//   while (true) {
//     // find an unvisited vertex
//     bool found = false;
//     for (v_it = current_pos; v_it != v_end; ++v_it)
//       if (mesh->property(comp_id, *v_it) == -1) {
//         found = true;
//         vh = *v_it;
//         mesh->property(comp_id, *v_it) = numOfComps;
//         current_pos = v_it;
//         break;
//       }

//     // if none was found -> finished
//     if (!found) break;

//     numOfComps++;

//     std::vector<Mesh::VertexHandle> handles;
//     handles.push_back(vh);

//     // grow from found vertex
//     while (handles.size() > 0) {
//       Mesh::VertexHandle current = handles.back();
//       handles.pop_back();

//       Mesh::VertexVertexIter vv_it;

//       for (vv_it = mesh->vv_iter(current); vv_it.is_valid(); ++vv_it)
//         if (mesh->property(comp_id, *vv_it) == -1) {
//           mesh->property(comp_id, *vv_it) = numOfComps - 1;
//           handles.push_back(*vv_it);
//         }
//     }
//   }
// }
void DeformationGraph::updateTransforms(
    const std::vector<Transformation> &newTransforms) {
  // TODO
  for (int i = 0; i < stPairs.size(); ++i) {
    int xid = stPairs[i].first;
    // X_T[xid] = newMatrix[i]  * X_T[xid] , update for all pairs
    X_T[xid] = newTransforms[i] * X_T[xid];
  }
}

void DeformationGraph::init(const std::vector<Vector3d> &v_pos, Mesh *mesh,
                            const std::vector<OpenMesh::VertexHandle> &vHandles,
                            const std::vector<int> &sample_ids,
                            const std::vector<unsigned int> &src_indices,
                            int numOfComps,
                            OpenMesh::VPropHandleT<int> comp_id) {
  // 1. get all graph nodes
  m_mesh = mesh;
  static const int k = 4;
  X.clear();
  X_edges.clear();
  X_normal.clear();
  X_T.clear();
  X.reserve(sample_ids.size());
  X_normal.reserve(sample_ids.size());
  for (int i = 0; i < sample_ids.size(); ++i) {
    X.emplace_back(v_pos[sample_ids[i]]);
    OpenMesh::Vec3f n(mesh->normal(vHandles[sample_ids[i]]));
    X_normal.emplace_back(n[0], n[1], n[2]);
  }
  X_edges.assign(X.size(), std::unordered_set<int>());

  X_T.assign(X.size(), Transformation());  // set all as identity matrix
  // 3. init kd_tree for each connected component
  std::vector<ClosestPoint> closest_meshes(numOfComps);
  std::vector<std::vector<Vector3d>> closest_points(numOfComps);
  std::vector<std::unordered_map<size_t, size_t>> closest_points_id_2_X_id(
      numOfComps);
  for (int i = 0; i < sample_ids.size(); ++i) {
    //
    OpenMesh::VertexHandle vh = vHandles[sample_ids[i]];
    int cid = mesh->property(comp_id, vh);
    closest_points[cid].emplace_back(v_pos[sample_ids[i]]);
    closest_points_id_2_X_id[cid][closest_points[cid].size() - 1] = i;
  }
  // TODO: could parallel, need profile
  //#pragma omp parallel for
  for (int i = 0; i < numOfComps; ++i) {
    closest_meshes[i].init(closest_points[i]);
  }

  // for each vertex, find its k=4 nearest points
  // TODO: could parallel, need profile
  mesh->add_property(m_weights);
  mesh->add_property(m_wXids);
  //#pragma omp parallel for
  Timer v_xtimer;
  for (int i = 0; i < v_pos.size(); ++i) {
    OpenMesh::VertexHandle vh = vHandles[i];
    int cid = mesh->property(comp_id, vh);
    // find k=4 nearest points id of X and weights
    int kids[5] = {-1, -1, -1, -1, -1};
    double dists[5];
    // maybe one connected component's pts size is less than k+1=5, need special
    int closest_n_pts = (int)closest_points[cid].size();
    int kp1_n_pts_min = std::min(closest_n_pts, k + 1);
    closest_meshes[cid].getKClosestPoint(v_pos[i], kp1_n_pts_min, kids, dists);

    //!!!!!!!!!!dists are length2!!!Fuck ANN!!!!!!!!!!!!!!!
    double dist_kp1 = sqrt(dists[kp1_n_pts_min - 1]);

    // change from dists to weight, normalize to 1
    OpenMesh::Vec4d weights;
    OpenMesh::Vec4i w_id(-1, -1, -1, -1);
    double normalizer = 0;
    for (int j = 0; j < kp1_n_pts_min - 1; ++j) {
      w_id[j] = closest_points_id_2_X_id[cid][kids[j]];
      double v_x_dis = length(v_pos[i] - v_pos[sample_ids[w_id[j]]]);
      double w = 1 - v_x_dis / dist_kp1;
      weights[j] = w;
      normalizer += w;
    }
    if (normalizer != 0) {
      normalizer = 1.0f / normalizer;
    }
    weights *= normalizer;
    mesh->property(m_weights, vh) = weights;
    mesh->property(m_wXids, vh) = w_id;

    // add edges in DG
    for (int j = 0; j < kp1_n_pts_min - 1; ++j) {
      int j_xid = w_id[j];
      for (int k = j + 1; k < kp1_n_pts_min - 1; ++k) {
        int k_xid = w_id[k];
        X_edges[j_xid].insert(k_xid);
        X_edges[k_xid].insert(j_xid);
      }
    }
  }
  printf("\t[DeformationGraph]: finish building DG");
  std::cout << ", took " << v_xtimer.elapsedString() << '\n';
  //       std::vector<std::vector<double>> geo_points(numOfComps,
  //                                                   std::vector<double>());
  //   std::vector<std::vector<unsigned int>> geo_indices(
  //       numOfComps, std::vector<unsigned int>());

  //   // interate all faces
  //   Mesh::ConstFaceIter f_it(mesh.faces_sbegin()), f_end(mesh.faces_end());
  //   Mesh::ConstFaceVertexIter fv_it;
  //   std::unordered_map<unsigned int, size_t> handle_idx_2_geoid;

  //   for (; f_it != f_end; ++f_it) {
  //     // find which components it belongs to
  //     fv_it = mesh.cfv_iter(*f_it);
  //     int cid = mesh.property(comp_id, *fv_it);
  //     for (; fv_it.is_valid(); ++fv_it) {
  //       // add face indices to geo_indices[cid] and vertex points to
  //       // geo_points[cid] if needed.
  //       unsigned int geo_ver_id = 0;
  //       if (handle_idx_2_geoid.find((*fv_it).idx()) ==
  //       handle_idx_2_geoid.end()) {
  //         // this vertex never visited, push it in
  //         Mesh::Point mesh_ver(mesh.point(*fv_it));
  //         geo_points[cid].push_back(mesh_ver[0]);
  //         geo_points[cid].push_back(mesh_ver[1]);
  //         geo_points[cid].push_back(mesh_ver[2]);
  //         handle_idx_2_geoid[(*fv_it).idx()] = geo_ver_id =
  //             geo_points[cid].size() / 3 - 1;
  //       } else {
  //         geo_ver_id = handle_idx_2_geoid[(*fv_it).idx()];
  //       }
  //       geo_indices[cid].push_back(geo_ver_id);
  //     }
  //   }
  // // only for
  // //
  // debug///////////////////////////////////////////////////////////////////////
  // #ifndef NDEBUG

  //   int geo_indices_num = 0;
  //   for (const std::vector<unsigned int> &geo_comp_id : geo_indices) {
  //     geo_indices_num += geo_comp_id.size();
  //   }
  //   size_t mesh_face_count = mesh.n_faces();
  //   assert(geo_indices_num == mesh.n_faces() * 3);

  //   int geo_points_num = 0;
  //   for (const std::vector<double> &geo_comp_points : geo_points) {
  //     geo_points_num += geo_comp_points.size();
  //   }
  //   assert(geo_points_num == mesh.n_vertices() * 3);
  // #endif
  //   ////////////////////////////////////////////////////////////////////////////////////////
  //   // std::vector<double> points;
  //   // std::vector<unsigned> ge_faces;
  //   // bool success = geodesic::read_mesh_from_file(
  //   //     "/home/immocat/Downloads/geodesic/flat_triangular_mesh.txt",
  //   points,
  //   //     ge_faces);
  //   // if (!success) {
  //   //   std::cout << "something is wrong with the input file" << std::endl;
  //   //   return;
  //   // }

  //   // geodesic::Mesh ge_mesh;
  //   // ge_mesh.initialize_mesh_data(
  //   //     points, ge_faces);  // create internal mesh data structure
  //   including
  //   //     edges
  //   ///////////////////////////////////////////////////////////////////////////////////////
  //   std::vector<geodesic::Mesh> geo_meshes(numOfComps);
  //   std::vector<geodesic::GeodesicAlgorithmExact> geo_algorithms;
  //   std::vector<std::vector<geodesic::SurfacePoint>> geo_sources(numOfComps);
  //   //#pragma omp parallel for
  //   for (int i = 0; i < numOfComps; ++i) {
  //     geo_meshes[i].initialize_mesh_data(geo_points[i], geo_indices[i]);
  //     geo_algorithms.emplace_back(&geo_meshes[i]);
  //   }
  //   // add all sources to each mesh
  //   for (int i = 0; i < sample_ids.size(); ++i) {
  //     // get the vertex handle
  //     OpenMesh::VertexHandle vh = vHandles[sample_ids[i]];
  //     int cid = mesh.property(comp_id, vh);
  //     int geo_id = handle_idx_2_geoid[vh.idx()];
  //     geo_sources[cid].emplace_back(&(geo_meshes[cid].vertices()[geo_id]));
  //   }
  // #ifndef NDEBUG
  //   int sourcesCount = 0;
  //   for (int i = 0; i < numOfComps; ++i) {
  //     sourcesCount += geo_sources[i].size();
  //   }
  //   assert(sourcesCount == X.size());
  // #endif
  //   // geo_sources[0].emplace_back(&(geo_meshes[0].vertices()[0]));
  //   // geo_sources[1].emplace_back(&(geo_meshes[1].vertices()[0]));

  //   // propagate over meshes
  //   //#pragma omp parallel for
  //   for (int i = 0; i < numOfComps; ++i) {
  //     geo_algorithms[i].propagate(geo_sources[i]);
  //   }

  // geodesic::GeodesicAlgorithmExact algorithm(&mesh);	//create exact
  // algorithm for the mesh

  // geodesic::Mesh geo_mesh;
  // assert(src_indices.size() == mesh.n_faces() * 3);
  // geo_mesh.initialize_mesh_data(geo_points, src_indices);
  // geodesic::GeodesicAlgorithmExact algorithm(
  //     &geo_mesh);  // create exact algorithm for the mesh

  // // push all sample points(DG nodes) into all_sources
  // std::vector<geodesic::SurfacePoint> all_sources;
  // for (int i = 0; i < sample_ids.size(); ++i) {
  //   geodesic::SurfacePoint source(
  //       &geo_mesh.vertices()[sample_ids[i]]);  // create source
  //   all_sources.emplace_back(source)
  // }
  // printf("\t[DeformationGraph]: start propagate\n");
  // ////////////////////////////////////////////////////////
  // // Perform geodesic algorithm
  // algorithm.propagate(all_sources);  // cover the whole mesh
  // ///////////////////////////////////////////////////////
  // printf("\t[DeformationGraph]: finish propagate, start finding k
  // nearest\n");

  // for (unsigned i = 0; i < geo_mesh.vertices().size(); ++i) {
  //   geodesic::SurfacePoint p(&geo_mesh.vertices()[i]);

  //   double distance;
  //   unsigned best_source = algorithm.best_source(
  //       p, distance);  // for a given surface point, find closets source and
  //                      // distance to this source

  //   std::cout << distance << " ";  // print geodesic distance for every
  //   vertex
  // }

  // create internal mesh data structure including edges
  // For each v_pos, find k = 4(in paper) nearest sources
  // Be careful, there may be less than 4 sources that is visible for v_pot[i]
  // Clean property from last step, Calculate k = 4 weights for each v_pos[i],
  // save in property of vertices of mesh

  // mesh->remove_property(comp_id);
}