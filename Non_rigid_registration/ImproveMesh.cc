#include <EventList.hh>
#include <RegistrationViewer.hh>
#include <cassert>
static void split_edge(OpenMesh::TriMesh_ArrayKernelT<>& mesh,
                       const OpenMesh::EdgeHandle& eh,
                       const OpenMesh::Vec3f& newPoint, EventList* eventList,
                       unsigned int frame_id) {
  // given mesh and edge, split it at newPoint
  // if eventList is not nullptr, add a new event to eventlist

  // not split boundary edge to keep topology
  if (mesh.is_boundary(eh)) return;

  OpenMesh::HalfedgeHandle h0 = mesh.halfedge_handle(eh, 0),
                           h1 = mesh.halfedge_handle(eh, 1),
                           h2 = mesh.next_halfedge_handle(h0),
                           h3 = mesh.next_halfedge_handle(h2),
                           h4 = mesh.next_halfedge_handle(h1),
                           h5 = mesh.next_halfedge_handle(h4);

  OpenMesh::VertexHandle v0 = mesh.from_vertex_handle(h0),
                         v1 = mesh.to_vertex_handle(h0),
                         v2 = mesh.to_vertex_handle(h2),
                         v3 = mesh.to_vertex_handle(h4);
  // assert(v1 == mesh.from_vertex_handle(h1));
  // assert(mesh.from_vertex_handle(h2) == v1);
  // assert(mesh.from_vertex_handle(h3) == v2);
  // assert(mesh.to_vertex_handle(h3) == v0);
  // assert(mesh.from_vertex_handle(h4) == v0);
  // assert(mesh.from_vertex_handle(h5) == v3);
  // assert(mesh.to_vertex_handle(h5) == v1);
  // assert(mesh.face_handle(h0) == mesh.face_handle(h2));
  // assert(mesh.face_handle(h2) == mesh.face_handle(h3));

  OpenMesh::VertexHandle vh = mesh.new_vertex(newPoint);
  // when set_vertex_handle is called, the opposite half edge's from vertex will
  // set automatically
  mesh.set_vertex_handle(h0, vh);
  assert(mesh.to_vertex_handle(h0) == vh);
  assert(mesh.from_vertex_handle(h1) == vh);
  assert(mesh.to_vertex_handle(h1) == v0);

  OpenMesh::HalfedgeHandle hh2 = mesh.new_edge(vh, v2);
  mesh.set_face_handle(hh2, mesh.face_handle(h3));
  mesh.set_next_halfedge_handle(hh2, h3);
  mesh.set_next_halfedge_handle(h0, hh2);
  assert(mesh.next_halfedge_handle(h3) == h0);

  OpenMesh::HalfedgeHandle h2h = mesh.opposite_halfedge_handle(hh2);
  // mesh.set_next_halfedge_handle(hh2, h2h);
  assert(mesh.from_vertex_handle(h2h) == v2);
  assert(mesh.to_vertex_handle(h2h) == vh);
  OpenMesh::FaceHandle f2 = mesh.new_face();
  OpenMesh::HalfedgeHandle hh1 = mesh.new_edge(vh, v1);

  ///
  // OpenMesh::VertexHandle v1test = mesh.from_vertex_handle(h1);
  ///

  mesh.set_halfedge_handle(f2, hh1);
  mesh.set_face_handle(hh1, f2);
  mesh.set_face_handle(h2, f2);
  mesh.set_face_handle(h2h, f2);
  mesh.set_next_halfedge_handle(hh1, h2);
  mesh.set_next_halfedge_handle(h2h, hh1);
  mesh.set_next_halfedge_handle(h2, h2h);

  if (v2 == v3) {
    printf("what?\n");
    return;
  }
  OpenMesh::HalfedgeHandle h3h = mesh.new_edge(v3, vh);
  // OpenMesh::FaceHandle f4 = mesh.face_handle(h4);
  mesh.set_face_handle(h3h, mesh.face_handle(h4));
  mesh.set_next_halfedge_handle(h3h, h1);
  mesh.set_next_halfedge_handle(h1, h4);
  mesh.set_next_halfedge_handle(h4, h3h);

  OpenMesh::HalfedgeHandle hh3 = mesh.opposite_halfedge_handle(h3h);
  OpenMesh::FaceHandle f5 = mesh.new_face();
  OpenMesh::HalfedgeHandle h1h = mesh.opposite_halfedge_handle(hh1);
  mesh.set_halfedge_handle(f5, h1h);
  mesh.set_face_handle(h1h, f5);
  mesh.set_face_handle(hh3, f5);
  mesh.set_face_handle(h5, f5);

  mesh.set_next_halfedge_handle(h1h, hh3);
  mesh.set_next_halfedge_handle(hh3, h5);
  mesh.set_next_halfedge_handle(h5, h1h);

  mesh.set_halfedge_handle(vh, hh2);

  // // Never forget this, when playing with the topology
  mesh.adjust_outgoing_halfedge(v0);
  mesh.adjust_outgoing_halfedge(v1);
  mesh.adjust_outgoing_halfedge(v2);
  mesh.adjust_outgoing_halfedge(v3);
  mesh.adjust_outgoing_halfedge(vh);
  // for (OpenMesh::TriMesh_ArrayKernelT<>::VertexFaceIter vf_it =
  //          mesh.vf_iter(vh);
  //      vf_it.is_valid(); ++vf_it) {
  //   std::cout << "vf test: " << *vf_it << std::endl;
  // }
  // for (OpenMesh::TriMesh_ArrayKernelT<>::VertexVertexIter vv_it =
  //          mesh.vv_iter(vh);
  //      vv_it.is_valid(); ++vv_it) {
  //   std::cout << "vv test: " << *vv_it << std::endl;
  // }
  // is eventList is valid, add new event to frame S_id - 2, since S_id is point
  // to next frame's mesh
  if (eventList) {
    Event newEvent;
    newEvent.etype = Event::EType::SPLIT_EDGE;
    newEvent.start[0] = v0.idx();
    newEvent.start[1] = v1.idx();
    newEvent.end = vh.idx();
    newEvent.alpha[0] = 0.5;
    eventList->events_at_frame[frame_id].emplace_back(newEvent);
  }
}
static void collapse_edge(OpenMesh::TriMesh_ArrayKernelT<>& mesh,
                          const OpenMesh::EdgeHandle& eh, EventList* eventList,
                          unsigned int frame_id) {
  //
  // not split boundary edge to keep topology
  // if (mesh.is_boundary(eh)) return;

  OpenMesh::HalfedgeHandle h0 = mesh.halfedge_handle(eh, 0),
                           h1 = mesh.halfedge_handle(eh, 1),
                           h2 = mesh.next_halfedge_handle(h0),
                           h4 = mesh.next_halfedge_handle(h1);

  OpenMesh::VertexHandle v0 = mesh.from_vertex_handle(h0),
                         v1 = mesh.to_vertex_handle(h0),
                         v2 = mesh.to_vertex_handle(h2),
                         v3 = mesh.to_vertex_handle(h4);

  OpenMesh::Vec3f midPoint = 0.5 * (mesh.point(v0) + mesh.point(v1));
  mesh.set_point(v0, midPoint);
  mesh.collapse(h1);
  mesh.adjust_outgoing_halfedge(v0);
  mesh.adjust_outgoing_halfedge(v2);
  mesh.adjust_outgoing_halfedge(v3);
  // is eventList is valid, add new event to frame S_id - 2, since S_id is point
  // to next frame's mesh
  if (eventList) {
    Event newEvent;
    newEvent.etype = Event::EType::COLLAPSE_EDGE;
    newEvent.start[0] = v0.idx();
    newEvent.start[1] = v1.idx();
    // v0 is moved to (v0 + v1)/2, v1 is gone
    newEvent.alpha[0] = 0.5;
    eventList->events_at_frame[frame_id].emplace_back(newEvent);
  }
}
//
// improve mesh triangle mesh mesh with split and collapse, save operation to
// eventList
void RegistrationViewer::improveMesh(Mesh& mesh, EventList* eventList,
                                     std::vector<unsigned int>& indices) {
  // 1. interate all edges, calculate average edge distance
  // OpenMesh::EPropHandleT<double> norm;
  // mesh.add_property(norm);
  //
  int edge_long = 0, edge_short = 0;
  float ave_edge_norm = 0;
  int edge_num = 0;
  for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end();
       ++e_it, ++edge_num) {
    OpenMesh::Vec3f p =
        mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(e_it, 0)));
    OpenMesh::Vec3f q =
        mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(e_it, 0)));
    float edgeLength = (p - q).norm();
    // mesh.property(norm)
    ave_edge_norm += edgeLength;
  }
  if (edge_num == 0) {
    printf("[Improve Mesh]: edge number is 0, return.\n");
    return;
  }
  ave_edge_norm /= (float)edge_num;
  // 2. calculate maximum and minimum edge length
  // min = 0.5 * ave, max = 1.5 * ave
  float min_edge_norm = 0.5 * ave_edge_norm;
  float min_edge_norm2 = min_edge_norm * min_edge_norm;
  float max_edge_norm = 1.5 * ave_edge_norm;
  float max_edge_norm2 = max_edge_norm * max_edge_norm;

  // 3. loop over all edges again, find invalid edges, modify mesh, record in
  // eventList
  Mesh::EdgeIter e_end = mesh.edges_end();

  for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != e_end; ++e_it) {
    HalfedgeHandle he = mesh.halfedge_handle(*e_it, 0);
    OpenMesh::VertexHandle v0 = mesh.to_vertex_handle(he);
    OpenMesh::VertexHandle v1 = mesh.from_vertex_handle(he);

    OpenMesh::Vec3f p = mesh.point(v0);
    OpenMesh::Vec3f q = mesh.point(v1);
    float edgeLength2 = (p - q).sqrnorm();
    if (edgeLength2 > max_edge_norm2) {
      // need edge split
      Mesh::Point newPoint = (p + q) * 0.5;
      split_edge(mesh, *e_it, newPoint, eventList, S_id - 2);
      // mesh.garbage_collection();
      edge_long++;
    }
  }
  e_end = mesh.edges_end();
  for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != e_end; ++e_it) {
    HalfedgeHandle he = mesh.halfedge_handle(*e_it, 0);
    OpenMesh::VertexHandle v0 = mesh.to_vertex_handle(he);
    OpenMesh::VertexHandle v1 = mesh.from_vertex_handle(he);

    OpenMesh::Vec3f p = mesh.point(v0);
    OpenMesh::Vec3f q = mesh.point(v1);
    float edgeLength2 = (p - q).sqrnorm();
    OpenMesh::FaceHandle f0 = mesh.face_handle(he);
    OpenMesh::FaceHandle f1 =
        mesh.face_handle(mesh.opposite_halfedge_handle(he));
    // OpenMesh::Vec3f normal0 = mesh.normal(f0);
    // OpenMesh::Vec3f normal1 = mesh.normal(f1);
    // dot product
    // float dihedral_cos_abs = std::abs(normal0 | normal1);
    if ((edgeLength2 < min_edge_norm2) && mesh.is_collapse_ok(he)) {
      // need edge collapse
      collapse_edge(mesh, *e_it, eventList, S_id - 2);
      edge_short++;
    }
  }
  mesh.garbage_collection();

  // check face inner angle, must greater than 10 degree
  Mesh::FaceIter f_end = mesh.faces_end();
  for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != f_end; ++f_it) {
    // got three vertices
    int ver_count = 0;
    OpenMesh::VertexHandle vers[3];
    OpenMesh::HalfedgeHandle hals[3];
    for (Mesh::FaceHalfedgeCCWIter fh_it = mesh.fh_ccwbegin(*f_it);
         fh_it.is_valid(); ++fh_it) {
      hals[ver_count++] = *fh_it;
    }
    assert(ver_count == 3);
    vers[0] = mesh.from_vertex_handle(hals[0]);
    vers[1] = mesh.to_vertex_handle(hals[0]);
    vers[2] = mesh.to_vertex_handle(hals[1]);
    // calculate inner angle
    OpenMesh::Vec3f p0 = mesh.point(vers[0]), p1 = mesh.point(vers[1]),
                    p2 = mesh.point(vers[2]);
    if (((p1 - p0) | (p2 - p0)) > 0.98480775301) {
      // inner angle less than 10 degree
      // collapse p1p2
      collapse_edge(mesh, mesh.edge_handle(hals[1]), eventList, S_id - 2);
      continue;
    }
    if (((p0 - p1) | (p2 - p1)) > 0.98480775301) {
      // inner angle less than 10 degree
      // collapse p1p2
      collapse_edge(mesh, mesh.edge_handle(hals[2]), eventList, S_id - 2);
      continue;
    }
    if (((p0 - p2) | (p1 - p2)) > 0.98480775301) {
      // inner angle less than 10 degree
      // collapse p1p2
      collapse_edge(mesh, mesh.edge_handle(hals[0]), eventList, S_id - 2);
      continue;
    }
  }
  mesh.garbage_collection();

  // normal updated here
  mesh.update_normals();

  // update indices for rendering
  update_face_indices(mesh, indices);
  printf("[ImproveMesh]: split %d edges, merge %d edges.\n", edge_long,
         edge_short);
}
