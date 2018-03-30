//== INCLUDES =================================================================

#include "QuadricT.hh"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <iostream>
#include <set>
//#include <unistd.h> //in linux version this header file was included. Removed for MSVC version
#include <float.h>
#include <vector>

//== IMPLEMENTATION ===========================================================

typedef OpenMesh::TriMesh_ArrayKernelT<>      Mesh;

Mesh                                          mesh;
OpenMesh::VPropHandleT<Quadricd>              vquadric;
OpenMesh::VPropHandleT<float>                 vprio;
OpenMesh::VPropHandleT<Mesh::HalfedgeHandle>  vtarget;

//-----------------------------------------------------------------------------

void  init();
bool  is_collapse_legal(Mesh::HalfedgeHandle _hh);
float priority(Mesh::HalfedgeHandle _heh);
void  decimate(unsigned int _n_vertices);
void  enqueue_vertex(Mesh::VertexHandle vh);

//-----------------------------------------------------------------------------

// access quadric of vertex _vh
Quadricd& quadric(Mesh::VertexHandle _vh) 
{ return mesh.property(vquadric, _vh); }

// access priority of vertex _vh
float& priority(Mesh::VertexHandle _vh)
{ return mesh.property(vprio, _vh); }

// access target halfedge of vertex _vh
Mesh::HalfedgeHandle& target(Mesh::VertexHandle _vh)
{ return mesh.property(vtarget, _vh); }

//-----------------------------------------------------------------------------

// compare functor for priority queue
struct VertexCmp
{
    bool operator()(Mesh::VertexHandle _v0, Mesh::VertexHandle _v1) const
    {
        // std::set needs UNIQUE keys -> handle equal priorities
        return ((priority(_v0) == priority(_v1)) ? 
            (_v0.idx() < _v1.idx()) :
        (priority(_v0) < priority(_v1)));
    }
};

std::set<Mesh::VertexHandle, VertexCmp>  queue;

//-----------------------------------------------------------------------------

int main(int argc, char **argv)
{
    if (argc < 4) 
    {
        std::cerr << "Usage: \n" 
            << argv[0] << " percentage  in.off  out.off\n\n";
        exit(1);
    }

    // add required properties
    mesh.request_vertex_status();
    mesh.request_edge_status();
    mesh.request_face_status();
    mesh.request_face_normals();
    mesh.add_property(vquadric);
    mesh.add_property(vprio);
    mesh.add_property(vtarget);


    // read mesh
    OpenMesh::IO::read_mesh(mesh, argv[2]);
    std::cout << "#vertices: " << mesh.n_vertices() << std::endl;

    // compute normals & quadrics
    init();

    // decimate
    decimate((int)(atof(argv[1])*mesh.n_vertices()));
    std::cout << "#vertices: " << mesh.n_vertices() << std::endl;

    // write mesh
    OpenMesh::IO::write_mesh(mesh, argv[3]);
}

//-----------------------------------------------------------------------------

void init()
{
    // compute face normals
    mesh.update_face_normals();

    Mesh::VertexIter  v_it, v_end = mesh.vertices_end();
    Mesh::Point       n;
    double              a, b, c, d;

    for (v_it=mesh.vertices_begin(); v_it != v_end; ++v_it)
    {
        priority(*v_it) = -1.0;
        quadric(*v_it).clear();

        // Exercise 5.1 --------------------------------------
        // INSERT CODE:
        // calc vertex quadrics from incident triangles
        // ---------------------------------------------------
        Mesh::Point v_p = mesh.point(*v_it);
        for(Mesh::VertexFaceIter vf_it = mesh.vf_begin(*v_it); vf_it.is_valid(); ++vf_it){
            // iterate all neighbor faces
            Mesh::Point face_normal = mesh.normal(*vf_it);
            float d = -OpenMesh::dot(face_normal, v_p);
            quadric(*v_it) = Quadricd(face_normal[0], face_normal[1], face_normal[2], d);
        }
    }
}

//-----------------------------------------------------------------------------

bool
is_collapse_legal(Mesh::HalfedgeHandle _hh)
{
    // collect vertices
    Mesh::VertexHandle v0, v1;
    v0 = mesh.from_vertex_handle(_hh);
    v1 = mesh.to_vertex_handle(_hh);

    // collect faces
    Mesh::FaceHandle fl = mesh.face_handle(_hh);
    Mesh::FaceHandle fr = mesh.face_handle(mesh.opposite_halfedge_handle(_hh));

    // backup point positions
    // Mesh::Point p0 = mesh.point(v0);
    Mesh::Point p1 = mesh.point(v1);

    // topological test
    if (!mesh.is_collapse_ok(_hh))
        return false;

    // test boundary stuff
    if (mesh.is_boundary(v0) && !mesh.is_boundary(v1))
        return false;

    // Exercise 5.2 -----------------------------------------------
    // INSERT CODE:
    // test normal flipping:
    //   if normal vector of a (non-degenerate) triangle changes by 
    //   more than pi/4 degrees, return false.
    // ------------------------------------------------------------
    static const float cos_piD4 = 0.70710678118;
    for(Mesh::VertexFaceIter v0f_it = mesh.vf_begin(v0); v0f_it.is_valid(); ++v0f_it){
        if((*v0f_it).idx() == fl.idx() || (*v0f_it).idx() == fr.idx())
            continue;
        // for each face, iterate all vertecies with ccw fashion
        Mesh::VertexHandle fv[3];
        Mesh::Point fv_p[3];
        int i = 0;
        for(Mesh::FaceVertexCCWIter v0fv_it = mesh.fv_ccwbegin(*v0f_it); v0fv_it.is_valid(); ++v0fv_it){
            assert(i < 3);
            fv[i] = *v0fv_it;
            fv_p[i] = mesh.point(*v0fv_it);
            ++i;
        }
        Mesh::Point origin_normal, changed_normal;
        if(fv[0].idx() == v0.idx()){
            origin_normal = OpenMesh::cross(fv_p[1] - fv_p[2], fv_p[0] - fv_p[2]).normalized();
            changed_normal = OpenMesh::cross(fv_p[1] - fv_p[2], p1 - fv_p[2]).normalized();
        }else if(fv[1].idx() == v0.idx()){
            origin_normal = OpenMesh::cross(fv_p[2] - fv_p[0], fv_p[1] - fv_p[0]).normalized();
            changed_normal = OpenMesh::cross(fv_p[2] - fv_p[0], p1 - fv_p[0]).normalized();
        }else{
            origin_normal = OpenMesh::cross(fv_p[0] - fv_p[1], fv_p[2] - fv_p[1]).normalized();
            changed_normal = OpenMesh::cross(fv_p[0] - fv_p[1], p1 - fv_p[1]).normalized();
        }
        if(OpenMesh::dot(origin_normal, changed_normal) < cos_piD4)
            return false;
    }
    // collapse passed all tests -> ok
    return true;
}

//-----------------------------------------------------------------------------

float priority(Mesh::HalfedgeHandle _heh)
{
    // Exercise 5.3 ----------------------------------------------
    // INSERT CODE:
    // return priority: the smaller the better
    // use quadrics to estimate approximation error
    // -----------------------------------------------------------
    Mesh::VertexHandle v1 = mesh.from_vertex_handle(_heh), v2 = mesh.to_vertex_handle(_heh);
    // v2 as vbar
    Quadricd q12 = quadric(v1);
    q12 += quadric(v2);
    return q12(mesh.point(v2));
    
}

//-----------------------------------------------------------------------------

void enqueue_vertex(Mesh::VertexHandle _vh)
{
    float                   prio, min_prio(FLT_MAX);
    Mesh::HalfedgeHandle  min_hh;


    // find best out-going halfedge
    for (Mesh::VOHIter vh_it(mesh, _vh); vh_it; ++vh_it)
    {
        if (is_collapse_legal(*vh_it))
        {
            prio = priority(*vh_it);
            if (prio != -1.0 && prio < min_prio)
            {
                min_prio = prio;
                min_hh   = *vh_it;
            }
        }
    }

    // update queue
    if (priority(_vh) != -1.0) 
    {
        queue.erase(_vh);
        priority(_vh) = -1.0;
    }

    if (min_hh.is_valid()) 
    {
        priority(_vh) = min_prio;
        target(_vh)   = min_hh;
        queue.insert(_vh);
    }
}

//-----------------------------------------------------------------------------

void  decimate(unsigned int _n_vertices)
{
    unsigned int nv(mesh.n_vertices());

    Mesh::HalfedgeHandle hh;
    Mesh::VertexHandle   to, from;
    Mesh::VVIter         vv_it;

    std::vector<Mesh::VertexHandle>            one_ring;
    std::vector<Mesh::VertexHandle>::iterator  or_it, or_end;

    // build priority queue
    Mesh::VertexIter  v_it  = mesh.vertices_begin(), 
        v_end = mesh.vertices_end();

    queue.clear();
    for (; v_it!=v_end; ++v_it)
        enqueue_vertex(*v_it);

    while (nv > _n_vertices && !queue.empty())
    {
        // Exercise 5.4 ----------------------------------------------
        // INSERT CODE:
        // Decimate using priority queue:
        //   1) take 1st element of queue
        //   2) collapse this halfedge
        //   3) update queue
        // -----------------------------------------------------------
        //   1) take 1st element of queue
        Mesh::VertexHandle v_min = *(queue.begin());
        Mesh::HalfedgeHandle he_min = mesh.property(vtarget, v_min);
        Mesh::VertexHandle v_min_to = mesh.to_vertex_handle(he_min);
        assert(mesh.from_vertex_handle(he_min).idx() == v_min.idx());
        //   update quadric error Qbar
        quadric(v_min_to) += quadric(v_min);
        std::vector<Mesh::VertexHandle> enqueue_candidate;
        for(Mesh::VertexVertexIter vv_it = mesh.vv_begin(v_min); vv_it.is_valid(); ++vv_it){
            enqueue_candidate.push_back(*vv_it);
        }
        mesh.collapse(he_min);
        mesh.adjust_outgoing_halfedge(v_min_to);
        queue.erase(v_min);
        for(Mesh::VertexHandle vh : enqueue_candidate){
            enqueue_vertex(vh);
        }
        --nv;
    }

    // clean up
    queue.clear();

    // now, delete the items marked to be deleted
    mesh.garbage_collection();
}

//-----------------------------------------------------------------------------
