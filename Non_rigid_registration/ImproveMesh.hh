// Chapter 4.1 of "Tracking Surfaces with Evolving Topology"
// Also in Chapter 4 of "Liquid Simulation With Mesh-Based Surface Traking"

#pragma once
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
// class OpenMesh::TriMesh_ArrayKernelT<>;
class EventList;
// improve mesh triangle mesh M with split and collapse, save operation to
// eventList
void improveMesh(OpenMesh::TriMesh_ArrayKernelT<> &M, EventList &eventList);
