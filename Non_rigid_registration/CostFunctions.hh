#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "Vector.hh"
// putting all Energy calculation structs here.
struct Dot3dResidual {
  Dot3dResidual() {}
  template <typename T>
  bool operator()(const T* const a, const T* const b, T* residual) const {
    residual[0] = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    return true;
  }
};

struct OneMinusSelfDotResidual {
  OneMinusSelfDotResidual() {}
  template <typename T>
  bool operator()(const T* const a, T* residual) const {
    residual[0] = (T)1 - (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    return true;
  }
};

struct SmoothResidual {
  SmoothResidual(const Vector3d& xi, const Vector3d& xj) : xj_m_xi(xj - xi) {}
  template <typename T>
  bool operator()(const T* const c1, const T* const c2, const T* const c3,
                  const T* const bi, const T* const bj, T* residual) const {
    T vec[3];
    vec[0] = c1[0] * xj_m_xi[0] + c2[0] * xj_m_xi[1] + c3[0] * xj_m_xi[2] -
             xj_m_xi[0] + bi[0] - bj[0];
    vec[1] = c1[1] * xj_m_xi[0] + c2[1] * xj_m_xi[1] + c3[1] * xj_m_xi[2] -
             xj_m_xi[1] + bi[1] - bj[1];
    vec[2] = c1[2] * xj_m_xi[0] + c2[2] * xj_m_xi[1] + c3[2] * xj_m_xi[2] -
             xj_m_xi[2] + bi[2] - bj[2];
    residual[0] = vec[0] * vec[0] + vec[1] * vec[1] +
                  vec[2] * vec[2];  // sqrt(std::max((T)0, ));
    return true;
  }

 private:
  const Vector3d xj_m_xi;
};

struct PointResidual {
  PointResidual(const Vector3d& xii, const Vector3d& cii) : xi(xii), ci(cii) {}
  template <typename T>
  bool operator()(const T* const c1, const T* const c2, const T* const c3,
                  const T* const bi, T* residual) const {
    T vec[3];
    vec[0] = c1[0] * xi[0] + c2[0] * xi[1] + c3[0] * xi[2] + bi[0] - ci[0];
    vec[1] = c1[1] * xi[0] + c2[1] * xi[1] + c3[1] * xi[2] + bi[1] - ci[1];
    vec[2] = c1[2] * xi[0] + c2[2] * xi[1] + c3[2] * xi[2] + bi[2] - ci[2];
    residual[0] = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
    return true;
  }

 private:
  const Vector3d& xi;
  const Vector3d& ci;
};
struct PointResidualFine {
  PointResidualFine(const OpenMesh::Vec3f& xii, const Vector3d& cii)
      : xi(xii), ci(cii) {}
  template <typename T>
  bool operator()(const T* const A, const T* const bi, T* residual) const {
    T vec[3];
    vec[0] = A[0] * (T)xi[0] + A[3] * (T)xi[1] + A[6] * (T)xi[2] + bi[0] - ci[0];
    vec[1] = A[1] * (T)xi[0] + A[4] * (T)xi[1] + A[7] * (T)xi[2] + bi[1] - ci[1];
    vec[2] = A[2] * (T)xi[0] + A[5] * (T)xi[1] + A[8] * (T)xi[2] + bi[2] - ci[2];
    residual[0] = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
    return true;
  }

 private:
  const OpenMesh::Vec3f& xi;
  const Vector3d& ci;
};
struct PlaneResidual {
  PlaneResidual(const Vector3d& xii, const Vector3d& cii, const Vector3d& nii)
      : xi(xii), ci(cii), ni(nii) {}
  template <typename T>
  bool operator()(const T* const c1, const T* const c2, const T* const c3,
                  const T* const bi, T* residual) const {
    T vec[3];
    vec[0] = c1[0] * xi[0] + c2[0] * xi[1] + c3[0] * xi[2] + bi[0] - ci[0];
    vec[1] = c1[1] * xi[0] + c2[1] * xi[1] + c3[1] * xi[2] + bi[1] - ci[1];
    vec[2] = c1[2] * xi[0] + c2[2] * xi[1] + c3[2] * xi[2] + bi[2] - ci[2];
    residual[0] = (vec[0] * ni[0] + vec[1] * ni[1] + vec[2] * ni[2]);
    return true;
  }

 private:
  const Vector3d& xi;
  const Vector3d& ci;
  const Vector3d& ni;
};
struct RigidResidualFine {
  RigidResidualFine(const OpenMesh::Vec3f& x00, const OpenMesh::Vec3f& x11,
                    const Vector3d& c00, const Vector3d& c11)
      : x0(x00), x1(x11), c0(c00), c1(c11) {}
  template <typename T>
  bool operator()(const T* const A0, const T* const b0, const T* const A1,
                  const T* const b1, T* residual) const {
    T veca[3];
    veca[0] = A0[0] * (T)x0[0] + A0[3] * (T)x0[1] + A0[6] * (T)x0[2] + b0[0] - c0[0];
    veca[1] = A0[1] * (T)x0[0] + A0[4] * (T)x0[1] + A0[7] * (T)x0[2] + b0[1] - c0[1];
    veca[2] = A0[2] * (T)x0[0] + A0[5] * (T)x0[1] + A0[8] * (T)x0[2] + b0[2] - c0[2];
              
    T vecb[3];
    vecb[0] = A1[0] * (T)x1[0] + A1[3] * (T)x1[1] + A1[6] * (T)x1[2] + b1[0] - c1[0];
    vecb[1] = A1[1] * (T)x1[0] + A1[4] * (T)x1[1] + A1[7] * (T)x1[2] + b1[1] - c1[1];
    vecb[2] = A1[2] * (T)x1[0] + A1[5] * (T)x1[1] + A1[8] * (T)x1[2] + b1[2] - c1[2];
    residual[0] = abs(veca[0] * veca[0] + veca[1] * veca[1] + veca[2] * veca[2]) 
    -(vecb[0] * vecb[0] + vecb[1] * vecb[1] + vecb[2] * vecb[2])  ;
    return true;
  }

 private:
  const OpenMesh::Vec3f& x0;
  const OpenMesh::Vec3f& x1;
  const Vector3d& c0;
  const Vector3d& c1;
};