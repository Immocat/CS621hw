#pragma once
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
    residual[0] = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];//sqrt(std::max((T)0, ));
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

struct PlaneResidual {
  PlaneResidual(const Vector3d& xii, const Vector3d& cii, const Vector3d& nii) : xi(xii), ci(cii),ni(nii) {}
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