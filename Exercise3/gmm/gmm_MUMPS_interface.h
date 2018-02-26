// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Generic Matrix Methods  (gmm)
// File    : gmm_MUMPS_interface.h : interface with MUMPS,
//           LU factorization and solve for sparse matrices.
// Date    : December 8, 2005.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2003-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================

/**@file gmm_MUMPS_interface.h
   @brief Interface with MUMPS (LU direct solver for sparse matrices).
*/
#if defined(GMM_USES_MUMPS)

#ifndef GMM_MUMPS_INTERFACE_H
#define GMM_MUMPS_INTERFACE_H

#include <gmm_kernel.h>


extern "C" {

#include <smumps_c.h>
#undef F_INT
#undef F_DOUBLE
#undef F_DOUBLE2
#include <dmumps_c.h>
#undef F_INT
#undef F_DOUBLE
#undef F_DOUBLE2
#include <cmumps_c.h>
#undef F_INT
#undef F_DOUBLE
#undef F_DOUBLE2
#include <zmumps_c.h>
#undef F_INT
#undef F_DOUBLE
#undef F_DOUBLE2

}

namespace gmm {

  template <typename T> struct ij_sparse_matrix {
    std::vector<int> irn;
    std::vector<int> jcn;
    std::vector<T>        a;
    
    template <typename L> void store(const L& l, size_type i) {
       typename linalg_traits<L>::const_iterator it = vect_const_begin(l),
	 ite = vect_const_end(l);
       for (; it != ite; ++it)
	 { irn.push_back(i + 1); jcn.push_back(it.index() + 1); a.push_back(*it); }
    }

    template <typename L> void build_from(const L& l, row_major) {
      for (size_type i = 0; i < mat_nrows(l); ++i)
	store(mat_const_row(l, i), i);
    }

    template <typename L> void build_from(const L& l, col_major) {
      for (size_type i = 0; i < mat_ncols(l); ++i)
	store(mat_const_col(l, i), i);
      irn.swap(jcn);
    }

    template <typename L> ij_sparse_matrix(const L& A) {
      size_type nz = nnz(A);
      irn.reserve(nz); jcn.reserve(nz); a.reserve(nz);
      build_from(A,  typename principal_orientation_type<typename
	       linalg_traits<L>::sub_orientation>::potype());
    }
  };

  /* ********************************************************************* */
  /*   MUMPS solve interface                                               */
  /* ********************************************************************* */


  template <typename T> struct mumps_interf {};

  template <> struct mumps_interf<float> {
    typedef SMUMPS_STRUC_C  MUMPS_STRUC_C;
    typedef float value_type;

    static void mumps_c(MUMPS_STRUC_C &id) { smumps_c(&id); }
  };

  template <> struct mumps_interf<double> {
    typedef DMUMPS_STRUC_C  MUMPS_STRUC_C;
    typedef double value_type;
    static void mumps_c(MUMPS_STRUC_C &id) { dmumps_c(&id); }
  };

  template <> struct mumps_interf<std::complex<float> > {
    typedef CMUMPS_STRUC_C  MUMPS_STRUC_C;
    typedef mumps_complex value_type;
    static void mumps_c(MUMPS_STRUC_C &id) { cmumps_c(&id); }
  };

  template <> struct mumps_interf<std::complex<double> > {
    typedef ZMUMPS_STRUC_C  MUMPS_STRUC_C;
    typedef mumps_double_complex value_type;
    static void mumps_c(MUMPS_STRUC_C &id) { zmumps_c(&id); }
  };


  /** MUMPS solve interface  
   *  Works only with sparse or skyline matrices
   */
  template <typename MAT, typename VECTX, typename VECTB>
  void MUMPS_solve(const MAT &A, const VECTX &X_, const VECTB &B) {
    VECTX &X = const_cast<VECTX &>(X_);

    typedef typename linalg_traits<MAT>::value_type T;
    typedef typename mumps_interf<T>::value_type MUMPS_T;
    if (gmm::mat_nrows(A) != gmm::mat_ncols(A))
      DAL_THROW(failure_error, "Non square matrix");
  
    std::vector<T> rhs(gmm::vect_size(B)); gmm::copy(B, rhs);
  
    ij_sparse_matrix<T> AA(A);
  
    const int JOB_INIT = -1;
    const int JOB_END = -2;
    const int USE_COMM_WORLD = -987654;

    typename mumps_interf<T>::MUMPS_STRUC_C id;

#ifdef GMM_USES_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    
    id.job = JOB_INIT;
    id.par = 1;
    id.sym = 0;
    id.comm_fortran = USE_COMM_WORLD;
    mumps_interf<T>::mumps_c(id);
    
#ifdef GMM_USES_MPI
    if (rank == 0) {
#endif
      id.n = gmm::mat_nrows(A);
      id.nz = AA.irn.size();
      id.irn = &(AA.irn[0]);
      id.jcn = &(AA.jcn[0]);
      id.a = (MUMPS_T*)(&(AA.a[0]));
      id.rhs = (MUMPS_T*)(&(rhs[0]));
#ifdef GMM_USES_MPI
    }
#endif

#define ICNTL(I) icntl[(I)-1]
#define INFO(I) info[(I)-1]
    id.ICNTL(1) = -1; id.ICNTL(2) = -1; id.ICNTL(3) = -1; id.ICNTL(4) = 0;
    id.job = 6;
    mumps_interf<T>::mumps_c(id);
    if (id.INFO(1) < 0) {
      switch (id.INFO(1)) {
	case -6 : case -10 : DAL_THROW(failure_error, "Solve with MUMPS failed: matrix is singular");
	case -13 : DAL_THROW(failure_error, "Solve with MUMPS failed: not enough memory");
	default :  DAL_THROW(failure_error, "Solve with MUMPS failed with error " << id.INFO(1));
      }
    }

    id.job = JOB_END;
    mumps_interf<T>::mumps_c(id);

    gmm::copy(rhs, X);

#undef ICNTL
#undef INFO

  }



  /** MUMPS solve interface for distributed matrices 
   *  Works only with sparse or skyline matrices
   */
  template <typename MAT, typename VECTX, typename VECTB>
  void MUMPS_distributed_matrix_solve(const MAT &A, const VECTX &X_,
				      const VECTB &B) {
    VECTX &X = const_cast<VECTX &>(X_);

    typedef typename linalg_traits<MAT>::value_type T;
    typedef typename mumps_interf<T>::value_type MUMPS_T;
    if (gmm::mat_nrows(A) != gmm::mat_ncols(A))
      DAL_THROW(failure_error, "Non square matrix");
  
    std::vector<T> rhs(gmm::vect_size(B)); gmm::copy(B, rhs);
  
    ij_sparse_matrix<T> AA(A);
  
    const int JOB_INIT = -1;
    const int JOB_END = -2;
    const int USE_COMM_WORLD = -987654;

    typename mumps_interf<T>::MUMPS_STRUC_C id;

#ifdef GMM_USES_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    
    id.job = JOB_INIT;
    id.par = 1;
    id.sym = 0;
    id.comm_fortran = USE_COMM_WORLD;
    mumps_interf<T>::mumps_c(id);
    
    id.n = gmm::mat_nrows(A);
    id.nz_loc = AA.irn.size();
    id.irn_loc = &(AA.irn[0]);
    id.jcn_loc = &(AA.jcn[0]);
    id.a_loc = (MUMPS_T*)(&(AA.a[0]));

#ifdef GMM_USES_MPI
    if (rank == 0) {
#endif
      id.rhs = (MUMPS_T*)(&(rhs[0]));
#ifdef GMM_USES_MPI
    }
#endif

#define ICNTL(I) icntl[(I)-1]
#define INFO(I) info[(I)-1]
    id.ICNTL(1) = -1; id.ICNTL(2) = 6; // id.ICNTL(2) = -1;
    id.ICNTL(3) = 6;
    // id.ICNTL(3) = -1; 
    id.ICNTL(4) = 2;
    id.ICNTL(5)=0; id.ICNTL(18)=3;
    id.job = 6;
    mumps_interf<T>::mumps_c(id);
    if (id.INFO(1) < 0) {
      switch (id.INFO(1)) {
	case -6 : case -10 : DAL_THROW(failure_error, "Solve with MUMPS failed: matrix is singular");
	case -13: DAL_THROW(failure_error, "Solve with MUMPS failed: not enough memory");
	default : DAL_THROW(failure_error, "Solve with MUMPS failed with error " << id.INFO(1));
      }
    }

    id.job = JOB_END;
    mumps_interf<T>::mumps_c(id);
#ifdef GMM_USES_MPI
    MPI_Bcast(&(rhs[0]),id.n,gmm::mpi_type(T()),0,MPI_COMM_WORLD);
#endif
    gmm::copy(rhs, X);

#undef ICNTL
#undef INFO

  }




}

  
#endif // GMM_MUMPS_INTERFACE_H

#endif // GMM_USES_MUMPS
