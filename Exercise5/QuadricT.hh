//=============================================================================
//                                                                            
//   Example code for the full-day course
//
//   M. Botsch, M. Pauly, C. Roessl, S. Bischoff, L. Kobbelt,
//   "Geometric Modeling Based on Triangle Meshes"
//   held at SIGGRAPH 2006, Boston, and Eurographics 2006, Vienna.
//
//   Copyright (C) 2006 by  Computer Graphics Laboratory, ETH Zurich, 
//                      and Computer Graphics Group,      RWTH Aachen
//
//                                                                            
//-----------------------------------------------------------------------------
//                                                                            
//                                License                                     
//                                                                            
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License
//   as published by the Free Software Foundation; either version 2
//   of the License, or (at your option) any later version.
//   
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//   
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin Street, Fifth Floor, 
//   Boston, MA  02110-1301, USA.
//                                                                            
//=============================================================================
//=============================================================================
//
//  CLASS QuadricT
//
//=============================================================================

#ifndef QUADRIC_HH
#define QUADRIC_HH

//== INCLUDES =================================================================

#include <OpenMesh/Core/Geometry/VectorT.hh>


//== CLASS DEFINITION =========================================================


/** /class QuadricT

    Stores a quadric as a 4x4 symmetrix matrix. Used by the
    error quadric based mesh decimation algorithms.
**/

template <class Scalar>  
class QuadricT
{
public:
  

  /// construct with upper triangle of symmetrix 4x4 matrix
  QuadricT(Scalar _a, Scalar _b, Scalar _c, Scalar _d,
	              Scalar _e, Scalar _f, Scalar _g,
	                         Scalar _h, Scalar _i,
	                                    Scalar _j)
    : a(_a), b(_b), c(_c), d(_d), 
             e(_e), f(_f), g(_g),    
                    h(_h), i(_i), 
                           j(_j)
  {}


  /// constructor from given plane equation: ax+by+cz+d=0
  QuadricT( Scalar _a=0.0, Scalar _b=0.0, Scalar _c=0.0, Scalar _d=0.0 )
    :  a(_a*_a), b(_a*_b),  c(_a*_c),  d(_a*_d),
                 e(_b*_b),  f(_b*_c),  g(_b*_d),
                            h(_c*_c),  i(_c*_d),
                                       j(_d*_d)
  {}


  /// set all entries to zero
  void clear()  { a = b = c = d = e = f = g = h = i = j = 0.0; }

  
  /// add quadrics
  QuadricT<Scalar>& operator+=( const QuadricT<Scalar>& _q )
  {
    a += _q.a;  b += _q.b;  c += _q.c;  d += _q.d;
                e += _q.e;  f += _q.f;  g += _q.g;
                            h += _q.h;  i += _q.i;
			                j += _q.j;
    return *this;
  }


  /// multiply by scalar
  QuadricT<Scalar>& operator*=( Scalar _s)
  {
    a *= _s;  b *= _s;  c *= _s;  d *= _s;
              e *= _s;  f *= _s;  g *= _s;
                        h *= _s;  i *= _s;
                                  j *= _s;
    return *this;
  }


  /// evaluate quadric Q at vector v: v*Q*v
  template <typename T>
  Scalar operator()(const OpenMesh::VectorT<T,3> _v) const
  {
    Scalar x(_v[0]), y(_v[1]), z(_v[2]);
    return a*x*x + 2.0*b*x*y + 2.0*c*x*z + 2.0*d*x
                 +     e*y*y + 2.0*f*y*z + 2.0*g*y
	                     +     h*z*z + 2.0*i*z
                                         +     j;
  }



private:

  Scalar a, b, c, d, 
            e, f, g, 
               h, i, 
                  j;
};


/// Quadric using floats
typedef QuadricT<float> Quadricf;

/// Quadric using double
typedef QuadricT<double> Quadricd;


//=============================================================================
#endif // QUADRIC_HH defined
//=============================================================================
