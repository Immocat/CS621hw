//=============================================================================

#include "ImplicitRBF.hh"

//== IMPLEMENTATION ==========================================================

ImplicitRBF::ImplicitRBF( 
    const std::vector<Vec3f>& _points, 
    const std::vector<Vec3f>& _normals )
{
    //////////////////////////////////////////////////////////////////////
    // INSERT CODE:
    // 1) collect constraints (on-surface and off-surface)
    // get bounding box, recalculate because do not want to change reconstruct.cc
    static float dis_threshold = 0.01;
    weights_.clear();
    centers_.clear();
    Vec3f bb_min( _points[0]), bb_max( _points[0]);

    for (unsigned int i=1; i<_points.size(); ++i)
    {
        bb_min.minimize( _points[i] );
        bb_max.maximize( _points[i] );
    }
    double eps  = (bb_max-bb_min).norm() * dis_threshold;

    // copy centers
    int N = _points.size();
    int NpN = N + N;
    centers_.reserve(NpN);
    weights_.assign(NpN, 0);
    for(int i = 0; i < N; ++i){
        centers_.emplace_back(Vec3d(_points[i][0], _points[i][1], _points[i][2]));
        //weights_.push_back(0.0);
    }
    for(int i = 0; i < N; ++i){
        Vec3f newPos = _points[i] + eps * _normals[i];
        centers_.emplace_back(Vec3d(newPos[0], newPos[1], newPos[2]));
        //weights_.push_back(0.0);
    }
    // 2) setup matrix
    
    gmmMatrix M(NpN, NpN);
    gmmVector d(NpN, 0.0);

    for(int i = 0; i < NpN; ++i){
        for(int j = 0; j < NpN; ++j){
            M(i, j) = kernel(centers_[i], centers_[j]);
        }
    }
    for(int i = N; i < NpN; ++i){
        d[i] = eps;
    }
    // for(int i = N; i < NpN; ++i){
    //     for(int j = 0; j < NpN; ++j){
    //         M(i, j) = kernel((Vec3d)(_points[i] + eps * _normals[i]), (Vec3d)_points[j]); 
    //     }
    // }
    // 3) solve linear system for weights_
    

    solve_linear_system(M, d, weights_);
    
    
    //////////////////////////////////////////////////////////////////////
}

//-----------------------------------------------------------------------------

void ImplicitRBF::solve_linear_system( 
    gmmMatrix& _M, 
    gmmVector& _b, 
    gmmVector& _x )
{
    // solve linear system by gmm's LU factorization
    unsigned int N = _b.size();
    _x.resize(N);
    std::vector< size_t >  ipvt(N);
    gmm::lu_factor( _M, ipvt );
    gmm::lu_solve( _M, ipvt, _x, _b );
}

//-----------------------------------------------------------------------------

double ImplicitRBF::operator()(const Vec3f& _p) const
{
    std::vector<Vec3d>::const_iterator  
        c_it(centers_.begin()),
        c_end(centers_.end());

    std::vector<double>::const_iterator   
        w_it(weights_.begin());

    const Vec3d p(_p);
    double f(0);

    for (; c_it!=c_end; ++c_it, ++w_it)
        f += *w_it * kernel(*c_it, p);

    return f;
}

//=============================================================================
