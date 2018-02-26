//=============================================================================

#ifndef HOPPE_HH
#define HOPPE_HH

//=============================================================================

#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <vector>
#include <float.h>

//=============================================================================

class ImplicitHoppe
{
public:

    typedef OpenMesh::Vec3f Vec3f;

    // fit RBF to given constraints
    ImplicitHoppe( 
        const std::vector<Vec3f>& _points, 
        const std::vector<Vec3f>& _normals )
        : points_(_points), normals_(_normals)
    {}

    // evaluate implicit at position _p
    float operator()(const Vec3f& _p) const
    {
        float dist(0);

        //////////////////////////////////////////////////////////////////////
        // INSERT CODE:
        // 1) find closest sample point
        // 2) compute distance to its plane
        //--- start strip
        
        
        //////////////////////////////////////////////////////////////////////
        return dist;
    }

private:

    const std::vector<Vec3f>&  points_;
    const std::vector<Vec3f>&  normals_;
};

//=============================================================================
#endif // RBF_HH defined
//=============================================================================
