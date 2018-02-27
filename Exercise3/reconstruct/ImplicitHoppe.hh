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
        //float dist(0);

        //////////////////////////////////////////////////////////////////////
        // INSERT CODE:
        // 1) find closest sample point
        float min_length = FLT_MAX;
        int min_id = -1;
        for(int i = 0; i < points_.size(); ++i){
            float dist = (_p - points_[i]).sqrnorm();
            if(dist < min_length){
                min_length = dist;
                min_id = i;
            }
        }
        // 2) compute distance to its plane
        const Vec3f& xi(points_[min_id]);
        const Vec3f& ni(normals_[min_id]);
        float dist = dot(_p - xi, ni);
        
        
        
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
