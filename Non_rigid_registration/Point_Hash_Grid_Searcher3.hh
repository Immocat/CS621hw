#pragma once
#include "Vector.hh"
#include <vector>


//!
//! \brief Hash grid-based 3-D point searcher.
//!
//! This class implements 3-D point searcher by using hash grid for its internal
//! acceleration data structure. Each point is recorded to its corresponding
//! bucket where the hashing function is 3-D grid mapping.
//!
class PointHashGridSearcher3{
 public:




  //!
  //! Returns true if there are any nearby points for given origin within
  //! radius.
  //!
  //! \param[in]  origin The origin.
  //! \param[in]  radius The radius.
  //!
  //! \return     True if has nearby point, false otherwise.
  //!
  bool hasNearbyPoint(const Vector3d& origin, double radius) const;


  //! Builds internal acceleration structure for given points list.
  void build(const std::vector<Vector3d>& points);
  //!
  //! \brief      Adds a single point to the hash grid.
  //!
  //! This function adds a single point to the hash grid for future queries.
  //! It can be used for a hash grid that is already built by calling function
  //! PointHashGridSearcher3::build.
  //!
  //! \param[in]  point The point to be added.
  //!
  void add(const Vector3d& point);

  //!
  //! \brief      Returns the internal bucket.
  //!
  //! A bucket is a list of point indices that has same hash value. This
  //! function returns the (immutable) internal bucket structure.
  //!
  //! \return     List of buckets.
  //!
  const std::vector<std::vector<size_t>>& buckets() const;

  //!
  //! Returns the hash value for given 3-D bucket index.
  //!
  //! \param[in]  bucketIndex The bucket index.
  //!
  //! \return     The hash key from bucket index.
  //!
  size_t getHashKeyFromBucketIndex(const Vector3i& bucketIndex) const;

  //!
  //! Gets the bucket index from a point.
  //!
  //! \param[in]  position The position of the point.
  //!
  //! \return     The bucket index.
  //!
  Vector3i getBucketIndex(const Vector3d& position) const;


 private:
  double _gridSpacing = 1.0;
  Vector3i _resolution = Vector3i(64, 64, 64);
  std::vector<Vector3d> _points;
  std::vector<std::vector<size_t>> _buckets;

  size_t getHashKeyFromPosition(const Vector3d& position) const;

  void getNearbyKeys(const Vector3d& position, size_t* bucketIndices) const;
};



