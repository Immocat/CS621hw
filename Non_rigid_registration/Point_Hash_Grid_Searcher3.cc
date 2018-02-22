#include "Point_Hash_Grid_Searcher3.hh"


bool PointHashGridSearcher3::hasNearbyPoint(const Vector3d& origin,
                                            double radius) const {
  if (_buckets.empty()) {
    return false;
  }

  size_t nearbyKeys[8];
  getNearbyKeys(origin, nearbyKeys);

  const double queryRadiusSquared = radius * radius;

  for (int i = 0; i < 8; i++) {
    const auto& bucket = _buckets[nearbyKeys[i]];
    size_t numberOfPointsInBucket = bucket.size();

    for (size_t j = 0; j < numberOfPointsInBucket; ++j) {
      size_t pointIndex = bucket[j];
      double rSquared = length2(_points[pointIndex] - origin);
      if (rSquared <= queryRadiusSquared) {
        return true;
      }
    }
  }

  return false;
}

void PointHashGridSearcher3::add(const Vector3d& point) {
  if (_buckets.empty()) {
    std::vector<Vector3d> arr = {point};
    build(arr);
  } else {
    size_t i = _points.size();
    _points.push_back(point);
    size_t key = getHashKeyFromPosition(point);
    _buckets[key].push_back(i);
  }
}
void PointHashGridSearcher3::build(
    const std::vector<Vector3d>& points) {
  _buckets.clear();
  _points.clear();

  // Allocate memory chuncks
  _buckets.resize(_resolution[0] * _resolution[1] * _resolution[2]);
  _points.resize(points.size());

  if (points.size() == 0) {
    return;
  }

  // Put points into buckets
  for (size_t i = 0; i < points.size(); ++i) {
    _points[i] = points[i];
    size_t key = getHashKeyFromPosition(points[i]);
    _buckets[key].push_back(i);
  }
}
const std::vector<std::vector<size_t>>& PointHashGridSearcher3::buckets()
    const {
  return _buckets;
}

Vector3i PointHashGridSearcher3::getBucketIndex(const Vector3d& position) const {
  Vector3i bucketIndex;
  bucketIndex[0] = static_cast<ssize_t>(std::floor(position[0] / _gridSpacing));
  bucketIndex[1] = static_cast<ssize_t>(std::floor(position[1] / _gridSpacing));
  bucketIndex[2] = static_cast<ssize_t>(std::floor(position[2] / _gridSpacing));
  return bucketIndex;
}

size_t PointHashGridSearcher3::getHashKeyFromPosition(
    const Vector3d& position) const {
  Vector3i bucketIndex = getBucketIndex(position);
  return getHashKeyFromBucketIndex(bucketIndex);
}

size_t PointHashGridSearcher3::getHashKeyFromBucketIndex(
    const Vector3i& bucketIndex) const {
  Vector3i wrappedIndex = bucketIndex;
  wrappedIndex[0] = bucketIndex[0] % _resolution[0];
  wrappedIndex[1] = bucketIndex[1] % _resolution[1];
  wrappedIndex[2] = bucketIndex[2] % _resolution[2];
  if (wrappedIndex[0] < 0) {
    wrappedIndex[0] += _resolution[0];
  }
  if (wrappedIndex[1] < 0) {
    wrappedIndex[1] += _resolution[1];
  }
  if (wrappedIndex[2] < 0) {
    wrappedIndex[2] += _resolution[2];
  }
  return static_cast<size_t>((wrappedIndex[2] * _resolution[1] + wrappedIndex[1]) *
                                 _resolution[0] +
                             wrappedIndex[0]);
}

void PointHashGridSearcher3::getNearbyKeys(const Vector3d& position,
                                           size_t* nearbyKeys) const {
  Vector3i originIndex = getBucketIndex(position), nearbyBucketIndices[8];

  for (int i = 0; i < 8; i++) {
    nearbyBucketIndices[i] = originIndex;
  }

  if ((originIndex[0] + 0.5f) * _gridSpacing <= position[0]) {
    nearbyBucketIndices[4][0] += 1;
    nearbyBucketIndices[5][0] += 1;
    nearbyBucketIndices[6][0] += 1;
    nearbyBucketIndices[7][0] += 1;
  } else {
    nearbyBucketIndices[4][0] -= 1;
    nearbyBucketIndices[5][0] -= 1;
    nearbyBucketIndices[6][0] -= 1;
    nearbyBucketIndices[7][0] -= 1;
  }

  if ((originIndex[1] + 0.5f) * _gridSpacing <= position[1]) {
    nearbyBucketIndices[2][1] += 1;
    nearbyBucketIndices[3][1] += 1;
    nearbyBucketIndices[6][1] += 1;
    nearbyBucketIndices[7][1] += 1;
  } else {
    nearbyBucketIndices[2][1] -= 1;
    nearbyBucketIndices[3][1] -= 1;
    nearbyBucketIndices[6][1] -= 1;
    nearbyBucketIndices[7][1] -= 1;
  }

  if ((originIndex[2] + 0.5f) * _gridSpacing <= position[2]) {
    nearbyBucketIndices[1][2] += 1;
    nearbyBucketIndices[3][2] += 1;
    nearbyBucketIndices[5][2] += 1;
    nearbyBucketIndices[7][2] += 1;
  } else {
    nearbyBucketIndices[1][2] -= 1;
    nearbyBucketIndices[3][2] -= 1;
    nearbyBucketIndices[5][2] -= 1;
    nearbyBucketIndices[7][2] -= 1;
  }

  for (int i = 0; i < 8; i++) {
    nearbyKeys[i] = getHashKeyFromBucketIndex(nearbyBucketIndices[i]);
  }
}

