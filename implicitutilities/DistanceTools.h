#ifndef DISTANCE_TOOLS_H
#define DISTANCE_TOOLS_H

#include "MathDefines.h"

class AlgorithmProgressTracker
{

public:

  virtual ~AlgorithmProgressTracker() = default;
  AlgorithmProgressTracker(const AlgorithmProgressTracker&) = delete;
  AlgorithmProgressTracker(AlgorithmProgressTracker&&) = delete;
  AlgorithmProgressTracker& operator=(const AlgorithmProgressTracker&) = delete;
  AlgorithmProgressTracker& operator=(AlgorithmProgressTracker&&) = delete;

  // Value should be in the range 0...100 inclusive
  virtual void setProgressIndicator( const int value ) = 0;

  // True if the user wants to terminate this run
  virtual bool algorithmCancelled() = 0;

protected:

  AlgorithmProgressTracker() = default;

};

class EdgeIntersectionCountCache final
{

public:

  EdgeIntersectionCountCache( const Array3u& dimensions );

  unsigned& getCount( const Array3u& edge );
  const unsigned& getCount( const Array3u& edge ) const;

  bool operator==( const EdgeIntersectionCountCache& rhs ) const;

private:

  const Array3u m_dimensions;
  ArrayXu m_counts;

};

namespace DistanceTools
{

  // TODO: Pass in signed distance field instead of lots of things
  void computeSignedDistanceField( const Vector3s& grid_start, const Array3u& num_samples, const Vector3s& step_width, const Matrix3Xsc& vertices, const Matrix3Xuc& triangles, VectorXs& sdf, AlgorithmProgressTracker* tracker = nullptr );

}

#endif
