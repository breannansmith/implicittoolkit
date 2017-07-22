#ifndef SURFACE_SAMPLING_H
#define SURFACE_SAMPLING_H

#include "MathDefines.h"

class SurfaceSampling final
{

public:

  SurfaceSampling();

  void computeSurfaceSampling( const Matrix3Xsc& vertices, const Matrix2Xuc& edges, const Matrix3Xuc& faces );

  const Matrix3Xsc& samples() const;

  bool empty() const;

  unsigned numSamples() const;
  unsigned numEdgeSubSamples() const;
  unsigned numFaceSubSamples() const;

  void setEdgeSubsamples( const unsigned num_samples );
  void setFaceSubsamples( const unsigned num_samples );

private:

  unsigned m_edge_subsamples;
  unsigned m_face_subsamples;
  Matrix3Xsc m_samples;

};

#endif
