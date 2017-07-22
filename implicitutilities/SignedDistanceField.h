#ifndef SIGNED_DISTANCE_FIELD_H
#define SIGNED_DISTANCE_FIELD_H

#include "MathDefines.h"

class SignedDistanceField final
{

public:

  SignedDistanceField();

  void clear();
  void clearValues();
  void initializeCenteredGrid( const Eigen::Matrix<scalar,6,1>& bounding_box, const Vector3s& cell_width, const Vector3s& grid_padding );

  bool empty() const;
  unsigned nx() const;
  unsigned ny() const;
  unsigned nz() const;
  unsigned numVals() const;
  const scalar& dx() const;
  const scalar& dy() const;
  const scalar& dz() const;
  const scalar& getValue( const unsigned i, const unsigned j, const unsigned k ) const;
  const Vector3u& gridDimensions() const;
  const Vector3s& gridStart() const;
  Vector3s gridEnd() const;
  const Vector3s& cellDelta() const;
  VectorXs& vals();
  const VectorXs& vals() const;

  bool sdfWasComputed() const;

  // True if point x lies inside or on the isosurface of given value
  bool pointIsWithinIsosurface( const Vector3s& x, const scalar& isosurface_value ) const;

  void evaluateGradients( const Matrix3Xsc& points, Matrix3Xsc& gradients ) const;

private:

  Vector3s evaluateGradient( const Vector3s& x ) const;

  // Compute the distance to the zero level set of the signed distance field at the point x using trilinear interpolation
  scalar computeDistanceToSurface( const Vector3s& x ) const;

  // Compute the value of the signed distance field at the point x using trilinear interpolation
  scalar trilinearInterpolation( const scalar& v000, const scalar& v100, const scalar& v010, const scalar& v110,
                                 const scalar& v001, const scalar& v101, const scalar& v011, const scalar& v111,
                                 const Vector3s& bc ) const;

  // Linear interpolation of v0 and v1 parameterized by alpha
  scalar linearInterpolation( const scalar& v0, const scalar& v1, const scalar& alpha ) const;

  void resetValues();
  
  Vector3u m_grid_dimensions;
  Vector3s m_cell_delta;
  Vector3s m_grid_start;
  VectorXs m_sdf_vals;

};

#endif
