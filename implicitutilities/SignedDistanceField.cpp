#include "SignedDistanceField.h"

#include <limits>

SignedDistanceField::SignedDistanceField()
: m_grid_dimensions( Vector3u::Zero() )
, m_cell_delta( Vector3s::Constant( std::numeric_limits<scalar>::signaling_NaN() ) )
, m_grid_start( Vector3s::Constant( std::numeric_limits<scalar>::signaling_NaN() ) )
, m_sdf_vals()
{}

void SignedDistanceField::clear()
{
  m_grid_dimensions.setConstant( 0 );
  m_cell_delta.setConstant( std::numeric_limits<scalar>::signaling_NaN() );
  m_grid_start.setConstant( std::numeric_limits<scalar>::signaling_NaN() );
  resetValues();
}

void SignedDistanceField::clearValues()
{
  resetValues();
}

void SignedDistanceField::initializeCenteredGrid( const Eigen::Matrix<scalar,6,1>& bounding_box, const Vector3s& cell_width, const Vector3s& grid_padding )
{
  const Vector3s grid_center = 0.5 * ( bounding_box.segment<3>( 0 ) + bounding_box.segment<3>( 3 ) );

  assert( ( grid_padding.array() >= 0.0 ).all() );
  Eigen::Matrix<scalar,6,1> inflated_bounding_box = bounding_box;
  inflated_bounding_box.segment<3>( 0 ) -= grid_padding;
  inflated_bounding_box.segment<3>( 3 ) += grid_padding;

  assert( ( cell_width.array() > 0.0 ).all() );
  m_grid_dimensions = ( ( inflated_bounding_box.segment<3>( 3 ) - inflated_bounding_box.segment<3>( 0 ) ).array() / cell_width.array() ).unaryExpr([](scalar s){return ceil(s);}).cast<unsigned>() + 1;;
  assert( ( m_grid_dimensions.array() >= 2 ).all() );

  // Compute the width of the grid
  const Vector3s grid_width = cell_width.array() * ( m_grid_dimensions.array() - 1 ).cast<scalar>().array();

  // Compute the first point of the grid
  m_grid_start = grid_center - 0.5 * grid_width;

  m_cell_delta = cell_width;

  resetValues();
}

bool SignedDistanceField::empty() const
{
  // return m_sdf_vals( 0 ) != m_sdf_vals( 0 );
  return (m_grid_dimensions.array() == 0).any();
}

unsigned SignedDistanceField::nx() const
{
  return m_grid_dimensions.x();
}

unsigned SignedDistanceField::ny() const
{
  return m_grid_dimensions.y();
}

unsigned SignedDistanceField::nz() const
{
  return m_grid_dimensions.z();
}

unsigned SignedDistanceField::numVals() const
{
  return nx() * ny() * nz();
}

const scalar& SignedDistanceField::dx() const
{
  return m_cell_delta.x();
}

const scalar& SignedDistanceField::dy() const
{
  return m_cell_delta.y();
}

const scalar& SignedDistanceField::dz() const
{
  return m_cell_delta.z();
}

const scalar& SignedDistanceField::getValue( const unsigned i, const unsigned j, const unsigned k ) const
{
  assert( i < nx() );
  assert( j < ny() );
  assert( k < nz() );
  assert( ( k * ny() + j ) * nx() + i < m_sdf_vals.size() );
  
  return m_sdf_vals( ( k * ny() + j ) * nx() + i );
}

const Vector3u& SignedDistanceField::gridDimensions() const
{
  return m_grid_dimensions;
}

const Vector3s& SignedDistanceField::gridStart() const
{
  return m_grid_start;
}

Vector3s SignedDistanceField::gridEnd() const
{
  assert( ( m_grid_dimensions.array() >= 1 ).all() );
  return m_grid_start + ( ( m_grid_dimensions.array() - 1 ).cast<scalar>() * m_cell_delta.array() ).matrix();
}

const Vector3s& SignedDistanceField::cellDelta() const
{
  return m_cell_delta;
}

VectorXs& SignedDistanceField::vals()
{
  return m_sdf_vals;
}

const VectorXs& SignedDistanceField::vals() const
{
  return m_sdf_vals;
}

void SignedDistanceField::resetValues()
{
  m_sdf_vals = VectorXs::Constant( numVals(), std::numeric_limits<scalar>::signaling_NaN() );
}

bool SignedDistanceField::sdfWasComputed() const
{
  return numVals() != 0 && !std::isnan( m_sdf_vals( 0 ) );
}

void SignedDistanceField::evaluateGradients( const Matrix3Xsc& points, Matrix3Xsc& gradients ) const
{
  if( empty() )
  {
    return;
  }

  gradients.resize( 3, points.cols() );
  for( int pnt_num = 0; pnt_num < points.cols(); ++pnt_num )
  {
    gradients.col( pnt_num ) = evaluateGradient( points.col( pnt_num ) );
  }
}

// TODO: Can probably avoid a bunch of dx/dy divisions here
Vector3s SignedDistanceField::evaluateGradient( const Vector3s& x ) const
{
  assert( ( x.array() >= m_grid_start.array() ).all() );
  assert( ( x.array() <= gridEnd().array() ).all() );

  // Determine which cell this piont lies within
  const Array3u indices = ( ( x - m_grid_start ).array() / m_cell_delta.array() ).unaryExpr( [](scalar s){return floor(s);} ).cast<unsigned>();
  assert( ( indices + 1 < m_grid_dimensions.array() ).all() );

  // Compute the 'barycentric' (and one minus) coordinates of the point in the cell
  const Vector3s bc = ( x.array() - ( m_grid_start.array() + indices.cast<scalar>().array() * m_cell_delta.array() ) ) / m_cell_delta.array();
  assert( ( bc.array() >= -1.0e-12 ).all() ); assert( ( bc.array() <= 1.0 + 1.0e-12 ).all() );
  const Vector3s bci = Vector3s::Ones() - bc;

  // Grab the value of the distance field at each grid point
  const scalar v000 = getValue( indices.x(),     indices.y(),     indices.z() );
  const scalar v100 = getValue( indices.x() + 1, indices.y(),     indices.z() );
  const scalar v010 = getValue( indices.x(),     indices.y() + 1, indices.z() );
  const scalar v110 = getValue( indices.x() + 1, indices.y() + 1, indices.z() );
  const scalar v001 = getValue( indices.x(),     indices.y(),     indices.z() + 1 );
  const scalar v101 = getValue( indices.x() + 1, indices.y(),     indices.z() + 1 );
  const scalar v011 = getValue( indices.x(),     indices.y() + 1, indices.z() + 1 );
  const scalar v111 = getValue( indices.x() + 1, indices.y() + 1, indices.z() + 1 );

  Vector3s gradient;

  gradient.x() = bci.z() * ( bci.y() * ( v100 - v000 ) + bc.y() * ( v110 - v010 ) )
                + bc.z() * ( bci.y() * ( v101 - v001 ) + bc.y() * ( v111 - v011 ) );

  gradient.y() = bci.z() * ( bci.x() * ( v010 - v000 ) + bc.x() * ( v110 - v100 ) )
                + bc.z() * ( bci.x() * ( v011 - v001 ) + bc.x() * ( v111 - v101 ) );

  gradient.z() = bci.y() * ( bci.x() * ( v001 - v000 ) + bc.x() * ( v101 - v100 ) )
                + bc.y() * ( bci.x() * ( v011 - v010 ) + bc.x() * ( v111 - v110 ) );

  gradient.array() /= m_cell_delta.array();
  gradient.normalize();

  return gradient;
}

bool SignedDistanceField::pointIsWithinIsosurface( const Vector3s& x, const scalar& isosurface_value ) const
{
  // If the point lies outside the grid, reject it
  if( ( x.array() < m_grid_start.array() ).any() )
  {
    return false;
  }
  if( ( x.array() > m_grid_start.array() + ( m_grid_dimensions.array() - 1 ).cast<scalar>() * m_cell_delta.array() ).any() )
  {
    return false;
  }
  
  // Compute the distance to the zero level set
  const scalar dist = computeDistanceToSurface( x );

  if( dist <= isosurface_value )
  {
    return true;
  }
  
  return false;
}

scalar SignedDistanceField::computeDistanceToSurface( const Vector3s& x ) const
{
  // Determine which cell this point lies in
  const Vector3u indices = ( ( x - m_grid_start ).array() / m_cell_delta.array() ).unaryExpr( [](scalar s){return floor(s);} ).cast<unsigned>();
  assert( ( indices.array() < m_grid_dimensions.array() ).all() );

  // Handle points on the 'right' boundary of the grid
  if( indices.x() + 1 == unsigned(nx()) )
  {
    return std::numeric_limits<scalar>::infinity();
  }
  if( indices.y() + 1 == unsigned(ny()) )
  {
    return std::numeric_limits<scalar>::infinity();
  }
  if( indices.z() + 1 == unsigned(nz()) )
  {
    return std::numeric_limits<scalar>::infinity();
  }

  // Compute the 'barycentric' coordinates of the point in the cell
  const Vector3s bc = ( x.array() - ( m_grid_start.array() + indices.cast<scalar>().array() * m_cell_delta.array() ) ) / m_cell_delta.array();
  assert( ( bc.array() >= 0.0 ).all() ); assert( ( bc.array() <= 1.0 ).all() );
  
  // Grab the value of the distance field at each grid point
  const scalar v000 = getValue( indices.x(),     indices.y(),     indices.z() );
  const scalar v100 = getValue( indices.x() + 1, indices.y(),     indices.z() );
  const scalar v010 = getValue( indices.x(),     indices.y() + 1, indices.z() );
  const scalar v110 = getValue( indices.x() + 1, indices.y() + 1, indices.z() );
  const scalar v001 = getValue( indices.x(),     indices.y(),     indices.z() + 1 );
  const scalar v101 = getValue( indices.x() + 1, indices.y(),     indices.z() + 1 );
  const scalar v011 = getValue( indices.x(),     indices.y() + 1, indices.z() + 1 );
  const scalar v111 = getValue( indices.x() + 1, indices.y() + 1, indices.z() + 1 );
  
  // Compute the value of the distance field at the sample point with trilinear interpolation
  return trilinearInterpolation( v000, v100, v010, v110, v001, v101, v011, v111, bc );
}

scalar SignedDistanceField::trilinearInterpolation( const scalar& v000, const scalar& v100, const scalar& v010, const scalar& v110,
                                                    const scalar& v001, const scalar& v101, const scalar& v011, const scalar& v111,
                                                    const Vector3s& bc ) const
{
  assert( ( bc.array() >= 0.0 ).all() );
  assert( ( bc.array() <= 1.0 ).all() );
  
  const scalar v00 = linearInterpolation( v000, v100, bc.x() );
  const scalar v10 = linearInterpolation( v010, v110, bc.x() );
  const scalar v01 = linearInterpolation( v001, v101, bc.x() );
  const scalar v11 = linearInterpolation( v011, v111, bc.x() );
  
  const scalar v0 = linearInterpolation( v00, v10, bc.y() );
  const scalar v1 = linearInterpolation( v01, v11, bc.y() );
  
  const scalar v = linearInterpolation( v0, v1, bc.z() );
  
  return v;
}

scalar SignedDistanceField::linearInterpolation( const scalar& v0, const scalar& v1, const scalar& alpha ) const
{
  assert( alpha >= 0.0 ); assert( alpha <= 1.0 );
  return ( 1.0 - alpha ) * v0 + alpha * v1;
}
