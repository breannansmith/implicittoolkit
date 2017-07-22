#include "DistanceTools.h"

#include <Eigen/Geometry>

EdgeIntersectionCountCache::EdgeIntersectionCountCache( const Array3u& dimensions )
: m_dimensions( dimensions )
, m_counts( ArrayXu::Zero( dimensions.x() * dimensions.y() * dimensions.z() ) )
{}

unsigned& EdgeIntersectionCountCache::getCount( const Array3u& edge )
{
  assert( ( edge < m_dimensions ).all() ); assert( m_counts.size() == m_dimensions.x() * m_dimensions.y() * m_dimensions.z() );
  assert( ( edge.z() * m_dimensions.y() + edge.y() ) * m_dimensions.x() + edge.x() < m_counts.size() );
  return m_counts( ( edge.z() * m_dimensions.y() + edge.y() ) * m_dimensions.x() + edge.x() );
}

const unsigned& EdgeIntersectionCountCache::getCount( const Array3u& edge ) const
{
  assert( ( edge < m_dimensions ).all() ); assert( m_counts.size() == m_dimensions.x() * m_dimensions.y() * m_dimensions.z() );
  assert( ( edge.z() * m_dimensions.y() + edge.y() ) * m_dimensions.x() + edge.x() < m_counts.size() );
  return m_counts( ( edge.z() * m_dimensions.y() + edge.y() ) * m_dimensions.x() + edge.x() );
}

bool EdgeIntersectionCountCache::operator==( const EdgeIntersectionCountCache& rhs ) const
{
  return ( m_dimensions == rhs.m_dimensions ).all() && ( m_counts == rhs.m_counts ).all();
}

static void computeSDFSignsGivenIntersectionCount( const Array3u& num_samples, const EdgeIntersectionCountCache& intersection_count, VectorXs& sdf )
{
  assert( num_samples.x() > 0 );
  // For each initial z coord
  for( unsigned k = 0; k < num_samples.z(); ++k )
  {
    // For each initial y coord
    for( unsigned j = 0; j < num_samples.y(); ++j )
    {
      // 'Ray trace' starting from a point outside the mesh to determine the inside/outside status of each grid point
      scalar sign = 1.0;
      // For each segment along the x direction
      for( unsigned i = 0; i < num_samples.x() - 1; ++i )
      {
        // If we intersect the mesh in an odd number of places, flip the sign
        if( intersection_count.getCount( Array3u( i, j, k ) ) % 2 == 1 )
        {
          sign *= -1.0;
        }
        // Set the sign on the distance computation
        sdf( ( k * num_samples.y() + j ) * num_samples.x() + i + 1 ) *= sign;
      }
    }
  }
}

// Given segment pq and triangle abc, returns whether segment intersects triangle
static bool intersectSegmentTriangle( const Vector3s& p0, const Vector3s& q0, const Vector3s& a, const Vector3s& b, const Vector3s& c )
{
  const Vector3s ab = b - a;
  const Vector3s ac = c - a;

  // TODO: Kind of a hack to get this to work but I'm short on time!
  Vector3s p = p0;
  Vector3s q = q0;
  if( ( p0 - q0 ).dot( ab.cross( ac ) ) < 0.0 )
  {
    q = p0;
    p = q0;
  }

  const Vector3s qp = p - q;

  // Compute triangle normal. Can be precalculated or cached if
  // intersecting multiple segments against the same triangle
  const Vector3s n = ab.cross( ac );

  // Compute denominator d. If d <= 0, segment is parallel to or points
  // away from triangle, so exit early
  const scalar d = qp.dot( n );
  if( d == 0.0 ) { return false; }
  assert( d > 0.0 );

  // Compute intersection t value of pq with plane of triangle. A ray
  // intersects iff 0 <= t. Segment intersects iff 0 <= t <= 1. Delay
  // dividing by d until intersection has been found to pierce triangle
  const Vector3s ap = p - a;
  const scalar t = ap.dot( n );
  if( t < 0.0 ) { return false; }
  if( t > d ) { return false; } // For segment; exclude this code line for a ray test

  // Compute barycentric coordinate components and test if within bounds
  const Vector3s e = qp.cross( ap );
  const scalar v = ac.dot( e );
  if( v < 0.0 || v > d ) { return false; }
  const scalar w = -ab.dot( e );
  if( w < 0.0 || v + w > d ) { return false; }

  return true;
}

static void computeEdgeTriangleIntersectionsGivenBoundingIndices( const Array3s& grid_start, const Array3s& cell_width, const Matrix3Xsc& vertices, const Matrix3Xuc& triangles, const unsigned tri_num, const Eigen::Array<unsigned,6,1>& bounding_indices, EdgeIntersectionCountCache& intersection_count )
{
  assert( ( bounding_indices.segment<3>( 0 ) <= bounding_indices.segment<3>( 3 ) ).all() );
  // For each z point in the given range
  for( unsigned z_idx = bounding_indices(2); z_idx <= bounding_indices(5); ++z_idx )
  {
    // For each y point in the given range
    for( unsigned y_idx = bounding_indices(1); y_idx <= bounding_indices(4); ++y_idx )
    {
      // For each x edge in the given range
      for( unsigned x_idx = bounding_indices(0); x_idx < bounding_indices(3); ++x_idx )
      {
        // TODO: Could save some computation here by re-using points
        // Compute the endpoints of the edge
        const Vector3s p0 = grid_start + Array3s( x_idx, y_idx, z_idx ) * cell_width;
        const Vector3s p1 = grid_start + Array3s( x_idx + 1, y_idx, z_idx ) * cell_width;

        // If the edge intersects the current triangle, record that intersection
        const bool intersected = intersectSegmentTriangle( p0, p1, vertices.col( triangles( 0, tri_num ) ), vertices.col( triangles( 1, tri_num ) ), vertices.col( triangles( 2, tri_num ) ) );
        if( intersected )
        {
          intersection_count.getCount( Array3u( x_idx, y_idx, z_idx ) ) += 1;
        }
      }
    }
  }
}

static void computeTriangleBoundingBox( const Matrix3Xsc& vertices, const Matrix3Xuc& triangles, const unsigned tri_num, Eigen::Array<scalar,6,1>& bounding_box )
{
  assert( tri_num < triangles.cols() );
  assert( ( triangles.col( tri_num ).array() < static_cast<unsigned int>(vertices.cols()) ).all() );

  bounding_box.segment<3>( 0 ).setConstant(  std::numeric_limits<scalar>::infinity() );
  bounding_box.segment<3>( 3 ).setConstant( -std::numeric_limits<scalar>::infinity() );

  // For each vertex of the triangle
  for( unsigned vrt_num = 0; vrt_num < 3; ++vrt_num )
  {
    bounding_box.segment<3>( 0 ) = vertices.col( triangles( vrt_num, tri_num ) ).array().min( bounding_box.segment<3>( 0 ) );
    bounding_box.segment<3>( 3 ) = vertices.col( triangles( vrt_num, tri_num ) ).array().max( bounding_box.segment<3>( 3 ) );
  }

  assert( ( bounding_box.segment<3>( 0 ) <= bounding_box.segment<3>( 3 ) ).all() );
}

#ifndef NDEBUG
static void computeEdgeTriangleIntersectionCountWithSpatialSubdivision( const Array3s& grid_start, const Array3u& num_samples, const Array3s& cell_width, const Matrix3Xsc& vertices, const Matrix3Xuc& triangles, EdgeIntersectionCountCache& intersection_count )
#else
static void computeEdgeTriangleIntersectionCountWithSpatialSubdivision( const Array3s& grid_start, const Array3u& /*num_samples*/, const Array3s& cell_width, const Matrix3Xsc& vertices, const Matrix3Xuc& triangles, EdgeIntersectionCountCache& intersection_count )
#endif
{
  // For each triangle
  for( int cur_tri = 0; cur_tri < triangles.cols(); ++cur_tri )
  {
    // Compute a bounding box around the curent triangle
    Eigen::Array<scalar,6,1> bounding_box;
    computeTriangleBoundingBox( vertices, triangles, static_cast<unsigned int>(cur_tri), bounding_box );
    // Inflate the bounding box to account for FPA errors
    bounding_box.segment<3>( 0 ) -= 1.0e-6;
    bounding_box.segment<3>( 3 ) += 1.0e-6;

    // Compute the grid indices that bound the bounding box
    Eigen::Array<unsigned,6,1> bounding_indices;
    bounding_indices.segment<3>( 0 ) = ( ( bounding_box.segment<3>( 0 ) - grid_start ) / cell_width ).unaryExpr( std::ptr_fun( floor ) ).cast<unsigned>();
    bounding_indices.segment<3>( 3 ) = ( ( bounding_box.segment<3>( 3 ) - grid_start ) / cell_width ).unaryExpr( std::ptr_fun( ceil ) ).cast<unsigned>();
    assert( ( bounding_indices.segment<3>( 0 ) <= bounding_indices.segment<3>( 3 ) ).all() );
    assert( ( bounding_indices.segment<3>( 0 ) < num_samples ).all() ); assert( ( bounding_indices.segment<3>( 3 ) < num_samples ).all() );
    computeEdgeTriangleIntersectionsGivenBoundingIndices( grid_start, cell_width, vertices, triangles, static_cast<unsigned int>(cur_tri), bounding_indices, intersection_count );
  }
}

// Computes the signs of a signed distance field of a triangle mesh
static void computeSDFSigns( const Array3s& grid_start, const Array3u& num_samples, const Array3s& cell_width, const Matrix3Xsc& vertices, const Matrix3Xuc& triangles, VectorXs& sdf )
{
  // For each edge in the grid in the x direction, determine how many triangles intersect that edge
  assert( num_samples.x() > 0 );
  EdgeIntersectionCountCache intersection_count( Array3u( num_samples.x() - 1, num_samples.y(), num_samples.z() ) );
  computeEdgeTriangleIntersectionCountWithSpatialSubdivision( grid_start, num_samples, cell_width, vertices, triangles, intersection_count );

  // Using the intersection count, compute the inside/outside sign at each grid point
  computeSDFSignsGivenIntersectionCount( num_samples, intersection_count, sdf );
}

// Given a triangle with vertices a, b, and c, computes the closest point to the input point p
static Vector3s closestPointPointTriangle( const Vector3s& p, const Vector3s& a, const Vector3s& b, const Vector3s& c )
{
  // Check if P in vertex region outside A
  const Vector3s ab = b - a;
  const Vector3s ac = c - a;
  const Vector3s ap = p - a;
  const scalar d1 = ab.dot( ap );
  const scalar d2 = ac.dot( ap );
  if( d1 <= 0.0 && d2 <= 0.0 ) { return a; } // barycentric coordinates (1,0,0)

  // Check if P in vertex region outside B
  const Vector3s bp = p - b;
  const scalar d3 = ab.dot( bp );
  const scalar d4 = ac.dot( bp );
  if( d3 >= 0.0 && d4 <= d3 ) { return b; } // barycentric coordinates (0,1,0)

  // Check if P in edge region of AB, if so return projection of P onto AB
  const scalar vc = d1 * d4 - d3 * d2;
  if( vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 )
  {
    const scalar v = d1 / ( d1 - d3 );
    return a + v * ab; // barycentric coordinates (1-v,v,0)
  }

  // Check if P in vertex region outside C
  const Vector3s cp = p - c;
  const scalar d5 = ab.dot( cp );
  const scalar d6 = ac.dot( cp );
  if( d6 >= 0.0 && d5 <= d6 ) { return c; } // barycentric coordinates (0,0,1)

  // Check if P in edge region of AC, if so return projection of P onto AC
  const scalar vb = d5 * d2 - d1 * d6;
  if( vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 )
  {
    const scalar w = d2 / ( d2 - d6 );
    return a + w * ac; // barycentric coordinates (1-w,0,w)
  }

  // Check if P in edge region of BC, if so return projection of P onto BC
  const scalar va = d3 * d6 - d5 * d4;
  if( va <= 0.0 && ( d4 - d3 ) >= 0.0 && ( d5 - d6 ) >= 0.0 )
  {
    const scalar w = ( d4 - d3 ) / ( ( d4 - d3 ) + ( d5 - d6 ) );
    return b + w * ( c - b ); // barycentric coordinates (0,1-w,w)
  }

  // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
  const scalar denom = 1.0 / ( va + vb + vc );
  const scalar v = vb * denom;
  const scalar w = vc * denom;
  return a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0f - v - w
}

// Given a triangle with vertices a, b, and c, computes the squared distance to point p
static scalar squaredDistancePointTriangle( const Vector3s& p, const Vector3s& a, const Vector3s& b, const Vector3s& c )
{
  const Vector3s closest_point = closestPointPointTriangle( p, a, b, c );
  return ( p - closest_point ).squaredNorm();
}

// Computes the distance between point p and a mesh using a brute force comparison
static scalar distancePointMeshBruteForce( const Vector3s& p, const Matrix3Xsc& vertices, const Matrix3Xuc& triangles )
{
  scalar min_squared_dist = std::numeric_limits<scalar>::infinity();
  for( int tri_num = 0; tri_num < triangles.cols(); ++tri_num )
  {
    assert( ( triangles.col( tri_num ).array() < static_cast<unsigned int>(vertices.cols()) ).all() );
    const scalar cur_squared_dist = squaredDistancePointTriangle( p, vertices.col( triangles( 0, tri_num ) ), vertices.col( triangles( 1, tri_num ) ), vertices.col( triangles( 2, tri_num ) ) );
    assert( cur_squared_dist >= 0.0 );
    min_squared_dist = std::min( cur_squared_dist, min_squared_dist );
  }
  return sqrt( min_squared_dist );
}

void DistanceTools::computeSignedDistanceField( const Vector3s& grid_start, const Array3u& num_samples, const Vector3s& step_width, const Matrix3Xsc& vertices, const Matrix3Xuc& triangles, VectorXs& sdf, AlgorithmProgressTracker* tracker )
{
  assert( sdf.size() == num_samples.x() * num_samples.y() * num_samples.z() );
  assert( ( num_samples > 0 ).all() );

  int next_percent = 1;
  // For each grid point
  for( unsigned k = 0; k < num_samples.z(); ++k ) for( unsigned j = 0; j < num_samples.y(); ++j ) for( unsigned i = 0; i < num_samples.x(); ++i )
  {
    // Compute the position of the grid point
    const Vector3s current_point = grid_start.array() + Array3s( i, j, k ) * step_width.array();
    // Compute the distance to the grid point
    const scalar dist = distancePointMeshBruteForce( current_point, vertices, triangles );
    assert( dist >= 0.0 );
    // Save the grid position
    const int sample_num = static_cast<int>(( k * num_samples.y() + j ) * num_samples.x() + i);
    assert( sample_num < sdf.size() );
    sdf( sample_num ) = dist;

    if( tracker != nullptr )
    {
      if( tracker->algorithmCancelled() ) return;
      const int percent = int( 100.0 * scalar( sample_num ) / scalar( sdf.size() ) );
      if( percent >= next_percent )
      {
        tracker->setProgressIndicator( percent );
        next_percent = percent + 1;
      }
    }
  }

  // Compute the sign of each grid point
  computeSDFSigns( grid_start, num_samples, step_width, vertices, triangles, sdf );
}
