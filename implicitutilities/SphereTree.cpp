#include "SphereTree.h"

#include <Eigen/Eigenvalues>
#include <iostream>
#include <numeric>
#include <cmath>
#include <random>
#include <algorithm>

SphereTree::SphereTree()
: m_tree_depth( 0 )
, m_x()
, m_r()
, m_leaves()
{}

unsigned SphereTree::depth() const
{
  return m_tree_depth;
}

const Matrix3Xsc& SphereTree::nodeCenters() const
{
  return m_x;
}

const VectorXs& SphereTree::nodeRadii() const
{
  return m_r;
}

void SphereTree::extractLeaves( VectorXu& node_ids, Matrix2Xuc& contents ) const
{
  const unsigned num_leaves = unsigned(m_leaves.size());
  node_ids.resize( num_leaves );
  contents.resize( 2, num_leaves );
  unsigned i = 0;
  for( const std::pair<const unsigned,Vector2u>& leaf : m_leaves )
  {
    node_ids( i ) = leaf.first;
    contents.col( i++ ) = leaf.second;
  }
}

static unsigned pow2( const unsigned n )
{
  int total = 1;
  for( unsigned j = 0; j < n; ++j )
  {
    total *= 2;
  }
  return unsigned(total);
}

void SphereTree::spheres( const unsigned depth, std::vector<Vector3s>& centers, std::vector<scalar>& radii ) const
{
  if( depth > m_tree_depth )
  {
    return;
  }
  // Starting index of this level in storage
  const unsigned start_idx = pow2( depth ) - 1;
  // Number of nodes in this level
  const unsigned num_nodes = start_idx + 1;
  // Grab the sphere centers and radii
  for( unsigned node_idx = start_idx; node_idx < start_idx + num_nodes; ++node_idx )
  {
    assert( node_idx < m_x.cols() ); assert( node_idx < m_r.size() );
    if( !std::isnan( m_r( node_idx ) ) )
    {
      assert( !std::isnan( m_x( 0, node_idx ) ) ); assert( !std::isnan( m_x( 1, node_idx ) ) ); assert( !std::isnan( m_x( 2, node_idx ) ) );
      centers.emplace_back( m_x.col( node_idx ) );
      radii.emplace_back( m_r( node_idx ) );
    }
    // Should only see nan on the lowest level
    #ifndef NDEBUG
    else
    {
      assert( depth == m_tree_depth );
    }
    #endif
  }
}

static unsigned treeDepth( const unsigned x )
{
  if( x <= 2 )
  {
    return 0;
  }
  return 1 + std::max( treeDepth( x / 2 + x % 2 ), treeDepth( x / 2 ) );
}

static unsigned numTreeNodes( const unsigned n )
{
  int total = 1;
  for( unsigned j = 0; j <= n; ++j )
  {
    total *= 2;
  }
  return unsigned(total - 1);
}

void SphereTree::buildTree( const Matrix3Xsc& vertices, const Vector3s& grid_start, const Vector3s& grid_end )
{
  const unsigned num_points = unsigned(vertices.cols());
  m_tree_depth = treeDepth( num_points );
  const unsigned num_tree_nodes = numTreeNodes( m_tree_depth );

  // Initialize space for the tree
  m_x = Matrix3Xsc::Constant( 3, num_tree_nodes, std::numeric_limits<scalar>::signaling_NaN() );
  m_r = VectorXs::Constant( num_tree_nodes, std::numeric_limits<scalar>::signaling_NaN() );
  m_leaves.clear();

  // Root of the sphere tree bounds the grid
  m_x.col( 0 ) = 0.5 * ( grid_end + grid_start );
  m_r( 0 ) = ( grid_end - m_x.col( 0 ) ).norm();
  assert( fabs( m_r( 0 ) - ( m_x.col( 0 ) - grid_start ).norm() ) <= 1.0e-6 );

  // If there are no points, we are done
  if( num_points == 0 )
  {
    return;
  }

  // Otherwise build the tree
  std::vector<unsigned> indices( num_points );
  std::iota( std::begin( indices ), std::end( indices ), 0 );
  splitNode( vertices, indices, 0 );

  #ifndef NDEBUG
  // Ensure that all points are accounted for
  {
    std::vector<bool> index_found( static_cast<unsigned long>(vertices.cols()), false );
    for( const std::pair<const unsigned,Vector2u>& entry : m_leaves )
    {
      index_found[ entry.second(0) ] = true;
      index_found[ entry.second(1) ] = true;
    }
    assert( std::all_of( std::begin(index_found), std::end(index_found), [](bool b){ return b; } ) );
  }
  // Ensure that points are not duplicated across leaves
  {
    int sum = 0;
    for( const std::pair<const unsigned,Vector2u>& entry : m_leaves )
    {
      sum += int(entry.second(0));
      if( entry.second(0) != entry.second(1) )
      {
        sum += int(entry.second(1));
      }
    }
    const int computed = int(( vertices.cols() * vertices.cols() - vertices.cols() ) / 2);
    assert( sum == computed );
  }
  #endif
}

static Vector3s averageVertices( const Matrix3Xsc& vertices, const std::vector<unsigned>& indices )
{
  Vector3s average = Vector3s::Zero();
  for( const unsigned index : indices )
  {
    assert( index < vertices.cols() );
    average += vertices.col( index );
  }
  return average / scalar( indices.size() );
}

static Matrix33s computeCovarianceMatrix( const Matrix3Xsc& vertices, const std::vector<unsigned>& indices, const Vector3s& average )
{
  Matrix33s covariance = Matrix33s::Zero();
  for( const unsigned index : indices )
  {
    assert( index < vertices.cols() );
    covariance += ( vertices.col( index ) - average ) * ( vertices.col( index ) - average ).transpose();
  }
  assert( ( covariance - covariance.transpose() ).lpNorm<Eigen::Infinity>() <= 1.0e-6 );
  return covariance / scalar( indices.size() );
}

// TODO: Sort the eigenvalues by magnitude
static Vector3s computeMaximumVarianceAxis( const Matrix33s& C )
{
  assert( ( C - C.transpose() ).lpNorm<Eigen::Infinity>() <= 1.0e-6 );

  // Compute the eigenvectors and eigenvalues of the input matrix
  const Eigen::SelfAdjointEigenSolver<Matrix33s> es( C );

  // Check for errors
  if( es.info() == Eigen::NumericalIssue )
  {
    std::cerr << "Warning, failed to compute eigenvalues of covariance matrix due to Eigen::NumericalIssue" << std::endl;
  }
  else if( es.info() == Eigen::NoConvergence )
  {
    std::cerr << "Warning, failed to compute eigenvalues of covariance matrix due to Eigen::NoConvergence" << std::endl;
  }
  else if( es.info() == Eigen::InvalidInput )
  {
    std::cerr << "Warning, failed to compute eigenvalues of covariance matrix due to Eigen::InvalidInput" << std::endl;
  }
  assert( es.info() == Eigen::Success );

  // Save the eigenvectors and eigenvalues
  #ifndef NDEBUG
  {
    const Vector3s eigs = es.eigenvalues();
    assert( ( eigs.array() >= -1.0e-15 ).all() );
    assert( eigs.x() <= eigs.y() );
    assert( eigs.y() <= eigs.z() );
  }
  #endif

  assert( fabs( es.eigenvectors().col(2).norm() - 1.0 ) <= 1.0e-6 );
  return es.eigenvectors().col(2).normalized();
}

static Vector3s computeMaximumSpreadAxis( const Matrix3Xsc& vertices, const std::vector<unsigned>& indices )
{
  // Compute the average of all current vertices
  const Vector3s average = averageVertices( vertices, indices );

  // Compute the covariance of all current vertices
  const Matrix33s covariance = computeCovarianceMatrix( vertices, indices, average );

  // Compute the axis along which points are maximally spread
  return computeMaximumVarianceAxis( covariance );
}

static void projectOntoAxis( const Matrix3Xsc& vertices, const std::vector<unsigned>& indices, const Vector3s& axis, std::vector<scalar>& projections )
{
  projections.reserve( indices.size() );
  for( const unsigned index : indices )
  {
    assert( index < vertices.cols() );
    projections.emplace_back( axis.dot( vertices.col( index ) ) );
  }
}

static void expandSphereToContainPoint( Vector3s& x, scalar& r, const Vector3s& point )
{
  const Vector3s d = point - x;
  scalar dist = d.squaredNorm();
  if( dist > r * r )
  {
    dist = sqrt( dist );
    const scalar new_radius = 0.5 * ( r + dist );
    const scalar k = ( new_radius - r ) / dist;
    r = new_radius;
    x += d * k;
  }
}

static void computeBoundingSphere( const Matrix3Xsc& vertices, const std::vector<unsigned>& indices, const std::vector<scalar>& projections_on_axis, Vector3s& center, scalar& radius )
{
  if( indices.size() == 1 )
  {
    center = vertices.col( indices.front() );
    radius = 0.0;
    return;
  }

  auto min_max = std::minmax_element( std::begin( projections_on_axis ), std::end( projections_on_axis ) );
  const unsigned min_idx = unsigned(min_max.first - std::begin( projections_on_axis ));
  const unsigned max_idx = unsigned(min_max.second - std::begin( projections_on_axis ));
  assert( projections_on_axis[min_idx] <= projections_on_axis[max_idx] );
  center = 0.5 * ( vertices.col( indices[min_idx] ) + vertices.col( indices[max_idx] ) );
  radius = 0.5 * ( vertices.col( indices[max_idx] ) - vertices.col( indices[min_idx] ) ).norm();
  for( unsigned sample_num = 0; sample_num < indices.size(); ++sample_num )
  {
    expandSphereToContainPoint( center, radius, vertices.col( indices[ sample_num ] ) );
  }

  // Iterative technique to reduce the sphere size. Doesn't help much, but doesn't hurt, either
  {
    std::vector<unsigned> indices_copy = indices;

    std::mt19937_64 mt( 0 );

    Vector3s center2 = center;
    scalar radius2 = radius;
    for( int k = 0; k < 8; ++k )
    {
      // Shrink sphere somewhat to make it an underestimate (not bound)
      radius2 = 0.95 * radius2;

      // Make sphere bound data again
      for( unsigned i = 0; i < indices_copy.size() - 1; ++i )
      {
        // Swap pt[i] with pt[j], where j randomly from interval [i+1,numPts-1]
        std::uniform_int_distribution<unsigned> swap_gen( i + 1, unsigned(indices_copy.size() - 1) );
        const unsigned j = swap_gen( mt );
        assert( i < indices_copy.size() ); assert( j < indices_copy.size() );
        std::swap( indices_copy[i], indices_copy[j] );
        expandSphereToContainPoint( center2, radius2, vertices.col( indices_copy[i] ) );
      }
      expandSphereToContainPoint( center2, radius2, vertices.col( indices_copy[indices_copy.size()-1] ) );

      // Update s whenever a tighter sphere is found
      if( radius2 < radius )
      {
        radius = radius2;
        center = center2;
      }
    }
  }

  radius += 1.0e-6; // To account for FPA errors
}

// TODO: Clean partitionPoints up
static void partitionPoints( const std::vector<unsigned>& indices, std::vector<scalar>& projections, std::vector<unsigned>& left, std::vector<unsigned>& right )
{
  std::vector<std::pair<unsigned,unsigned>> indices_to_use;
  for( unsigned i = 0; i < indices.size(); ++i )
  {
    indices_to_use.emplace_back( std::make_pair( i, indices[i] ) );
  }
  std::vector<std::pair<unsigned,unsigned>>::iterator new_first = indices_to_use.begin();
  std::vector<std::pair<unsigned,unsigned>>::iterator new_last = indices_to_use.end();
  std::vector<std::pair<unsigned,unsigned>>::iterator new_middle = new_first + (new_last - new_first) / 2;
  auto lambda = [&projections]( const std::pair<unsigned,unsigned>& i1, const std::pair<unsigned,unsigned>& i2 )
  {
    assert( i1.first < projections.size() );
    assert( i2.first < projections.size() );
    return projections[i1.first] < projections[i2.first];
  };
  std::nth_element( new_first, new_middle, new_last, lambda );
  if( indices.size() % 2 == 1 ) new_middle++;
  for( auto itr = new_first; itr != new_middle; ++itr )
  {
    left.emplace_back( itr->second );
  }
  for( auto itr = new_middle; itr != new_last; ++itr )
  {
    right.emplace_back( itr->second );
  }
}

void SphereTree::splitNode( const Matrix3Xsc& vertices, std::vector<unsigned>& indices, const unsigned node_index )
{
  // Must be a valid node
  assert( node_index < m_x.cols() );
  // Empty case shouldn't be reached
  assert( !indices.empty() );

  std::vector<unsigned> left;
  std::vector<unsigned> right;

  {
    // Compute the axis along which points are maximally spread
    const Vector3s max_spread_axis = computeMaximumSpreadAxis( vertices, indices );

    // Project points onto the axis of maximum spread
    std::vector<scalar> projections_on_axis;
    projectOntoAxis( vertices, indices, max_spread_axis, projections_on_axis );

    // Create a bounding sphere for the current points
    if( node_index != 0 )
    {
      Vector3s center;
      scalar radius;
      computeBoundingSphere( vertices, indices, projections_on_axis, center, radius );
      m_x.col( node_index ) = center;
      m_r( node_index ) = radius;
    }

    // All input points should be within the current bounding sphere
    #ifndef NDEBUG
    for( unsigned sample_num = 0; sample_num < indices.size(); ++sample_num )
    {
      assert( ( vertices.col( indices[ sample_num ] ) - m_x.col( node_index ) ).squaredNorm() <= m_r( node_index ) * m_r( node_index ) );
    }
    #endif

    // If there are 1 or 2 points, we have hit a leaf node
    if( indices.size() <= 2 )
    {
      // TODO: Check return value here
      m_leaves.insert( std::make_pair( node_index, Vector2u( indices[ 0 ], ( indices.size() == 1 ) ? indices[ 0 ] : indices[ 1 ] ) ) );
      return;
    }

    // Partition the points into two sets along the axis of maximum spread
    partitionPoints( indices, projections_on_axis, left, right );

    assert( left.size() + right.size() == indices.size() );
    assert( left.size() == right.size() || left.size() - 1 == right.size() );
    #ifndef NDEBUG
    // Maximum value of the 'left' parition should be less than or equal to the minimum value of the 'right'
    {
      scalar left_max = -std::numeric_limits<scalar>::infinity();
      for( const unsigned idx : left )
      {
        if( max_spread_axis.dot( vertices.col( idx ) ) > left_max )
        {
          left_max = max_spread_axis.dot( vertices.col( idx ) );
        }
      }
      scalar right_min = std::numeric_limits<scalar>::infinity();
      for( const unsigned idx : right )
      {
        if( max_spread_axis.dot( vertices.col( idx ) ) < right_min )
        {
          right_min = max_spread_axis.dot( vertices.col( idx ) );
        }
      }
      assert( left_max <= right_min );
    }
    {
      unsigned left_sum = 0;
      std::accumulate( std::begin(left), std::end(left), left_sum );
      unsigned right_sum = 0;
      std::accumulate( std::begin(right), std::end(right), right_sum );
      unsigned total_sum = 0;
      std::accumulate( std::begin(indices), std::end(indices), total_sum );
      assert( left_sum + right_sum == total_sum );
    }
    #endif
  }

  indices.clear();

  splitNode( vertices, left, 2 * node_index + 1 );
  splitNode( vertices, right, 2 * node_index + 2 );
}
