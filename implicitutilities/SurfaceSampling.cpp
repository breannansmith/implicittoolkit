#include "SurfaceSampling.h"

SurfaceSampling::SurfaceSampling()
: m_edge_subsamples( 0 )
, m_face_subsamples( 0 )
, m_samples()
{}

static int intPower( const int x, const int power )
{
  int total = 1;
  for( int count = 1; count <= power; ++count )  total *= x;
  return total;
}

static void subsampleFace( const unsigned recursion_level, const Vector3s& a, const Vector3s& b, const Vector3s& c, unsigned& current_sample, Matrix3Xsc& samples )
{
  if( recursion_level == 0 )
  {
    return;
  }
  const Vector3s d = ( a + b + c ) / 3.0;
  samples.col( current_sample++ ) = d;
  subsampleFace( recursion_level - 1, a, b, d, current_sample, samples );
  subsampleFace( recursion_level - 1, b, c, d, current_sample, samples );
  subsampleFace( recursion_level - 1, c, a, d, current_sample, samples );
}

void SurfaceSampling::computeSurfaceSampling( const Matrix3Xsc& vertices, const Matrix2Xuc& edges, const Matrix3Xuc& faces )
{
  const unsigned num_edge_samples = unsigned(m_edge_subsamples * edges.cols());
  const unsigned num_face_samples = unsigned(( unsigned( intPower( 3, int(m_face_subsamples) ) - 1 ) / 2 ) * faces.cols());
  const unsigned num_samples = unsigned(vertices.cols() + num_edge_samples + num_face_samples);

  m_samples.resize( 3, num_samples );
  m_samples.block( 0, 0, 3, vertices.cols() ) = vertices;

  unsigned current_sample = unsigned(vertices.cols());

  // Subsample edges
  {
    // Parametric spacing between edge subsamples
    const scalar dalpha = 1.0 / scalar( m_edge_subsamples + 1 );

    // For each edge
    for( int edge_num = 0; edge_num < edges.cols(); ++edge_num )
    {
      // Grab the endpoints of the edge
      const Vector3s x0 = vertices.col( edges( 0, edge_num ) );
      const Vector3s x1 = vertices.col( edges( 1, edge_num ) );
      // For each subsample
      for( unsigned i = 1; i <= m_edge_subsamples; ++i )
      {
        const scalar alpha = i * dalpha;
        assert( alpha > 0.0 ); assert( alpha < 1.0 );
        // Generate the sample point
        const Vector3s new_sample = ( 1.0 - alpha ) * x0 + alpha * x1;
        m_samples.col( current_sample++ ) = new_sample;
      }
    }
  }

  // For each triangle
  for( int tri_num = 0; tri_num < faces.cols(); ++tri_num )
  {
    const Vector3s x0 = vertices.col( faces( 0, tri_num ) );
    const Vector3s x1 = vertices.col( faces( 1, tri_num ) );
    const Vector3s x2 = vertices.col( faces( 2, tri_num ) );
    subsampleFace( m_face_subsamples, x0, x1, x2, current_sample, m_samples );
  }
}

const Matrix3Xsc& SurfaceSampling::samples() const
{
  return m_samples;
}

bool SurfaceSampling::empty() const
{
  return m_samples.cols() == 0;
}

unsigned SurfaceSampling::numSamples() const
{
  return unsigned( m_samples.cols() );
}

unsigned SurfaceSampling::numEdgeSubSamples() const
{
  return m_edge_subsamples;
}

unsigned SurfaceSampling::numFaceSubSamples() const
{
  return m_face_subsamples;
}

void SurfaceSampling::setEdgeSubsamples( const unsigned num_samples )
{
  m_edge_subsamples = num_samples;
}

void SurfaceSampling::setFaceSubsamples( const unsigned num_samples )
{
  m_face_subsamples = num_samples;
}
