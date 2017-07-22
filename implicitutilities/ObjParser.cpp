#include "ObjParser.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include <set>

#include "StringUtilities.h"

// Temporary storage used for detecting duplicates
class Edge final
{

public:

  Edge( const Vector2u& edge );

  bool operator<( const Edge& other ) const;

  const Vector2u& indices() const;

private:

  const Vector2u m_edge;

};

Edge::Edge( const Vector2u& edge )
: m_edge( edge )
{}

bool Edge::operator<( const Edge& other ) const
{
  return std::tie( m_edge.x(), m_edge.y() ) < std::tie( other.m_edge.x(), other.m_edge.y() );
}

const Vector2u& Edge::indices() const
{
  return m_edge;
}

static void parsevCommand( std::istringstream& command_stream, std::vector<Vector3s>& vertices )
{
  Vector3s vert;
  if( !( command_stream >> vert.x() >> vert.y() >> vert.z() ) )
  {
    throw std::string{ "Invalid vertex in obj file, vertices must have three real coordinates" };
  }
  vertices.push_back( vert );
}

static void extractVertex( const std::string& str, unsigned& idx )
{
  // Read the first vertex of this face
  std::vector<std::string> vertex_spec;
  StringUtilities::tokenize( str, '/', vertex_spec );
  if( vertex_spec.size() < 1 || vertex_spec.size() > 3 )
  {
    throw std::string{ "Invalid face in obj file" };
  }

  // Attempt to extract the vertex index
  std::istringstream xstream( vertex_spec[0] );
  if( !( xstream >> idx ) )
  {
    throw std::string{ "Invalid face in obj file" };
  }
}

static void parsefCommand( std::istringstream& command_stream, std::vector<Vector3u>& triangles )
{
  Vector3u triangle;

  std::vector<std::string> vert_strngs;
  std::string vert_cmmnd;
  while( command_stream >> vert_cmmnd )
  {
    vert_strngs.push_back( vert_cmmnd );
  }

  if( vert_strngs.size() != 3 )
  {
    throw std::string{ "Detected face with invalid number of vertices, only triangle meshes are supported" };
  }

  extractVertex( vert_strngs[0], triangle.x() );
  extractVertex( vert_strngs[1], triangle.y() );
  extractVertex( vert_strngs[2], triangle.z() );

  triangle.x() -= 1;
  triangle.y() -= 1;
  triangle.z() -= 1;

  triangles.push_back(triangle);
}

void ObjParser::loadMesh( const std::string& file_name, Matrix3Xsc& vertices, Matrix2Xuc& edges, Matrix3Xuc& faces )
{
  // Open and read the mesh file
  {
    // Attempt to open the user specified file
    std::ifstream obj_file;
    obj_file.open( file_name );
    if( !obj_file.is_open() )
    {
      throw std::string{ "Failed to open obj file" };
    }

    // Temporary storage for vertices and faces
    std::vector<Vector3s> local_vertices;
    std::vector<Vector3u> local_faces;

    // Parse any vertex or face lines
    std::string obj_command;
    while( !obj_file.eof() )
    {
      // Read a single line at a time
      getline( obj_file, obj_command );
  
      // Use a string stream for easy tokenizing
      std::istringstream commandstream( obj_command );

      // First element of a command is the command's name
      std::string command;
      commandstream >> command;

      // Vertex command
      if( command == "v" )
      {
        parsevCommand( commandstream, local_vertices );
      }
      // Face command
      else if( command == "f" )
      {
        parsefCommand( commandstream, local_faces );
      }
    }

    // Copy the the vertices to their final storage
    vertices.resize( 3, long(local_vertices.size()) );
    for( std::vector<Vector3s>::size_type vrt_idx = 0; vrt_idx < local_vertices.size(); ++vrt_idx )
    {
      vertices.col( long(vrt_idx) ) = local_vertices[ vrt_idx ];
    }

    // Copy the the faces to their final storage
    faces.resize( 3, long(local_faces.size()) );
    for( std::vector<Vector3u>::size_type fce_idx = 0; fce_idx < local_faces.size(); ++fce_idx )
    {
      faces.col( long(fce_idx) ) = local_faces[ fce_idx ];
    }
  }

  // Ensure all vertices of faces correspond to actual vertices
  assert( ( faces.array() < unsigned( vertices.size() ) ).all() );

  // Compute edges without duplicates
  {
    std::set<Edge> edge_set;
    for( int fce_idx = 0; fce_idx < faces.cols(); ++fce_idx )
    {
      const unsigned idx0 = faces( 0, fce_idx );
      const unsigned idx1 = faces( 1, fce_idx );
      const unsigned idx2 = faces( 2, fce_idx );
      // 0 - 1
      edge_set.insert( Edge( Vector2u( std::min( idx0, idx1 ), std::max( idx0, idx1 ) ) ) );
      // 1 - 2
      edge_set.insert( Edge( Vector2u( std::min( idx1, idx2 ), std::max( idx1, idx2 ) ) ) );
      // 2 - 0
      edge_set.insert( Edge( Vector2u( std::min( idx0, idx2 ), std::max( idx0, idx2 ) ) ) );
    }
    edges.resize( 2, long(edge_set.size()) );
    unsigned edge_num = 0;
    for( const Edge& edge : edge_set )
    {
      edges.col( edge_num++ ) = edge.indices();
    }
  }

  #ifndef NDEBUG
  for( unsigned edge_num0 = 0; edge_num0 < edges.cols(); ++edge_num0 )
  {
    for( unsigned edge_num1 = edge_num0 + 1; edge_num1 < edges.cols(); ++edge_num1 )
    {
      assert( edges.col( edge_num0 ) != edges.col( edge_num1 ) );
    }
  }
  #endif
}
