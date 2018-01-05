#include "parse_params.h"
#include "implicitutilities/MathDefines.h"
#include "implicitutilities/ObjParser.h"
#include "implicitutilities/DistanceTools.h"
#include "implicitutilities/SignedDistanceField.h"
#include "implicitutilities/SurfaceSampling.h"
#include "implicitutilities/Moments.h"
#include "implicitutilities/SphereTree.h"
#include "marchingcubes/MarchingCubes.h"

#include "implicitutilities/MeshTools.h"
#include "implicitutilities/ProcessedMeshWriter.h"
#include "implicitutilities/Qhull.h"


#include <string>
#include <iostream>

// Mostly copied from implicittoolkitqt4/GLWidget.cpp
int main(int argc, const char * argv[])
{
  std::string USAGE = R"(USAGE:

./implicittoolkit_cli [options] input.obj output.h5

-w  followed by "cell width" {1}
-p  followed by "grid padding" {1}

For example:

./implicittoolkit_cli -w 0.2 -p 2 bunny.{obj,h5}
)";
  scalar cell_width = 1;
  scalar grid_padding = 1;
  std::string flags = "";
  int argi = 1;
  if(!igl::parse_params(argc,argv,argi,flags,
    'p',&grid_padding,
    'w',&cell_width))
  {
    std::cout<<"Error: failed to parse params"<<std::endl;
    std::cerr<<USAGE<<std::endl;
    return EXIT_FAILURE;
  }
  if(argi>argc-1)
  {
    std::cout<<"Error: failed to find input/output names"<<std::endl;
    std::cerr<<USAGE<<std::endl;
    return EXIT_FAILURE;
  }

  Matrix3Xsc new_vertices;
  Matrix2Xuc new_edges;
  Matrix3Xuc new_triangles;
  std::string input_file_name(argv[argi]);
  std::string output_file_name(argv[argi+1]);
  try
  {
    ObjParser::loadMesh(input_file_name, new_vertices, new_edges, new_triangles );
  }
  catch( const std::string& error )
  {
    std::cerr << "Failed to load " << input_file_name << ": " << error << std::endl;
    return EXIT_FAILURE;
  }

  // Compute moments
  Moments moments( new_vertices, new_triangles );
  // Translate and rotate the mesh to align the principal axes with the Cartesian axes at the origin
  new_vertices = moments.R().transpose() * ( new_vertices.colwise() - moments.x() );

  //// If everything checked out, save the mesh
  //m_input_vertices.swap( new_vertices );
  //m_input_edges.swap( new_edges );
  //m_input_triangles.swap( new_triangles );
  //// ... and save the moments
  //m_moments = std::move( moments );

  // Compute the vertices that compose the convex hull of the input vertices
  Matrix3Xsc convex_hull;
  // This printing some annoying things to the screen
  Qhull::computeConvexHull( new_vertices, convex_hull );
  SphereTree sphere_tree;
  SignedDistanceField sdf;

  const Eigen::Matrix<scalar,6,1> bbox = MeshTools::computeBoundingBox( convex_hull );
  scalar mesh_max_bbox_width  = ( bbox.segment<3>( 3 ) - bbox.segment<3>( 0 ) ).maxCoeff();

  Matrix3Xuc sdf_triangles;
  Matrix3Xsc sdf_vertices;
  sdf_triangles.resize( 3, 0 );
  sdf_vertices.resize( 3, 0 );
  sdf.clear();
  SurfaceSampling surface_sampling;
  Matrix3Xsc sample_normals;
  sample_normals.setZero();
  // Compute the grid for the implicit function
  assert( ( new_vertices.cols() == 0 && new_triangles.cols() == 0 ) || ( new_vertices.cols() != 0 && new_triangles.cols() != 0 ) );
  if( new_vertices.cols() == 0 ) return EXIT_FAILURE;
  const Eigen::Matrix<scalar,6,1> bounding_box = MeshTools::computeBoundingBox( new_vertices );
  sdf.initializeCenteredGrid( bounding_box, Vector3s::Constant( cell_width ), Vector3s::Constant( grid_padding ) );
  
  // Compute the initial surface sampling
  surface_sampling.computeSurfaceSampling( new_vertices, new_edges, new_triangles );
  sphere_tree.buildTree( surface_sampling.samples(), sdf.gridStart(), sdf.gridEnd() );
  sdf.evaluateGradients( surface_sampling.samples(), sample_normals );

  DistanceTools::computeSignedDistanceField( 
    sdf.gridStart(),
    sdf.gridDimensions(),
    sdf.cellDelta(),
    new_vertices,
    new_triangles,
    sdf.vals());

  ProcessedMeshWriter::write( 
    output_file_name,
    sdf,
    surface_sampling,
    moments,
    new_vertices,
    new_edges,
    new_triangles,
    convex_hull,
    cell_width,
    grid_padding,
    surface_sampling.numEdgeSubSamples(),
    surface_sampling.numFaceSubSamples(),
    sphere_tree );


  return EXIT_SUCCESS;
}
