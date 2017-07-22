#include "ProcessedMeshWriter.h"

#include "HDF5File.h"
#include "SignedDistanceField.h"
#include "SurfaceSampling.h"
#include "Moments.h"
#include "SphereTree.h"

void ProcessedMeshWriter::write( const std::string& file_name, const SignedDistanceField& sdf, const SurfaceSampling& surface_sampling, const Moments& moments, const Matrix3Xsc& vertices, const Matrix2Xuc& edges, const Matrix3Xuc& faces, const Matrix3Xsc& convex_hull_vertices, const scalar& cell_width, const scalar& grid_padding, const unsigned edge_subsamples, const unsigned face_subsamples, const SphereTree& sphere_tree )
{
  HDF5File output( file_name, HDF5AccessType::READ_WRITE );
  // Create a group for the signed distance field
  output.findOrCreateGroup( "sdf" );
  output.write( "sdf/grid_dimensions", sdf.gridDimensions() );
  output.write( "sdf/cell_delta", sdf.cellDelta() );
  output.write( "sdf/grid_origin", sdf.gridStart() );
  output.write( "sdf/signed_distance", sdf.vals() );
  // Create a group for the surface samples
  output.findOrCreateGroup( "surface_samples" );
  output.write( "surface_samples/samples", surface_sampling.samples() );
  // Create a group for the convex hull vertices
  output.findOrCreateGroup( "convex_hull" );
  output.write( "convex_hull/vertices", convex_hull_vertices );
  // Create a group for the moments
  output.findOrCreateGroup( "moments" );
  output.write( "moments/volume", moments.volume() );
  output.write( "moments/I_on_rho", moments.I_on_rho() );
  output.write( "moments/x", moments.x() );
  output.write( "moments/R", moments.R() );
  // Create a group for the input mesh
  output.findOrCreateGroup( "mesh" );
  output.write( "mesh/vertices", vertices );
  output.write( "mesh/edges", edges );
  output.write( "mesh/faces", faces );
  // Output settings used to generate this sdf and sampling
  output.findOrCreateGroup( "settings" );
  output.write( "settings/cell_width", cell_width );
  output.write( "settings/grid_padding", grid_padding );
  output.write( "settings/edge_subsamples", edge_subsamples );
  output.write( "settings/face_subsamples", face_subsamples );
  // Output the sphere tree
  output.findOrCreateGroup( "sphere_tree" );
  output.write( "sphere_tree/depth", sphere_tree.depth() );
  output.write( "sphere_tree/node_centers", sphere_tree.nodeCenters() );
  output.write( "sphere_tree/node_radii", sphere_tree.nodeRadii() );
  {
    VectorXu node_ids;
    Matrix2Xuc contents;
    sphere_tree.extractLeaves( node_ids, contents );
    output.write( "sphere_tree/leaf_node_ids", node_ids );
    output.write( "sphere_tree/leaf_contents", contents );
  }
}
