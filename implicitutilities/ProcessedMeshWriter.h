#ifndef PROCESSED_MESH_WRITER_H
#define PROCESSED_MESH_WRITER_H

#include "MathDefines.h"

class SignedDistanceField;
class SurfaceSampling;
class Moments;
class SphereTree;

namespace ProcessedMeshWriter
{
  void write( const std::string& file_name, const SignedDistanceField& sdf, const SurfaceSampling& surface_sampling, const Moments& moments, const Matrix3Xsc& vertices, const Matrix2Xuc& edges, const Matrix3Xuc& faces, const Matrix3Xsc& convex_hull_vertices, const scalar& cell_width, const scalar& grid_padding, const unsigned edge_subsamples, const unsigned face_subsamples, const SphereTree& sphere_tree );
}

#endif
