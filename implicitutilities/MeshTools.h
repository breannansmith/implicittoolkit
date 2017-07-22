#ifndef MESH_TOOLS_H
#define MESH_TOOLS_H

#include "MathDefines.h"

namespace MeshTools
{

  // Computes the bounding box of a set of vertices
  Eigen::Matrix<scalar,6,1> computeBoundingBox( const Matrix3Xsc& vertices );

  // Compute the Euler characteristic of a mesh
  int eulerCharacteristic( const Matrix3Xsc& vertices, const Matrix2Xuc& edges, const Matrix3Xuc& faces );

}

#endif
