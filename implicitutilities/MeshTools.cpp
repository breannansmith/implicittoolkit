#include "MeshTools.h"

Eigen::Matrix<scalar,6,1> MeshTools::computeBoundingBox( const Matrix3Xsc& vertices )
{
  Eigen::Matrix<scalar,6,1> bounding_box;
  bounding_box.segment<3>( 0 ) = vertices.rowwise().minCoeff();
  bounding_box.segment<3>( 3 ) = vertices.rowwise().maxCoeff();
  assert( ( bounding_box.segment<3>( 0 ).array() <= bounding_box.segment<3>( 3 ).array() ).all() );
  return bounding_box;
}

int MeshTools::eulerCharacteristic( const Matrix3Xsc& vertices, const Matrix2Xuc& edges, const Matrix3Xuc& faces )
{
  return int(vertices.cols() - edges.cols() + faces.cols());
}
