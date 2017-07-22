#ifndef SPHERE_TREE_H
#define SPHERE_TREE_H

#include "implicitutilities/MathDefines.h"
#include <unordered_map>

class SphereTree final
{

public:

  SphereTree();

  unsigned depth() const;
  const Matrix3Xsc& nodeCenters() const;
  const VectorXs& nodeRadii() const;

  void extractLeaves( VectorXu& node_ids, Matrix2Xuc& contents ) const;

  // Retrieves spheres at given depth
  void spheres( const unsigned depth, std::vector<Vector3s>& centers, std::vector<scalar>& radii ) const;

  void buildTree( const Matrix3Xsc& vertices, const Vector3s& grid_start, const Vector3s& grid_end );

private:

  void splitNode( const Matrix3Xsc& vertices, std::vector<unsigned>& indices, const unsigned node_index );

  unsigned m_tree_depth;
  Matrix3Xsc m_x;
  VectorXs m_r;
  std::unordered_map<unsigned,Vector2u> m_leaves;

};

#endif
