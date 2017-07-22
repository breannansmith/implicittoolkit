#ifndef Q_HULL_H
#define Q_HULL_H

#include "MathDefines.h"

// TODO: Figure out how to silence qhull status output
// TODO: Requries static linking, fix up build system to allow for dynamic
namespace Qhull
{

  void computeConvexHull( const Matrix3Xsc& vertices, Matrix3Xsc& convex_hull );

}

#endif
