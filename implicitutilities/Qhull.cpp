#include "Qhull.h"

#include <iostream>

#ifdef USE_QHULL
extern "C"
{
#include "libqhull/libqhull.h"
}
#endif

void Qhull::computeConvexHull( const Matrix3Xsc& vertices, Matrix3Xsc& convex_hull )
{
#ifdef USE_QHULL
  static_assert( std::is_same<coordT, scalar>::value, "Error, Qhull coordT must be a scalar type." );
  static_assert( std::is_same<coordT, double>::value, "Error, Qhull coordT must be a double type." );

  char flags[]{ "qhull s Tcv" };
  const int exit_code= qh_new_qhull( 3, int(vertices.cols()), const_cast<scalar*>( vertices.data() ), false, flags, nullptr, stderr );
  if( exit_code == 0 )
  {
    convex_hull.resize( 3, qh num_vertices );
    unsigned vrt_num = 0;
    vertexT* vertex;
    FORALLvertices
    {
      convex_hull.col( vrt_num++ ) = Eigen::Map<Vector3s>( vertex->point );
    }
  }
  else
  {
    std::cerr << "Error, failed to compute convex hull" << std::endl;
  }
  qh_freeqhull( !qh_ALL );
  int curlong;
  int totlong;
  qh_memfreeshort( &curlong, &totlong );
  if( curlong != 0 || totlong != 0 )
  {
    std::cout << "qhull internal warning (Qhull::computeConvexHull): did not free " << curlong << " bytes of long memory (" << totlong << " pieces)" << std::endl;
  }
#else
  std::cerr << "Error, please rebuild with Qhull (http://www.qhull.org) support before executing Qhull::computeConvexHull." << std::endl;
  std::exit( EXIT_FAILURE );
#endif
}
