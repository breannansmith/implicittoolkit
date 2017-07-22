#ifndef MATH_DEFINES_H
#define MATH_DEFINES_H

#include <Eigen/Core>
#include <Eigen/StdVector>

namespace MathDefines
{
  template<typename T>
  constexpr T PI = T(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624);
}

using scalar = double;
using Vector2u = Eigen::Matrix<unsigned,2,1>;
using Vector3s = Eigen::Matrix<scalar,3,1>;
using Array3s = Eigen::Array<scalar,3,1>;
using Vector3u = Eigen::Matrix<unsigned,3,1>;
using Array3u = Eigen::Array<unsigned,3,1>;
using Vector3i = Eigen::Matrix<int,3,1>;
using VectorXs = Eigen::Matrix<scalar,Eigen::Dynamic,1>;
using ArrayXu = Eigen::Array<unsigned,Eigen::Dynamic,1>;
using VectorXu = Eigen::Matrix<unsigned,Eigen::Dynamic,1>;
using Matrix33s = Eigen::Matrix<scalar,3,3>;
using Matrix3Xsc = Eigen::Matrix<scalar,3,Eigen::Dynamic,Eigen::ColMajor>;
using Matrix2Xuc = Eigen::Matrix<unsigned,2,Eigen::Dynamic,Eigen::ColMajor>;
using Matrix3Xuc = Eigen::Matrix<unsigned,3,Eigen::Dynamic,Eigen::ColMajor>;

#endif
