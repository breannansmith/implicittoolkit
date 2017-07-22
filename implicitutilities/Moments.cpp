#include "Moments.h"

#include <Eigen/Eigenvalues>
#include <iostream>

static void diagonalizeInertiaTensor( const Matrix33s& I, Matrix33s& R0, Vector3s& I0 )
{
  // Inertia tensor should by symmetric
  assert( ( I - I.transpose() ).lpNorm<Eigen::Infinity>() <= 1.0e-6 );
  // Inertia tensor should have positive determinant
  assert( I.determinant() > 0.0 );

  // Compute the eigenvectors and eigenvalues of the input matrix
  const Eigen::SelfAdjointEigenSolver<Matrix33s> es( I );

  // Check for errors
  if( es.info() == Eigen::NumericalIssue )
  {
    std::cerr << "Warning, failed to compute eigenvalues of inertia tensor due to Eigen::NumericalIssue" << std::endl;
  }
  else if( es.info() == Eigen::NoConvergence )
  {
    std::cerr << "Warning, failed to compute eigenvalues of inertia tensor due to Eigen::NoConvergence" << std::endl;
  }
  else if( es.info() == Eigen::InvalidInput )
  {
    std::cerr << "Warning, failed to compute eigenvalues of inertia tensor due to Eigen::InvalidInput" << std::endl;
  }
  assert( es.info() == Eigen::Success );

  // Save the eigenvectors and eigenvalues
  I0 = es.eigenvalues();
  assert( ( I0.array() > 0.0 ).all() ); assert( I0.x() <= I0.y() ); assert( I0.y() <= I0.z() );
  R0 = es.eigenvectors();
  assert( fabs( fabs( R0.determinant() ) - 1.0 ) <= 1.0e-6 );

  // TODO: Double check that this is always legit
  // Ensure that we have an orientation preserving transform
  if( R0.determinant() < 0.0 )
  {
    R0.col( 0 ) *= -1.0;
  }
}

// Polyhedral Mass Properties (Revisited) by David Eberly
// *** Caller must rescale mass and I by density ***
static void computeMoments( const Matrix3Xsc& vertices, const Matrix3Xuc& indices, scalar& volume, Vector3s& I_on_rho, Vector3s& center_of_mass, Matrix33s& R )
{
  assert( ( indices.array() < unsigned(vertices.cols()) ).all() );

  constexpr scalar oneDiv6 = 1.0 / 6.0;
  constexpr scalar oneDiv24 = 1.0 / 24.0;
  constexpr scalar oneDiv60 = 1.0 / 60.0;
  constexpr scalar oneDiv120 = 1.0 / 120.0;

  // order:  1, x, y, z, x^2, y^2, z^2, xy, yz, zx
  VectorXs integral = VectorXs::Zero( 10 );

  for( int i = 0; i < indices.cols(); ++i )
  {
    // Copy the vertices of triangle i
    const Vector3s v0 = vertices.col( indices(0,i) );
    const Vector3s v1 = vertices.col( indices(1,i) );
    const Vector3s v2 = vertices.col( indices(2,i) );

    // Compute a normal for the current triangle
    const Vector3s N = ( v1 - v0 ).cross( v2 - v0 );

    // Compute the integral terms
    scalar tmp0, tmp1, tmp2;
    scalar f1x, f2x, f3x, g0x, g1x, g2x;
    tmp0 = v0.x() + v1.x();
    f1x = tmp0 + v2.x();
    tmp1 = v0.x()*v0.x();
    tmp2 = tmp1 + v1.x()*tmp0;
    f2x = tmp2 + v2.x()*f1x;
    f3x = v0.x()*tmp1 + v1.x()*tmp2 + v2.x()*f2x;
    g0x = f2x + v0.x()*(f1x + v0.x());
    g1x = f2x + v1.x()*(f1x + v1.x());
    g2x = f2x + v2.x()*(f1x + v2.x());

    scalar f1y, f2y, f3y, g0y, g1y, g2y;
    tmp0 = v0.y() + v1.y();
    f1y = tmp0 + v2.y();
    tmp1 = v0.y()*v0.y();
    tmp2 = tmp1 + v1.y()*tmp0;
    f2y = tmp2 + v2.y()*f1y;
    f3y = v0.y()*tmp1 + v1.y()*tmp2 + v2.y()*f2y;
    g0y = f2y + v0.y()*(f1y + v0.y());
    g1y = f2y + v1.y()*(f1y + v1.y());
    g2y = f2y + v2.y()*(f1y + v2.y());

    scalar f1z, f2z, f3z, g0z, g1z, g2z;
    tmp0 = v0.z() + v1.z();
    f1z = tmp0 + v2.z();
    tmp1 = v0.z()*v0.z();
    tmp2 = tmp1 + v1.z()*tmp0;
    f2z = tmp2 + v2.z()*f1z;
    f3z = v0.z()*tmp1 + v1.z()*tmp2 + v2.z()*f2z;
    g0z = f2z + v0.z()*(f1z + v0.z());
    g1z = f2z + v1.z()*(f1z + v1.z());
    g2z = f2z + v2.z()*(f1z + v2.z());

    // Update integrals.
    integral[0] += N.x()*f1x;
    integral[1] += N.x()*f2x;
    integral[2] += N.y()*f2y;
    integral[3] += N.z()*f2z;
    integral[4] += N.x()*f3x;
    integral[5] += N.y()*f3y;
    integral[6] += N.z()*f3z;
    integral[7] += N.x()*(v0.y()*g0x + v1.y()*g1x + v2.y()*g2x);
    integral[8] += N.y()*(v0.z()*g0y + v1.z()*g1y + v2.z()*g2y);
    integral[9] += N.z()*(v0.x()*g0z + v1.x()*g1z + v2.x()*g2z);
  }

  integral[0] *= oneDiv6;
  integral[1] *= oneDiv24;
  integral[2] *= oneDiv24;
  integral[3] *= oneDiv24;
  integral[4] *= oneDiv60;
  integral[5] *= oneDiv60;
  integral[6] *= oneDiv60;
  integral[7] *= oneDiv120;
  integral[8] *= oneDiv120;
  integral[9] *= oneDiv120;

  // mass
  volume = integral[0];

  // center of mass
  center_of_mass = Vector3s( integral[1], integral[2], integral[3] ) / volume;

  // inertia relative to world origin
  R(0,0) = integral[5] + integral[6];
  R(0,1) = -integral[7];
  R(0,2) = -integral[9];
  R(1,0) = R(0,1);
  R(1,1) = integral[4] + integral[6];
  R(1,2) = -integral[8];
  R(2,0) = R(0,2);
  R(2,1) = R(1,2);
  R(2,2) = integral[4] + integral[5];

  // Comptue the inertia relative to the center of mass
  R(0,0) -= volume*(center_of_mass.y()*center_of_mass.y() + center_of_mass.z()*center_of_mass.z());
  R(0,1) += volume*center_of_mass.x()*center_of_mass.y();
  R(0,2) += volume*center_of_mass.z()*center_of_mass.x();
  R(1,0) = R(0,1);
  R(1,1) -= volume*(center_of_mass.z()*center_of_mass.z() + center_of_mass.x()*center_of_mass.x());
  R(1,2) += volume*center_of_mass.y()*center_of_mass.z();
  R(2,0) = R(0,2);
  R(2,1) = R(1,2);
  R(2,2) -= volume*(center_of_mass.x()*center_of_mass.x() + center_of_mass.y()*center_of_mass.y());

  // Diagonalize the inertia tensor
  Matrix33s R0;
  diagonalizeInertiaTensor( R, R0, I_on_rho );

  // Check that we actually diagonalized the inertia tensor
  assert( ( R0 * Matrix33s( I_on_rho.asDiagonal() ) * R0.transpose() - R ).lpNorm<Eigen::Infinity>() <= 1.0e-6 );
  assert( ( Matrix33s( I_on_rho.asDiagonal() ) - R0.transpose() * R * R0 ).lpNorm<Eigen::Infinity>() <= 1.0e-6 );
  R = R0;

  // All inertias should be positive
  assert( ( I_on_rho.array() > 0.0 ).all() );
  // Check that we have an orthonormal transformation
  assert( ( R * R.transpose() - Matrix33s::Identity() ).lpNorm<Eigen::Infinity>() <= 1.0e-6 );
  assert( fabs( R.determinant() - 1.0 ) <= 1.0e-6 );
}

Moments::Moments()
: m_volume( 0.0 )
, m_I_on_rho( Vector3s::Zero() )
, m_center_of_mass( Vector3s::Zero() )
, m_R( Matrix33s::Identity() )
{}

Moments::Moments( const Matrix3Xsc& vertices, const Matrix3Xuc& faces )
{
  computeMoments( vertices, faces, m_volume, m_I_on_rho, m_center_of_mass, m_R );
}

const scalar& Moments::volume() const
{
  return m_volume;
}

const Vector3s& Moments::I_on_rho() const
{
  return m_I_on_rho;
}

const Vector3s& Moments::x() const
{
  return m_center_of_mass;
}

const Matrix33s& Moments::R() const
{
  return m_R;
}
