#ifndef MOMENTS_H
#define MOMENTS_H

#include "MathDefines.h"

class Moments final
{

public:

  Moments();
  Moments( const Matrix3Xsc& vertices, const Matrix3Xuc& faces );

  const scalar& volume() const;
  const Vector3s& I_on_rho() const;
  const Vector3s& x() const;
  const Matrix33s& R() const;

private:

  scalar m_volume;
  Vector3s m_I_on_rho;
  Vector3s m_center_of_mass;
  Matrix33s m_R;

};

#endif
