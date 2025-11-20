#pragma once

#include "ADTimeKernel.h"

class MatCoefTimeDerivativeCoupled : public ADTimeKernel
{
public:
  MatCoefTimeDerivativeCoupled(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // Coupled variable time derivative
  // For our model our coupled variable is Pore pressure
  const ADVariableValue & _P_f_dot;

  const ADMaterialProperty<Real> & _alpha; // Biot's
  const ADMaterialProperty<Real> & _K_d;   // Drained modulus
};
