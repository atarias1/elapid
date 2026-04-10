#pragma once

#include "ADTimeKernel.h"

class NDElapidSolidElasticFluidPressure : public ADTimeKernel
{
public:
  NDElapidSolidElasticFluidPressure(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // Coupled variable time derivative
  // For our model our coupled variable is Pore pressure
  const ADVariableValue & _P_f_dot;

  const ADMaterialProperty<Real> & _Pi2; // Scaling factor Pstar alphar/K_dr

  const ADMaterialProperty<Real> & _alpha_hat; // Scaled Biot's
  const ADMaterialProperty<Real> & _K_d_hat;   // Scaled Drained modulus
};
