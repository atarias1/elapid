#pragma once

#include "ADTimeKernel.h"

class NDElapidPoroElasticFluidPressure : public ADTimeKernel
{
public:
  NDElapidPoroElasticFluidPressure(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  const ADMaterialProperty<Real> & _Pi4;

  const ADVariableValue & _P_f_dot;

  const ADMaterialProperty<Real> & _K_phi_hat;
};
