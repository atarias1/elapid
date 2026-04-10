#pragma once

#include "ADTimeKernel.h"

class NDElapidHydroElasticTotalPressure : public ADTimeKernel
{
public:
  NDElapidHydroElasticTotalPressure(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  const ADMaterialProperty<Real> & _Pi2; // Scaling factor Pstar alphar/K_dr

  const ADVariableValue & _P_tot_dot;

  const ADMaterialProperty<Real> & _K_d_hat;
  const ADMaterialProperty<Real> & _alpha_hat;
};
