#pragma once

#include "ADTimeKernel.h"

class NDElapidPoroElasticTotalPressure : public ADTimeKernel
{
public:
  NDElapidPoroElasticTotalPressure(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  const ADMaterialProperty<Real> & _Pi4; // Scaling factor Pstar /K_phir

  const ADVariableValue & _P_tot_dot;

  const ADMaterialProperty<Real> & _K_phi_hat;
};
