#pragma once

#include "ADKernel.h"

class ElapidViscousStress2D : public ADKernel
{
public:
  ElapidViscousStress2D(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  const unsigned int _component; // Provided to the kernel to decide which velocity to solve for

  const ADVariableValue & _V_x_s;
  const ADVariableValue & _V_y_s;
  const ADVariableValue & _P_tot;

  const ADVariableGradient & _grad_V_x_s;
  const ADVariableGradient & _grad_V_y_s;
  const ADVariableGradient & _grad_P_tot;

  const ADMaterialProperty<Real> & _eta_s; // consumed property needs const
};
