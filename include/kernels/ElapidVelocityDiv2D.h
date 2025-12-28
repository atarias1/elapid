#pragma once

#include "ADKernel.h"

class ElapidVelocityDiv2D : public ADKernel
{
public:
  ElapidVelocityDiv2D(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  const ADVariableValue & _V_x_s;
  const ADVariableValue & _V_y_s;
  const ADVariableGradient & _grad_V_x_s;
  const ADVariableGradient & _grad_V_y_s;
};
