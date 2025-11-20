#pragma once

#include "ADTimeKernel.h"

class Phi_FTimeDerivative : public ADKernel
{
public:
  Phi_FTimeDerivative(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // We have to be sneaky here and use the time derivatives of the tracked phases
  const ADVariableValue & _phi_ol_dot;
  const ADVariableValue & _phi_atg_dot;
};
