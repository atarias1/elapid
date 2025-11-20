#pragma once

#include "ADTimeKernel.h"

class PoroMechanoMatCoefTime : public ADTimeKernel
{
public:
  PoroMechanoMatCoefTime(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // no coupled values for this class

  const ADVariableValue & _P_tot_dot;

  const ADMaterialProperty<Real> & _K_phi;
};
