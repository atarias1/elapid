#pragma once

#include "ADTimeKernel.h"

class HydroMatCoefTimeCoupled : public ADTimeKernel
{
public:
  HydroMatCoefTimeCoupled(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // no coupled values for this class
  const ADVariableValue & _P_tot_dot;

  const ADMaterialProperty<Real> & _K_d;
  const ADMaterialProperty<Real> & _alpha;
};
