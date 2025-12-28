#pragma once

#include "ADTimeKernel.h"

class ElapidHydroElasticTotalPressure : public ADTimeKernel
{
public:
  ElapidHydroElasticTotalPressure(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // no coupled values for this class
  const ADVariableValue & _P_tot_dot;

  const ADMaterialProperty<Real> & _K_d;
  const ADMaterialProperty<Real> & _alpha;
};
