#pragma once

#include "ADTimeKernel.h"

class ElapidPoroElasticFluidPressure : public ADTimeKernel
{
public:
  ElapidPoroElasticFluidPressure(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // no coupled values for this class

  const ADVariableValue & _P_f_dot;

  const ADMaterialProperty<Real> & _K_phi;
};
