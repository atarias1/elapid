#pragma once

#include "ADTimeKernel.h"

class ElapidHydroElasticFluidPressure : public ADTimeKernel
{
public:
  ElapidHydroElasticFluidPressure(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // no coupled values for this class

  const ADMaterialProperty<Real> & _K_d;
  const ADMaterialProperty<Real> & _alpha;
  const ADMaterialProperty<Real> & _Skempton;
};
