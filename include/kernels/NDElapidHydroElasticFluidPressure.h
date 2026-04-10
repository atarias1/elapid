#pragma once

#include "ADTimeKernel.h"

class NDElapidHydroElasticFluidPressure : public ADTimeKernel
{
public:
  NDElapidHydroElasticFluidPressure(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // no coupled values for this class

  const ADMaterialProperty<Real> & _Pi3; // Scaling factor (Pstar alphar)/(K_dr Br)

  const ADMaterialProperty<Real> & _K_d_hat;
  const ADMaterialProperty<Real> & _alpha_hat;
  const ADMaterialProperty<Real> & _Skempton_hat;
};
