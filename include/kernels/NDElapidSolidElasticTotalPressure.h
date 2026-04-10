#pragma once

#include "ADTimeKernel.h"

class NDElapidSolidElasticTotalPressure : public ADTimeKernel
{
public:
  NDElapidSolidElasticTotalPressure(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // no coupled values for this class

  const ADMaterialProperty<Real> & _Pi1;     // Scaling factor Pstar/K_dr
  const ADMaterialProperty<Real> & _K_d_hat; // scaled K_d
};
