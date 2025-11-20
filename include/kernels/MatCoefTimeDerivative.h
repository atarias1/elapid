#pragma once

#include "ADTimeKernel.h"

class MatCoefTimeDerivative : public ADTimeKernel
{
public:
  MatCoefTimeDerivative(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // no coupled values for this class

  const ADMaterialProperty<Real> & _K_d;
};
