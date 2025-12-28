#pragma once

#include "ADKernel.h"

/**
 * Kernel = grad(test) * darcy_velocity * phi_f inspired by PorousFlowBasicAdvection kernel
 */

class ElapidDarcy : public ADKernel
{
public:
  ElapidDarcy(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // constant parameters
  const Real & _mu;
  const RealVectorValue & _gravity;
  const Real & _rho_f;

  // Material params
  const ADMaterialProperty<Real> & _k;
};
