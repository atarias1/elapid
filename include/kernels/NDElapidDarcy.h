#pragma once

#include "ADKernel.h"

/**
 * Kernel = grad(test) * darcy_velocity * phi_f inspired by PorousFlowBasicAdvection kernel
 */

class NDElapidDarcy : public ADKernel
{
public:
  NDElapidDarcy(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // constant parameters
  const RealVectorValue & _g_star;
  const Real & _rho_f;
  const Real & _phi_r;

  // couplings
  const ADVariableValue & _phi_f;

  const ADMaterialProperty<Real> & _rho_s;
};
