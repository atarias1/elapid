#pragma once

#include "ADKernel.h"

class ElapidSolidViscous : public ADKernel
{
public:
  ElapidSolidViscous(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // Coupled variables -- for this kernel we apply it to solve for P_tot so we do not need to add it
  const ADVariableValue & _P_f;
  const ADVariableValue & _phi_f;

  // Material Properties
  const ADMaterialProperty<Real> & _eta_compact;
};
