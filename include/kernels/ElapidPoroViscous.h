#pragma once

#include "ADKernel.h"

class ElapidPoroViscous : public ADKernel
{
public:
  ElapidPoroViscous(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // no coupled values for this class

  const ADVariableValue & _P_f;
  const ADVariableValue & _P_tot;

  const ADMaterialProperty<Real> & _eta_compact;
};
