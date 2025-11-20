#pragma once

#include "ADKernel.h"

class PoroViscous : public ADKernel
{
public:
  PoroViscous(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  // no coupled values for this class

  const ADVariableValue & _P_f;
  const ADVariableValue & _P_tot;

  const ADMaterialProperty<Real> & _comp_eta;
};
