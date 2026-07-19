#pragma once

#include "ADKernel.h"

class NDElapidPoroViscous : public ADKernel
{
public:
  NDElapidPoroViscous(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  const Real & _zeta;  // zeta links the eta_s0 to eta_c0
  const Real & _phi_r; // the reference phi that pins all scaled vars

  const ADVariableValue & _P_f;
  const ADVariableValue & _P_tot;
  const ADMaterialProperty<Real> & _eta_s;
};
