#pragma once

#include "ADKernel.h"

class NDElapidHydroViscous : public ADKernel
{
public:
  NDElapidHydroViscous(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  const Real & _zeta;  // zeta links the eta_s0 to eta_c0
  const Real & _phi_r; // the reference phi that pins all scaled vars

  // Coupled variables -- for this kernel we apply it to solve for P_tot so we do not need to add it
  const ADVariableValue & _P_tot;
  const ADVariableValue & _phi_f;
};
