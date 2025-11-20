#include "Phi_FTimeDerivative.h"

registerMooseObject("ElapidApp", Phi_FTimeDerivative);

InputParameters
Phi_FTimeDerivative::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription(
      "Here we find the time derivative of phi_f using the time derivatives of the solid species.");
  params.addRequiredCoupledVar("phi_ol", "Olivine vol fraction");
  params.addRequiredCoupledVar("phi_atg", "Antigorite vol fraction");

  return params;
}

Phi_FTimeDerivative::Phi_FTimeDerivative(const InputParameters & parameters)
  : ADKernel(parameters), _phi_ol_dot(adCoupledDot("phi_ol")), _phi_atg_dot(adCoupledDot("phi_atg"))
{
}

ADReal
Phi_FTimeDerivative::computeQpResidual()
{
  return _test[_i][_qp] * -(_phi_ol_dot[_qp] + _phi_atg_dot[_qp]);
}
