#include "PoroMechanoMatCoefTime.h"

registerMooseObject("ElapidApp", PoroMechanoMatCoefTime);

InputParameters
PoroMechanoMatCoefTime::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The time derivative operator with material property coefficient");
  params.addRequiredCoupledVar("P_tot", "Total pressure - this will be the time derivative");
  return params;
}

PoroMechanoMatCoefTime::PoroMechanoMatCoefTime(const InputParameters & parameters)
  : ADTimeKernel(parameters),
    _P_tot_dot(adCoupledDot("P_tot")),
    _K_phi(getADMaterialProperty<Real>("K_phi"))
{
}

ADReal
PoroMechanoMatCoefTime::computeQpResidual()
{
  return _test[_i][_qp] * (((1.0-_u[_qp])/_K_phi[_qp]) * _P_tot_dot[_qp]);
}
