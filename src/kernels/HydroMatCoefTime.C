#include "HydroMatCoefTime.h"

registerMooseObject("ElapidApp", HydroMatCoefTime);

InputParameters
HydroMatCoefTime::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The time derivative operator with material property coefficient");
  return params;
}

HydroMatCoefTime::HydroMatCoefTime(const InputParameters & parameters)
  : ADTimeKernel(parameters),
    _K_d(getADMaterialProperty<Real>("K_d")),
    _alpha(getADMaterialProperty<Real>("alpha")),
    _Skempton(getADMaterialProperty<Real>("alpha"))
{
}

ADReal
HydroMatCoefTime::computeQpResidual()
{
  return _test[_i][_qp] * ((_alpha[_qp] / (_K_d[_qp] * _Skempton[_qp])) * _u_dot[_qp]);
}
