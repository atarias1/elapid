#include "ElapidSolidElasticTotalPressure.h"

registerMooseObject("ElapidApp", ElapidSolidElasticTotalPressure);

InputParameters
ElapidSolidElasticTotalPressure::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The time derivative operator with material property coefficient");
  return params;
}

ElapidSolidElasticTotalPressure::ElapidSolidElasticTotalPressure(const InputParameters & parameters)
  : ADTimeKernel(parameters), _K_d(getADMaterialProperty<Real>("K_d"))
{
}

ADReal
ElapidSolidElasticTotalPressure::computeQpResidual()
{
  return _test[_i][_qp] * ((1.0 / _K_d[_qp]) * _u_dot[_qp]);
}
