#include "NDElapidSolidElasticTotalPressure.h"

registerMooseObject("ElapidApp", NDElapidSolidElasticTotalPressure);

InputParameters
NDElapidSolidElasticTotalPressure::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The time derivative operator with material property coefficient");

  return params;
}

NDElapidSolidElasticTotalPressure::NDElapidSolidElasticTotalPressure(
    const InputParameters & parameters)
  : ADTimeKernel(parameters),

    _Pi1(getADMaterialProperty<Real>("p1")),
    _K_d_hat(getADMaterialProperty<Real>("K_d_hat"))
{
}

ADReal
NDElapidSolidElasticTotalPressure::computeQpResidual()
{
  return _test[_i][_qp] * ((_Pi1[_qp] / _K_d_hat[_qp]) * _u_dot[_qp]);
}
