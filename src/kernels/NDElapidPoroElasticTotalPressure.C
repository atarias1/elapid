#include "NDElapidPoroElasticTotalPressure.h"

registerMooseObject("ElapidApp", NDElapidPoroElasticTotalPressure);

InputParameters
NDElapidPoroElasticTotalPressure::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The time derivative operator with material property coefficient");

  params.addRequiredCoupledVar("P_tot", "Total pressure - this will be the time derivative");
  return params;
}

NDElapidPoroElasticTotalPressure::NDElapidPoroElasticTotalPressure(
    const InputParameters & parameters)
  : ADTimeKernel(parameters),

    _Pi4(getADMaterialProperty<Real>("p4")),
    _P_tot_dot(adCoupledDot("P_tot")),
    _K_phi_hat(getADMaterialProperty<Real>("K_phi_hat"))
{
}

ADReal
NDElapidPoroElasticTotalPressure::computeQpResidual()
{
  return _test[_i][_qp] * (((_Pi4[_qp] * (1.0 - _u[_qp])) / _K_phi_hat[_qp]) * _P_tot_dot[_qp]);
}
