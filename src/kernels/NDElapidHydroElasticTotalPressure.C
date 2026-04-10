#include "NDElapidHydroElasticTotalPressure.h"

registerMooseObject("ElapidApp", NDElapidHydroElasticTotalPressure);

InputParameters
NDElapidHydroElasticTotalPressure::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The time derivative operator with material property coefficient");

  params.addRequiredCoupledVar("P_tot", "Total Pressure.");
  return params;
}

NDElapidHydroElasticTotalPressure::NDElapidHydroElasticTotalPressure(
    const InputParameters & parameters)
  : ADTimeKernel(parameters),

    _Pi2(getADMaterialProperty<Real>("p2")),
    _P_tot_dot(adCoupledDot("P_tot")),
    _K_d_hat(getADMaterialProperty<Real>("K_d_hat")),
    _alpha_hat(getADMaterialProperty<Real>("alpha_hat"))
{
}

ADReal
NDElapidHydroElasticTotalPressure::computeQpResidual()
{
  return _test[_i][_qp] * -((_Pi2[_qp] * _alpha_hat[_qp] / _K_d_hat[_qp]) * _P_tot_dot[_qp]);
}
