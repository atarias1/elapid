#include "ElapidHydroElasticTotalPressure.h"

registerMooseObject("ElapidApp", ElapidHydroElasticTotalPressure);

InputParameters
ElapidHydroElasticTotalPressure::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The time derivative operator with material property coefficient");
  params.addRequiredCoupledVar("P_tot", "Total Pressure.");
  return params;
}

ElapidHydroElasticTotalPressure::ElapidHydroElasticTotalPressure(const InputParameters & parameters)
  : ADTimeKernel(parameters),
    _P_tot_dot(adCoupledDot("P_tot")),
    _K_d(getADMaterialProperty<Real>("K_d")),
    _alpha(getADMaterialProperty<Real>("alpha"))
{
}

ADReal
ElapidHydroElasticTotalPressure::computeQpResidual()
{
  return _test[_i][_qp] * -((_alpha[_qp] / _K_d[_qp]) * _P_tot_dot[_qp]);
}
