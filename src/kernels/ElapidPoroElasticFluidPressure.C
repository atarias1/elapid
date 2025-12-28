#include "ElapidPoroElasticFluidPressure.h"

registerMooseObject("ElapidApp", ElapidPoroElasticFluidPressure);

InputParameters
ElapidPoroElasticFluidPressure::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The time derivative operator with material property coefficient");
  params.addRequiredCoupledVar("P_f", "Fluid pressure - this will be the time derivative");
  return params;
}

ElapidPoroElasticFluidPressure::ElapidPoroElasticFluidPressure(const InputParameters & parameters)
  : ADTimeKernel(parameters),
    _P_f_dot(adCoupledDot("P_f")),
    _K_phi(getADMaterialProperty<Real>("K_phi"))
{
}

ADReal
ElapidPoroElasticFluidPressure::computeQpResidual()
{
  return _test[_i][_qp] * -(((1.0 - _u[_qp]) / _K_phi[_qp]) * _P_f_dot[_qp]);
}
