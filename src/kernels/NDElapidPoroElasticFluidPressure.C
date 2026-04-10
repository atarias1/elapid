#include "NDElapidPoroElasticFluidPressure.h"

registerMooseObject("ElapidApp", NDElapidPoroElasticFluidPressure);

InputParameters
NDElapidPoroElasticFluidPressure::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The time derivative operator with material property coefficient");
  params.addRequiredCoupledVar("P_f", "Fluid pressure - this will be the time derivative");
  return params;
}

NDElapidPoroElasticFluidPressure::NDElapidPoroElasticFluidPressure(
    const InputParameters & parameters)
  : ADTimeKernel(parameters),

    _Pi4(getADMaterialProperty<Real>("p4")),
    _P_f_dot(adCoupledDot("P_f")),
    _K_phi_hat(getADMaterialProperty<Real>("K_phi_hat"))
{
}

ADReal
NDElapidPoroElasticFluidPressure::computeQpResidual()
{
  return _test[_i][_qp] * -(((_Pi4[_qp] * (1.0 - _u[_qp])) / _K_phi_hat[_qp]) * _P_f_dot[_qp]);
}
