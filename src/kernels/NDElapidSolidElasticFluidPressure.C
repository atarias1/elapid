#include "NDElapidSolidElasticFluidPressure.h"

registerMooseObject("ElapidApp", NDElapidSolidElasticFluidPressure);

InputParameters
NDElapidSolidElasticFluidPressure::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription(
      "The time derivative operator with material property coefficient. NEGATIVELY SIGNED.");

  params.addRequiredCoupledVar("P_f", "Fluid pressure - this will be the time derivative");

  return params;
}

NDElapidSolidElasticFluidPressure::NDElapidSolidElasticFluidPressure(
    const InputParameters & parameters)
  : ADTimeKernel(parameters),

    _P_f_dot(adCoupledDot("P_f")),

    _Pi2(getADMaterialProperty<Real>("p2")),

    _K_d_hat(getADMaterialProperty<Real>("K_d_hat")),

    _alpha_hat(getADMaterialProperty<Real>("alpha_hat"))

{
}

ADReal
NDElapidSolidElasticFluidPressure::computeQpResidual()
{
  // here we add a negative sign in front of alpha
  return _test[_i][_qp] * -(((_Pi2[_qp] * _alpha_hat[_qp]) / _K_d_hat[_qp]) * _P_f_dot[_qp]);
}
