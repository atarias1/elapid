#include "NDElapidHydroElasticFluidPressure.h"

registerMooseObject("ElapidApp", NDElapidHydroElasticFluidPressure);

InputParameters
NDElapidHydroElasticFluidPressure::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The time derivative operator with material property coefficient");

  return params;
}

NDElapidHydroElasticFluidPressure::NDElapidHydroElasticFluidPressure(
    const InputParameters & parameters)
  : ADTimeKernel(parameters),

    _Pi3(getADMaterialProperty<Real>("p3")),
    _K_d_hat(getADMaterialProperty<Real>("K_d_hat")),
    _alpha_hat(getADMaterialProperty<Real>("alpha_hat")),
    _Skempton_hat(getADMaterialProperty<Real>("Skempton_hat"))
{
}

ADReal
NDElapidHydroElasticFluidPressure::computeQpResidual()
{
  return _test[_i][_qp] *
         (((_Pi3[_qp] * _alpha_hat[_qp]) / (_K_d_hat[_qp] * _Skempton_hat[_qp])) * _u_dot[_qp]);
}
