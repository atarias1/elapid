#include "MatCoefTimeDerivative.h"

registerMooseObject("ElapidApp", MatCoefTimeDerivative);

InputParameters
MatCoefTimeDerivative::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The time derivative operator with material property coefficient");
  return params;
}

MatCoefTimeDerivative::MatCoefTimeDerivative(const InputParameters & parameters)
  : ADTimeKernel(parameters), _K_d(getADMaterialProperty<Real>("K_d"))
{
}

ADReal
MatCoefTimeDerivative::computeQpResidual()
{
  return _test[_i][_qp] * ((1.0 / _K_d[_qp]) * _u_dot[_qp]);
}
