#include "MatCoefTimeDerivativeCoupled.h"

registerMooseObject("ElapidApp", MatCoefTimeDerivativeCoupled);

InputParameters
MatCoefTimeDerivativeCoupled::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription(
      "The time derivative operator with material property coefficient. NEGATIVELY SIGNED.");
  params.addRequiredCoupledVar("P_f", "Fluid pressure - this will be the time derivative");
  return params;
}

MatCoefTimeDerivativeCoupled::MatCoefTimeDerivativeCoupled(const InputParameters & parameters)
  : ADTimeKernel(parameters),
    _P_f_dot(adCoupledDot("P_f")),
    _K_d(getADMaterialProperty<Real>("K_d")),
    _alpha(getADMaterialProperty<Real>("alpha"))

{
}

ADReal
MatCoefTimeDerivativeCoupled::computeQpResidual()
{
  // here we add a negative sign in front of alpha
  return _test[_i][_qp] * -((_alpha[_qp] / _K_d[_qp]) * _P_f_dot[_qp]);
}
