#include "ElapidVelocityDiv2D.h"

registerMooseObject("ElapidApp", ElapidVelocityDiv2D);

InputParameters
ElapidVelocityDiv2D::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addRequiredCoupledVar("V_x_s", "x velocity");
  params.addRequiredCoupledVar("V_y_s", "y velocity");
  return params;
}

ElapidVelocityDiv2D::ElapidVelocityDiv2D(const InputParameters & parameters)
  : ADKernel(parameters),
    _V_x_s(adCoupledValue("V_x_s")),
    _V_y_s(adCoupledValue("V_y_s")),
    _grad_V_x_s(adCoupledGradient("V_x_s")),
    _grad_V_y_s(adCoupledGradient("V_y_s"))
{
}

ADReal
ElapidVelocityDiv2D::computeQpResidual()
{
  return -(_grad_test[_i][_qp](0) * _V_x_s[_qp] + _grad_test[_i][_qp](1) * _V_y_s[_qp]);
}
