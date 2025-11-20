#include "VelocityDiv2D.h"

registerMooseObject("ElapidApp", VelocityDiv2D);

InputParameters
VelocityDiv2D::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addRequiredCoupledVar("V_x_s", "x velocity");
  params.addRequiredCoupledVar("V_y_s", "y velocity");
  return params;
}

VelocityDiv2D::VelocityDiv2D(const InputParameters & parameters)
  : ADKernel(parameters),
    _V_x_s(adCoupledValue("V_x_s")),
    _V_y_s(adCoupledValue("V_y_s")),
    _grad_V_x_s(adCoupledGradient("V_x_s")),
    _grad_V_y_s(adCoupledGradient("V_y_s"))
{
}

ADReal
VelocityDiv2D::computeQpResidual()
{
  // should be negative if using grad_test
  return -(_grad_test[_i][_qp](0) * _V_x_s[_qp] + _grad_test[_i][_qp](1) * _V_y_s[_qp]);
  // return (_grad_V_x_s[_qp](0) + _grad_V_y_s[_qp](1)) * _test[_i][_qp];
}
