#include "NDElapidSolidViscous.h"

registerMooseObject("ElapidApp", NDElapidSolidViscous);

InputParameters
NDElapidSolidViscous::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("This kernel defines the viscous compaction rheology.");

  params.addRequiredParam<Real>("zeta", "Constant linking eta_s0 to eta_c0");
  params.addRequiredParam<Real>("phi_r", "Reference porosity to which all scales are pinned");

  params.addRequiredCoupledVar("P_f", "Fluid pressure.");
  params.addRequiredCoupledVar("phi_f", "Fluid pressure.");

  return params;
}

NDElapidSolidViscous::NDElapidSolidViscous(const InputParameters & parameters)
  : ADKernel(parameters),

    _zeta(getParam<Real>("zeta")),
    _phi_r(getParam<Real>("phi_r")),

    _P_f(adCoupledValue("P_f")),
    _phi_f(adCoupledValue("phi_f"))

{
}

ADReal
NDElapidSolidViscous::computeQpResidual()
{
  return _test[_i][_qp] *
         ((_u[_qp] - _P_f[_qp]) / ((1 - _phi_f[_qp]) * _zeta * (_phi_r / _phi_f[_qp])));
}
