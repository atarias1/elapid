#include "NDElapidPoroViscous.h"

registerMooseObject("ElapidApp", NDElapidPoroViscous);

InputParameters
NDElapidPoroViscous::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("This kernel defines the viscous compaction rheology.");

  params.addRequiredParam<Real>("zeta", "Constant linking eta_s0 to eta_c0");
  params.addRequiredParam<Real>("phi_r", "Reference porosity to which all scales are pinned");

  params.addRequiredCoupledVar("P_f", "Fluid pressure.");
  params.addRequiredCoupledVar("P_tot", "Total pressure.");

  return params;
}

NDElapidPoroViscous::NDElapidPoroViscous(const InputParameters & parameters)
  : ADKernel(parameters),

    _zeta(getParam<Real>("zeta")),
    _phi_r(getParam<Real>("phi_r")),

    _P_f(adCoupledValue("P_f")),
    _P_tot(adCoupledValue("P_tot")),
    _eta_s(getADMaterialProperty<Real>("eta_s"))

{
}

ADReal
NDElapidPoroViscous::computeQpResidual()
{
  // return _test[_i][_qp] *
  //        -(((1.0 - _u[_qp]) * (_P_f[_qp] - _P_tot[_qp])) / (_zeta * (_phi_r / _u[_qp])));

  return _test[_i][_qp] * -((_P_f[_qp] - _P_tot[_qp]) / (_zeta * _eta_s[_qp] * (_phi_r / _u[_qp])));
}
// Try removing 1 - phi_f
