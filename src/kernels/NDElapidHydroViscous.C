#include "NDElapidHydroViscous.h"

registerMooseObject("ElapidApp", NDElapidHydroViscous);

InputParameters
NDElapidHydroViscous::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("This kernel defines the viscous compaction rheology.");

  params.addRequiredParam<Real>("zeta", "Constant linking eta_s0 to eta_c0");
  params.addRequiredParam<Real>("phi_r", "Reference porosity to which all scales are pinned");

  params.addRequiredCoupledVar("P_tot", "Total Pressure.");
  params.addRequiredCoupledVar("phi_f", "Total Pressure.");

  return params;
}

NDElapidHydroViscous::NDElapidHydroViscous(const InputParameters & parameters)
  : ADKernel(parameters),

    _zeta(getParam<Real>("zeta")),
    _phi_r(getParam<Real>("phi_r")),

    _P_tot(adCoupledValue("P_tot")),
    _phi_f(adCoupledValue("phi_f")),

    _eta_s(getADMaterialProperty<Real>("eta_s"))

{
}

ADReal
NDElapidHydroViscous::computeQpResidual()
{
  return _test[_i][_qp] * -((_P_tot[_qp] - _u[_qp]) /
                            ((1 - _phi_f[_qp]) * _zeta * _eta_s[_qp] * (_phi_r / _phi_f[_qp])));
}
