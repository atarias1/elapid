#include "ElapidSolidViscous.h"

registerMooseObject("ElapidApp", ElapidSolidViscous);

InputParameters
ElapidSolidViscous::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("This kernel defines the viscous compaction rheology.");

  params.addRequiredCoupledVar("P_f", "Fluid pressure.");
  params.addRequiredCoupledVar("phi_f", "Fluid pressure.");

  return params;
}

ElapidSolidViscous::ElapidSolidViscous(const InputParameters & parameters)
  : ADKernel(parameters),

    _P_f(adCoupledValue("P_f")),
    _phi_f(adCoupledValue("phi_f")),

    _eta_compact(getADMaterialProperty<Real>("eta_compact"))

{
}

ADReal
ElapidSolidViscous::computeQpResidual()
{
  return _test[_i][_qp] * ((_u[_qp] - _P_f[_qp]) / ((1 - _phi_f[_qp]) * _eta_compact[_qp]));
}
