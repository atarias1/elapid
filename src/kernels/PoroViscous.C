#include "PoroViscous.h"

registerMooseObject("ElapidApp", PoroViscous);

InputParameters
PoroViscous::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("This kernel defines the viscous compaction rheology.");

  params.addRequiredCoupledVar("P_f", "Fluid pressure.");
  params.addRequiredCoupledVar("P_tot", "Total pressure.");

  return params;
}

PoroViscous::PoroViscous(const InputParameters & parameters)
  : ADKernel(parameters),

  _P_f(adCoupledValue("P_f")),
  _P_tot(adCoupledValue("P_tot")),
  _comp_eta(getADMaterialProperty<Real>("comp_eta"))

{
}

ADReal
PoroViscous::computeQpResidual()
{
  return _test[_i][_qp] * -(((1.0 - _u[_qp])*(_P_f[_qp] - _P_tot[_qp])) / _comp_eta[_qp]);
}
