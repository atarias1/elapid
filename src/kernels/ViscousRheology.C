#include "ViscousRheology.h"

registerMooseObject("ElapidApp", ViscousRheology);

InputParameters
ViscousRheology::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("This kernel defines the viscous compaction rheology.");

  params.addRequiredCoupledVar("P_f", "Fluid pressure.");
  params.addRequiredCoupledVar("phi_f", "Fluid pressure.");

  return params;
}

ViscousRheology::ViscousRheology(const InputParameters & parameters)
  : ADKernel(parameters),

    _P_f(adCoupledValue("P_f")),
    _phi_f(adCoupledValue("phi_f")),

    _comp_eta(getADMaterialProperty<Real>("comp_eta"))

{
}

ADReal
ViscousRheology::computeQpResidual()
{
  return _test[_i][_qp] * ((_u[_qp] - _P_f[_qp]) / ((1 - _phi_f[_qp]) * _comp_eta[_qp]));
}
