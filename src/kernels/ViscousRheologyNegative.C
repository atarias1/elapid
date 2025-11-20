#include "ViscousRheologyNegative.h"

registerMooseObject("ElapidApp", ViscousRheologyNegative);

InputParameters
ViscousRheologyNegative::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("This kernel defines the viscous compaction rheology.");

  params.addRequiredCoupledVar("P_tot", "Total Pressure.");
  params.addRequiredCoupledVar("phi_f", "Total Pressure.");


  return params;
}

ViscousRheologyNegative::ViscousRheologyNegative(const InputParameters & parameters)
  : ADKernel(parameters),

    _P_tot(adCoupledValue("P_tot")),
    _phi_f(adCoupledValue("phi_f")),

    _comp_eta(getADMaterialProperty<Real>("comp_eta"))

{
}

ADReal
ViscousRheologyNegative::computeQpResidual()
{
  return _test[_i][_qp] * -((_P_tot[_qp] - _u[_qp]) / ((1 - _phi_f[_qp]) * _comp_eta[_qp]));
}
