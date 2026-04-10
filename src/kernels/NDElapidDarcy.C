#include "NDElapidDarcy.h"

registerMooseObject("ElapidApp", NDElapidDarcy);

InputParameters
NDElapidDarcy::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Here the porosity is advected by Darcy velocity (Non-Dim)");

  params.addRequiredParam<RealVectorValue>("g_star", "gravity direction 0, 1, or -1");
  params.addRequiredParam<Real>(
      "rho_f", "Constant fluid density, ths density defines the non-dimensional pressure");

  params.addRequiredParam<Real>("phi_r",
                                "Reference porosity that determines the non-dimensional values");

  params.addRequiredCoupledVar("phi_f", "Porosity.");

  return params;
}

NDElapidDarcy::NDElapidDarcy(const InputParameters & parameters)
  : ADKernel(parameters),

    _g_star(getParam<RealVectorValue>("g_star")),
    _rho_f(getParam<Real>("rho_f")),
    _phi_r(getParam<Real>("phi_r")),

    _phi_f(adCoupledValue("phi_f")),

    _rho_s(getADMaterialProperty<Real>("rho_s"))

{
}

ADReal
NDElapidDarcy::computeQpResidual()
{
  // Not coupled anymore, this is used in the hydraulic equation, original negative at front
  const ADRealVectorValue darcy_velocity =
      -std::pow((_phi_f[_qp] / _phi_r), 3) *
      (_grad_u[_qp] + (_rho_f / (_rho_s[_qp] - _rho_f)) * _g_star);
  // should be negative?
  return -_grad_test[_i][_qp] * (darcy_velocity);
}
