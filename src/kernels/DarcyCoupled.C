#include "DarcyCoupled.h"

registerMooseObject("ElapidApp", DarcyCoupled);

InputParameters
DarcyCoupled::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Here the porosity is advected by Darcy velocity");

  params.addRequiredParam<Real>("mu", "Fluid viscosity");
  params.addRequiredParam<RealVectorValue>("gravity", "Gravity");
  params.addRequiredParam<Real>("rho_f", "Fluid density");

  return params;
}

DarcyCoupled::DarcyCoupled(const InputParameters & parameters)
  : ADKernel(parameters),

    _mu(getParam<Real>("mu")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _rho_f(getParam<Real>("rho_f")),

    _k(getADMaterialProperty<Real>("k"))

{
}

ADReal
DarcyCoupled::computeQpResidual()
{
  // Not coupled anymore, this is used in the hydraulic equation, original negative at front
  const ADRealVectorValue darcy_velocity = -(_k[_qp] / _mu) * (_grad_u[_qp] + _rho_f * _gravity);
  // should be negative?
  return -_grad_test[_i][_qp] * (darcy_velocity);
}
