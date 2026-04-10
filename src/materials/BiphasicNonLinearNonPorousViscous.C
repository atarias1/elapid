#include "BiphasicNonLinearNonPorousViscous.h"
#include "libmesh/utility.h"

registerMooseObject("ElapidApp", BiphasicNonLinearNonPorousViscous);

InputParameters
BiphasicNonLinearNonPorousViscous::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<Real>("rho_x1", "Phase 1 constant density.");

  params.addRequiredParam<Real>("rho_x2", "Phase 2 constant density.");

  params.addRequiredParam<Real>("max_eta_s", "Maximum viscosity in case of low strain rate");

  params.addRequiredParam<Real>("x1_eta_0", "Phase 1 viscosity at strain rate of 1/s.");

  params.addRequiredParam<Real>("x2_eta_0", "Phase 2 viscosity at strain rate 1/s.");

  params.addRequiredParam<Real>("nsigma_x1", "Phase 1 stress exponent.");

  params.addRequiredParam<Real>("nsigma_x2", "Phase 2 stress exponent.");

  params.addRequiredCoupledVar("phi_x2", "Phase x2 volume fraction.");

  params.addRequiredCoupledVar("v_x", "X velocity to find strain rate.");

  params.addRequiredCoupledVar("v_y", "Y velocity to find strain rate.");

  return params;
}

BiphasicNonLinearNonPorousViscous::BiphasicNonLinearNonPorousViscous(
    const InputParameters & parameters)
  : Material(parameters),

    // Inputs
    _rho_x1(getParam<Real>("rho_x1")),
    _rho_x2(getParam<Real>("rho_x2")),
    _max_eta_s(getParam<Real>("max_eta_s")),
    _x1_eta_0(getParam<Real>("x1_eta_0")),
    _x2_eta_0(getParam<Real>("x2_eta_0")),
    _nsigma_x1(getParam<Real>("nsigma_x1")),
    _nsigma_x2(getParam<Real>("nsigma_x2")),

    // Get reference to coupled variables
    _phi_x2(adCoupledValue("phi_x2")),
    _v_x(adCoupledValue("v_x")),
    _v_y(adCoupledValue("v_y")),
    _grad_v_x(adCoupledGradient("v_x")),
    _grad_v_y(adCoupledGradient("v_y")),

    // To calculate

    // Declared properties
    // Densities & volume fractions
    _phi_x1(declareADProperty<Real>("phi_x1")),
    _rho_T(declareADProperty<Real>("rho_T")),

    // Viscous solid flow
    _eta_x1_eff(declareADProperty<Real>("eta_x1_eff")),
    _eta_x2_eff(declareADProperty<Real>("eta_x2_eff")),

    _eta_s(declareADProperty<Real>("eta_s"))

{
}

void
BiphasicNonLinearNonPorousViscous::computeQpProperties()
{

  // Determine densities and volume fractions first

  _phi_x1[_qp] = 1.0 - _phi_x2[_qp];

  _rho_T[_qp] = (_rho_x1 * _phi_x1[_qp]) + (_rho_x2 * _phi_x2[_qp]);

  // Solid viscosity

  // Deviatoric strain rate tensor
  // Velocity gradient components

  const ADReal div_v = _grad_v_x[_qp](0) + _grad_v_y[_qp](1);

  // Effective strain rate from Dev strain rate
  // make sure nums are floats even in the expressions
  const ADReal Dxx2 = std::pow(_grad_v_x[_qp](0) - (1.0 / 3.0) * div_v, 2.0);
  const ADReal Dyy2 = std::pow(_grad_v_y[_qp](1) - (1.0 / 3.0) * div_v, 2.0);
  const ADReal Dzz2 = std::pow(-(1.0 / 3.0) * div_v, 2.0);
  const ADReal Dxy2 = std::pow(0.5 * (_grad_v_x[_qp](1) + _grad_v_y[_qp](0)), 2.0);

  // need to add small value to prevent NoN
  const ADReal eff_strain_rate = std::sqrt(0.5 * (Dxx2 + Dyy2 + Dzz2 + 2 * Dxy2) + 1.0e-32);

  // effective viscosities for each phase - for biphasic materials with one linear phase you can
  // also use a modified MetaSerpentinite material - there, the eff viscosity of olivine is linear
  // (Diffusion creep @ constant grain size)

  _eta_x1_eff[_qp] = _x1_eta_0 * std::pow(eff_strain_rate, (1.0 / _nsigma_x1) - 1.0);
  _eta_x2_eff[_qp] = _x2_eta_0 * std::pow(eff_strain_rate, (1.0 / _nsigma_x2) - 1.0);

  // non-porous effective bulk viscosity -> maybe later we can try different mixing rules but here
  // we just do volume weighted average
  const ADReal eta_s_0 = ((_eta_x1_eff[_qp] * _phi_x1[_qp]) + (_eta_x2_eff[_qp] * _phi_x2[_qp]));

  // cap the solid viscosity and modify viscosity according to the local porosity

  if (eta_s_0 < _max_eta_s)
  {
    _eta_s[_qp] = eta_s_0;
  }

  else
  { // Here we use _max_eta_s if the predicted viscosity gets too high
    _eta_s[_qp] = _max_eta_s;
  }
}
