#include "SinglePhaseLinearViscoElastic.h"
#include "libmesh/utility.h"

registerMooseObject("ElapidApp", SinglePhaseLinearViscoElastic);

InputParameters
SinglePhaseLinearViscoElastic::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<Real>("rho_f", "Fluid constant density.");

  params.addRequiredParam<Real>("rho_s", "Solid constant density.");

  params.addRequiredParam<Real>(
      "k_ref", "A reference permeability value that correlates with the reference porosity.");

  params.addRequiredParam<Real>(
      "phi_ref", "A reference porosity that correlates with the reference permeability.");

  params.addRequiredParam<Real>("nk", "The Kozeny-Carman exponent.");

  params.addRequiredParam<Real>("phi_0", "The background porosity.");

  params.addRequiredParam<Real>("mu", "Viscosity of fluid phase.");

  params.addRequiredParam<Real>("eta_s_0", "Constant nonporous solid viscosity.");

  params.addRequiredParam<Real>(
      "a_eta", "Coefficient to modify the effects of porosity on solid viscosity.");

  params.addRequiredParam<Real>("zeta", "Coefficient to control the compaction viscosity.");

  params.addRequiredParam<Real>("K", "Solid bulk modulus.");
  params.addRequiredParam<Real>("G", "Solid shear modulus.");

  params.addRequiredParam<Real>("aspect_ratio",
                                "Effective aspect ratio of pores in porous aggregates.");

  params.addRequiredParam<Real>("fluid_K", "Fluid bulk modulus.");

  params.addRequiredCoupledVar("phi_f", "Fluid pressure used to determine thermodynamic values.");

  return params;
}

SinglePhaseLinearViscoElastic::SinglePhaseLinearViscoElastic(const InputParameters & parameters)
  : Material(parameters),

    // Inputs
    _rho_f(getParam<Real>("rho_f")),
    _rho_s(getParam<Real>("rho_s")),
    _k_ref(getParam<Real>("k_ref")),
    _phi_ref(getParam<Real>("phi_ref")),
    _nk(getParam<Real>("nk")),
    _phi_0(getParam<Real>("phi_0")),
    _mu(getParam<Real>("mu")),
    _eta_s_0(getParam<Real>("eta_s_0")),
    _a_eta(getParam<Real>("a_eta")),
    _zeta(getParam<Real>("zeta")),
    _K(getParam<Real>("K")),
    _G(getParam<Real>("G")),
    _aspect_ratio(getParam<Real>("aspect_ratio")),
    _fluid_K(getParam<Real>("fluid_K")),

    // Get reference to coupled variables
    _phi_f(adCoupledValue("phi_f")),

    // To calculate

    // Declared properties
    // Densities & volume fractions
    _phi_s(declareADProperty<Real>("phi_s")),
    _rho_T(declareADProperty<Real>("rho_T")),

    // Explicitly calculated
    _k(declareADProperty<Real>("k")),

    // Viscous solid flow
    _eta_s(declareADProperty<Real>("eta_s")),
    _eta_compact(declareADProperty<Real>("eta_compact")),

    // Solid elasticity
    _K_d(declareADProperty<Real>("K_d")),
    _alpha(declareADProperty<Real>("alpha")),

    _Skempton(declareADProperty<Real>("Skempton")),
    _K_phi(declareADProperty<Real>("K_phi"))

{
}

void
SinglePhaseLinearViscoElastic::computeQpProperties()
{

  // Determine densities and volume fractions

  _phi_s[_qp] = 1.0 - _phi_f[_qp];

  _rho_T[_qp] = (_rho_f * _phi_f[_qp]) + (_rho_s * _phi_s[_qp]);

  // Fluid flow - find local permeability

  _k[_qp] = _k_ref * std::pow(_phi_f[_qp] / _phi_ref, _nk);

  // Porous solid viscosity

  _eta_s[_qp] = _eta_s_0 * std::exp(-_a_eta * (_phi_f[_qp] / _phi_0 - 1.0));

  _eta_compact[_qp] = _zeta * _eta_s[_qp] * (_phi_0 / _phi_f[_qp]);

  // Poro-elasticity

  // for Mori-Tanaka
  const ADReal beta_s = _G * (3.0 * _K + _G) / (3.0 * _K + 4.0 * _G);

  _K_d[_qp] =
      _K -
      _phi_f[_qp] * (std::pow(_K, 2.0) / (pi * _aspect_ratio * beta_s)) *
          std::pow((1 - _phi_f[_qp] + ((_phi_f[_qp] * _K) / (pi * _aspect_ratio * beta_s))), -1.0);

  _alpha[_qp] = 1.0 - (_K_d[_qp] / _K);

  _Skempton[_qp] = ((1.0 / _K_d[_qp]) - (1.0 / _K)) /
                   ((1.0 / _K_d[_qp]) - (1.0 / _K) + _phi_f[_qp] * ((1.0 / _fluid_K) - (1.0 / _K)));

  _K_phi[_qp] = 1.0 / (((1 - _phi_f[_qp]) / _K_d[_qp]) - (1 / _K));
}
