#include "BiphasicLinearViscoElastic.h"
#include "libmesh/utility.h"

registerMooseObject("ElapidApp", BiphasicLinearViscoElastic);

InputParameters
BiphasicLinearViscoElastic::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<Real>("rho_f", "Fluid constant density.");

  params.addRequiredParam<Real>("rho_x1", "Phase 1 constant density.");

  params.addRequiredParam<Real>("rho_x2", "Phase 2 constant density.");

  params.addRequiredParam<Real>(
      "k_ref", "A reference permeability value that correlates with the reference porosity.");

  params.addRequiredParam<Real>(
      "phi_ref", "A reference porosity that correlates with the reference permeability.");

  params.addRequiredParam<Real>("nk", "The Kozeny-Carman exponent.");

  params.addRequiredParam<Real>("phi_0", "The initial ambient porosity.");

  params.addRequiredParam<Real>("mu", "Viscosity of water.");

  params.addRequiredParam<Real>("x1_eta_0", "Phase 1 constant viscosity.");

  params.addRequiredParam<Real>("x2_eta_0", "Phase 2 constant viscosity.");

  params.addRequiredParam<Real>("a_eta",
                                "Coefficient to modify the porosity effects on solid viscosity.");

  params.addRequiredParam<Real>("zeta", "Coefficient to control the compaction viscosity.");

  params.addRequiredParam<Real>("x1_K", "Phase 1 bulk modulus.");
  params.addRequiredParam<Real>("x2_K", "Phase 2 bulk modulus.");
  params.addRequiredParam<Real>("x1_G", "Phase 1 shear modulus.");
  params.addRequiredParam<Real>("x2_G", "Phase 2 shear modulus.");

  params.addRequiredParam<Real>("x1_aspect",
                                "Effective aspect ratio of pores in porous phase 1 aggregates.");
  params.addRequiredParam<Real>("x2_aspect",
                                "Effective aspect ratio of pores in porous phase 2 aggregates.");

  params.addRequiredParam<Real>("fluid_K", "Water bulk modulus.");

  params.addRequiredCoupledVar("phi_f", "Fluid volume fraction (porosity).");

  params.addRequiredCoupledVar("phi_x2", "Phase x2 volume fraction.");

  return params;
}

BiphasicLinearViscoElastic::BiphasicLinearViscoElastic(const InputParameters & parameters)
  : Material(parameters),

    // Inputs
    _rho_f(getParam<Real>("rho_f")),
    _rho_x1(getParam<Real>("rho_x1")),
    _rho_x2(getParam<Real>("rho_x2")),
    _k_ref(getParam<Real>("k_ref")),
    _phi_ref(getParam<Real>("phi_ref")),
    _nk(getParam<Real>("nk")),
    _phi_0(getParam<Real>("phi_0")),
    _mu(getParam<Real>("mu")),
    _x1_eta_0(getParam<Real>("x1_eta_0")),
    _x2_eta_0(getParam<Real>("x2_eta_0")),
    _a_eta(getParam<Real>("a_eta")),
    _zeta(getParam<Real>("zeta")),
    _x1_K(getParam<Real>("x1_K")),
    _x1_G(getParam<Real>("x1_G")),
    _x2_K(getParam<Real>("x2_K")),
    _x2_G(getParam<Real>("x2_G")),
    _x1_aspect(getParam<Real>("x1_aspect")),
    _x2_aspect(getParam<Real>("x2_aspect")),
    _fluid_K(getParam<Real>("fluid_K")),

    // Get reference to coupled variables
    _phi_f(adCoupledValue("phi_f")),
    _phi_x2(adCoupledValue("phi_x2")),

    // To calculate

    // Declared properties
    // Densities & volume fractions
    _phi_x1(declareADProperty<Real>("phi_x1")),
    _rho_T(declareADProperty<Real>("rho_T")),

    // Explicitly calculated
    _k(declareADProperty<Real>("k")),

    // Viscous solid flow

    _eta_s(declareADProperty<Real>("eta_s")),
    _eta_compact(declareADProperty<Real>("eta_compact")),

    // Solid elasticity
    _x1_K_d(declareADProperty<Real>("x1_K_d")),
    _x2_K_d(declareADProperty<Real>("x2_K_d")),
    _K_d(declareADProperty<Real>("K_d")),

    _x1_alpha(declareADProperty<Real>("x1_alpha")),
    _x2_alpha(declareADProperty<Real>("x2_alpha")),
    _alpha(declareADProperty<Real>("alpha")),

    _Skempton(declareADProperty<Real>("Skempton")),
    _K_phi(declareADProperty<Real>("K_phi"))

{
}

void
BiphasicLinearViscoElastic::computeQpProperties()
{

  // Determine densities and volume fractions first

  _phi_x1[_qp] = 1.0 - _phi_f[_qp] - _phi_x2[_qp];

  _rho_T[_qp] = (_rho_f * _phi_f[_qp]) + (_rho_x1 * _phi_x1[_qp]) + (_rho_x2 * _phi_x2[_qp]);

  // Fluid flow - permeability
  _k[_qp] = _k_ref * std::pow(_phi_f[_qp] / _phi_ref, _nk);

  // Solid viscosity

  // non-porous effective bulk viscosity -> maybe later we can try different mixing rules but here
  // we just do volume weighted average
  const ADReal eta_s_0 =
      ((_x1_eta_0 * _phi_x1[_qp]) + (_x2_eta_0 * _phi_x2[_qp])) / (1.0 - _phi_f[_qp]);

  // cap the solid viscosity and modify viscosity according to the local porosity

  _eta_s[_qp] = eta_s_0 * std::exp(-_a_eta * (_phi_f[_qp] / _phi_0 - 1.0));

  _eta_compact[_qp] = _zeta * _eta_s[_qp] * (_phi_0 / _phi_f[_qp]);

  // Solid elasticity -- here we use simple volume weighted averaging for solid phase mixes
  // (fluid-solid mixing uses Mori-Tanaka)

  // for Mori-Tanaka
  const ADReal beta_x1 = _x1_G * (3.0 * _x1_K + _x1_G) / (3.0 * _x1_K + 4.0 * _x1_G);

  const ADReal beta_x2 = _x2_G * (3.0 * _x2_K + _x2_G) / (3.0 * _x2_K + 4.0 * _x2_G);

  _x1_K_d[_qp] =
      _x1_K -
      _phi_f[_qp] * (std::pow(_x1_K, 2.0) / (pi * _x1_aspect * beta_x1)) *
          std::pow((1 - _phi_f[_qp] + ((_phi_f[_qp] * _x1_K) / (pi * _x1_aspect * beta_x1))), -1.0);

  _x2_K_d[_qp] =
      _x2_K -
      _phi_f[_qp] * (std::pow(_x2_K, 2.0) / (pi * _x2_aspect * beta_x2)) *
          std::pow((1.0 - _phi_f[_qp] + ((_phi_f[_qp] * _x2_K) / (pi * _x2_aspect * beta_x2))),
                   -1.0);

  _K_d[_qp] = ((_x1_K_d[_qp] * _phi_x1[_qp]) + (_x2_K_d[_qp] * _phi_x2[_qp])) / (1.0 - _phi_f[_qp]);

  _x1_alpha[_qp] = 1.0 - (_x1_K_d[_qp] / _x1_K);
  _x2_alpha[_qp] = 1.0 - (_x2_K_d[_qp] / _x2_K);

  // Rather than finding the averaged fully solid bulk modulus we just average the individual biots
  _alpha[_qp] =
      ((_x1_alpha[_qp] * _phi_x1[_qp]) + (_x2_alpha[_qp] * _phi_x2[_qp])) / (1.0 - _phi_f[_qp]);

  const ADReal K_s = _K_d[_qp] / (1.0 - _alpha[_qp]);

  _Skempton[_qp] =
      ((1.0 / _K_d[_qp]) - (1.0 / K_s)) /
      ((1.0 / _K_d[_qp]) - (1.0 / K_s) + _phi_f[_qp] * ((1.0 / _fluid_K) - (1.0 / K_s)));

  _K_phi[_qp] =
      1.0 / (((1 - _phi_f[_qp]) / _K_d[_qp]) -
             (1 / (((_phi_x1[_qp] * _x1_K) + (_phi_x2[_qp] * _x2_K)) / (1.0 - _phi_f[_qp]))));
}
