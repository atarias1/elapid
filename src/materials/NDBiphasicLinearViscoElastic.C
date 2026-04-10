#include "NDBiphasicLinearViscoElastic.h"
#include "libmesh/utility.h"

registerMooseObject("ElapidApp", NDBiphasicLinearViscoElastic);

InputParameters
NDBiphasicLinearViscoElastic::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<Real>("P_star", "Pressure scale");

  params.addRequiredParam<Real>("rho_f", "Fluid constant density.");

  params.addRequiredParam<Real>("rho_x1", "Phase 1 constant density.");

  params.addRequiredParam<Real>("rho_x2", "Phase 2 constant density.");

  params.addRequiredParam<Real>(
      "phi_ref", "A reference porosity that correlates with the reference permeability.");

  params.addRequiredParam<Real>("x1_eta_0", "Phase 1 constant viscosity.");

  params.addRequiredParam<Real>("x2_eta_0", "Phase 2 constant viscosity.");

  params.addRequiredParam<Real>("a_eta",
                                "Coefficient to modify the porosity effects on solid viscosity.");

  params.addRequiredParam<Real>("x1_K", "Phase 1 bulk modulus.");
  params.addRequiredParam<Real>("x2_K", "Phase 2 bulk modulus.");
  params.addRequiredParam<Real>("x1_G", "Phase 1 shear modulus.");
  params.addRequiredParam<Real>("x2_G", "Phase 2 shear modulus.");

  params.addRequiredParam<Real>("x1_aspect",
                                "Effective aspect ratio of pores in porous phase 1 aggregates.");
  params.addRequiredParam<Real>("x2_aspect",
                                "Effective aspect ratio of pores in porous phase 2 aggregates.");

  params.addRequiredParam<Real>("fluid_K", "Water bulk modulus.");

  params.addRequiredParam<Real>("x1_K_dr", "K_d scaling");
  params.addRequiredParam<Real>("x2_K_dr", "K_d scaling");
  params.addRequiredParam<Real>("x1_alpha_r", "alpha scaling");
  params.addRequiredParam<Real>("x2_alpha_r", "alpha scaling");
  params.addRequiredParam<Real>("x1_K_phir", "K_phi scaling");
  params.addRequiredParam<Real>("x2_K_phir", "K_phi scaling");
  params.addRequiredParam<Real>("x1_B_r", "Skempton scaling");
  params.addRequiredParam<Real>("x2_B_r", "Skempton scaling");

  params.addRequiredCoupledVar("phi_f", "Fluid volume fraction (porosity).");

  params.addRequiredCoupledVar("phi_x2", "Phase x2 volume fraction.");

  return params;
}

NDBiphasicLinearViscoElastic::NDBiphasicLinearViscoElastic(const InputParameters & parameters)
  : Material(parameters),

    // Inputs
    _P_star(getParam<Real>("P_star")),
    _rho_f(getParam<Real>("rho_f")),
    _rho_x1(getParam<Real>("rho_x1")),
    _rho_x2(getParam<Real>("rho_x2")),
    _phi_ref(getParam<Real>("phi_ref")),
    _x1_eta_0(getParam<Real>("x1_eta_0")),
    _x2_eta_0(getParam<Real>("x2_eta_0")),
    _a_eta(getParam<Real>("a_eta")),
    _x1_K(getParam<Real>("x1_K")),
    _x1_G(getParam<Real>("x1_G")),
    _x2_K(getParam<Real>("x2_K")),
    _x2_G(getParam<Real>("x2_G")),
    _x1_aspect(getParam<Real>("x1_aspect")),
    _x2_aspect(getParam<Real>("x2_aspect")),
    _fluid_K(getParam<Real>("fluid_K")),

    _x1_K_dr(getParam<Real>("x1_K_dr")),
    _x1_alpha_r(getParam<Real>("x1_alpha_r")),
    _x1_K_phir(getParam<Real>("x1_K_phir")),
    _x1_B_r(getParam<Real>("x1_B_r")),
    _x2_K_dr(getParam<Real>("x2_K_dr")),
    _x2_alpha_r(getParam<Real>("x2_alpha_r")),
    _x2_K_phir(getParam<Real>("x2_K_phir")),
    _x2_B_r(getParam<Real>("x2_B_r")),

    // Get reference to coupled variables
    _phi_f(adCoupledValue("phi_f")),
    _phi_x2(adCoupledValue("phi_x2")),

    // To calculate

    // Declared properties
    // Densities & volume fractions
    _phi_x1(declareADProperty<Real>("phi_x1")),
    _rho_T(declareADProperty<Real>("rho_T")),
    _rho_s(declareADProperty<Real>("rho_s")),
    _rho_scaled(declareADProperty<Real>("rho_scaled")),

    // Viscous solid flow

    _eta_s(declareADProperty<Real>("eta_s")),

    // Solid elasticity
    _x1_K_d(declareADProperty<Real>("x1_K_d")),
    _x2_K_d(declareADProperty<Real>("x2_K_d")),

    _K_d(declareADProperty<Real>("K_d")),
    _K_d_hat(declareADProperty<Real>("K_d_hat")),

    _x1_alpha(declareADProperty<Real>("x1_alpha")),
    _x2_alpha(declareADProperty<Real>("x2_alpha")),

    _alpha(declareADProperty<Real>("alpha")),
    _alpha_hat(declareADProperty<Real>("alpha_hat")),

    _Skempton_hat(declareADProperty<Real>("Skempton_hat")),
    _K_phi_hat(declareADProperty<Real>("K_phi_hat")),

    _p1(declareADProperty<Real>("p1")),
    _p2(declareADProperty<Real>("p2")),
    _p3(declareADProperty<Real>("p3")),
    _p4(declareADProperty<Real>("p4"))

{
}

void
NDBiphasicLinearViscoElastic::computeQpProperties()
{

  // Determine densities and volume fractions first

  _phi_x1[_qp] = 1.0 - _phi_f[_qp] - _phi_x2[_qp];

  _rho_T[_qp] = (_rho_f * _phi_f[_qp]) + (_rho_x1 * _phi_x1[_qp]) + (_rho_x2 * _phi_x2[_qp]);

  const ADReal x1_svf = (_phi_x1[_qp] / (_phi_x1[_qp] + _phi_x2[_qp]));
  _rho_s[_qp] = (_rho_x1 * x1_svf) + (_rho_x2 * (1 - x1_svf));

  _rho_scaled[_qp] = _rho_T[_qp] / (_rho_s[_qp] - _rho_f);

  // Solid viscosity

  // just get proper volume fraction
  const ADReal ave_eta_0 = (_x1_eta_0 * x1_svf) + (_x2_eta_0 * (1 - x1_svf));

  // Since we scale according to matrix viscosity we need to append the
  // volume averaged background viscosity normalized by the matrix viscosity

  // cap the solid viscosity and modify viscosity according to the local porosity

  _eta_s[_qp] = (ave_eta_0 / _x1_eta_0) * std::exp(-_a_eta * (_phi_f[_qp] / _phi_ref - 1.0));

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

  const ADReal K_dr = (_x1_K_dr * x1_svf) + (_x2_K_dr * (1 - x1_svf));

  _K_d[_qp] = ((_x1_K_d[_qp] * _phi_x1[_qp]) + (_x2_K_d[_qp] * _phi_x2[_qp])) / (1.0 - _phi_f[_qp]);

  _K_d_hat[_qp] = _K_d[_qp] / K_dr;

  _x1_alpha[_qp] = 1.0 - (_x1_K_d[_qp] / _x1_K);
  _x2_alpha[_qp] = 1.0 - (_x2_K_d[_qp] / _x2_K);

  // Rather than finding the averaged fully solid bulk modulus we just average the individual biots
  _alpha[_qp] =
      ((_x1_alpha[_qp] * _phi_x1[_qp]) + (_x2_alpha[_qp] * _phi_x2[_qp])) / (1.0 - _phi_f[_qp]);

  const ADReal alpha_r = (_x1_alpha_r * x1_svf) + (_x2_alpha_r * (1 - x1_svf));

  _alpha_hat[_qp] = _alpha[_qp] / alpha_r;

  const ADReal K_s = _K_d[_qp] / (1.0 - _alpha[_qp]);
  const ADReal B_r = (_x1_B_r * x1_svf) + (_x2_B_r * (1 - x1_svf));

  _Skempton_hat[_qp] =
      (((1.0 / _K_d[_qp]) - (1.0 / K_s)) /
       ((1.0 / _K_d[_qp]) - (1.0 / K_s) + _phi_f[_qp] * ((1.0 / _fluid_K) - (1.0 / K_s)))) /
      B_r;

  const ADReal K_phir = (_x1_K_phir * x1_svf) + (_x2_K_phir * (1 - x1_svf));

  _K_phi_hat[_qp] =
      (1.0 / (((1 - _phi_f[_qp]) / _K_d[_qp]) -
              (1 / (((_phi_x1[_qp] * _x1_K) + (_phi_x2[_qp] * _x2_K)) / (1.0 - _phi_f[_qp]))))) /
      K_phir;

  _p1[_qp] = _P_star / K_dr;
  _p2[_qp] = _P_star * alpha_r / K_dr;
  _p3[_qp] = (_P_star * alpha_r) / (K_dr * B_r);
  _p4[_qp] = _P_star / K_phir;
}
