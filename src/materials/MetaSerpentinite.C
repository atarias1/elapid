#include "MetaSerpentinite.h"
#include "libmesh/utility.h"

registerMooseObject("ElapidApp", MetaSerpentinite);

InputParameters
MetaSerpentinite::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<Real>("rho_f", "Fluid constant density (water).");

  params.addRequiredParam<Real>("rho_atg", "Antigorite constant density.");

  params.addRequiredParam<Real>("rho_ol", "Olivine constant density.");

  params.addRequiredParam<Real>(
      "k_ref", "A reference permeability value that correlates to a reference porosity.");

  params.addRequiredParam<Real>(
      "phi_ref", "A reference porosity that correlates to a reference permeability.");

  params.addRequiredParam<Real>("nk", "The Kozemy-Carman exponent.");

  params.addRequiredParam<Real>("phi_0", "The initial ambient porosity.");

  params.addRequiredParam<Real>("mu", "Viscosity of water.");

  params.addRequiredParam<Real>("max_eta_s", "Maximum viscosity in case of low strain rate");

  params.addRequiredParam<Real>("atg_eta_0", "Viscosity of antigorite at strain rate of 1/s.");

  params.addRequiredParam<Real>(
      "ol_eta_0", "Viscosity of olivine at strain rate 1/s (diffusion creep here is constant).");

  params.addRequiredParam<Real>("a_eta",
                                "Coefficient to modify the porosity effects on solid viscosity.");

  params.addRequiredParam<Real>("zeta", "Coefficient to control the compaction viscosity.");

  params.addRequiredParam<Real>("ol_K", "Olivine bulk modulus.");
  params.addRequiredParam<Real>("atg_K", "Antigorite bulk modulus.");
  params.addRequiredParam<Real>("ol_G", "Olivine shear modulus.");
  params.addRequiredParam<Real>("atg_G", "Antigorite shear modulus.");

  params.addRequiredParam<Real>("ol_aspect",
                                "Effective aspect ratio of pores in porous olivine aggregates.");
  params.addRequiredParam<Real>("atg_aspect",
                                "Effective aspect ratio of pores in porous antigorite aggregates.");

  params.addRequiredParam<Real>("water_K", "Water bulk modulus.");

  params.addRequiredCoupledVar("phi_f", "Fluid pressure used to determine thermodynamic values.");

  params.addRequiredCoupledVar("phi_ol", "Fluid pressure used to determine thermodynamic values.");

  params.addRequiredCoupledVar("v_x", "X velocity to find strain rate.");

  params.addRequiredCoupledVar("v_y", "Y velocity to find strain rate.");

  return params;
}

MetaSerpentinite::MetaSerpentinite(const InputParameters & parameters)
  : Material(parameters),

    // Inputs
    _rho_f(getParam<Real>("rho_f")),
    _rho_atg(getParam<Real>("rho_atg")),
    _rho_ol(getParam<Real>("rho_ol")),
    _k_ref(getParam<Real>("k_ref")),
    _phi_ref(getParam<Real>("phi_ref")),
    _nk(getParam<Real>("nk")),
    _phi_0(getParam<Real>("phi_0")),
    _mu(getParam<Real>("mu")),
    _max_eta_s(getParam<Real>("max_eta_s")),
    _atg_eta_0(getParam<Real>("atg_eta_0")),
    _ol_eta_0(getParam<Real>("ol_eta_0")),
    _a_eta(getParam<Real>("a_eta")),
    _zeta(getParam<Real>("zeta")),
    _ol_K(getParam<Real>("ol_K")),
    _ol_G(getParam<Real>("ol_G")),
    _atg_K(getParam<Real>("atg_K")),
    _atg_G(getParam<Real>("atg_G")),
    _ol_aspect(getParam<Real>("ol_aspect")),
    _atg_aspect(getParam<Real>("atg_aspect")),
    _water_K(getParam<Real>("water_K")),

    // Get reference to coupled variables
    _phi_f(adCoupledValue("phi_f")),
    _phi_ol(adCoupledValue("phi_ol")),
    _v_x(adCoupledValue("v_x")),
    _v_y(adCoupledValue("v_y")),
    _grad_v_x(adCoupledGradient("v_x")),
    _grad_v_y(adCoupledGradient("v_y")),

    // To calculate

    // Declared properties
    // Densities & volume fractions
    _phi_atg(declareADProperty<Real>("phi_atg")),
    _rho_T(declareADProperty<Real>("rho_T")),

    // Explicitly calculated
    _k(declareADProperty<Real>("k")),

    // Viscous solid flow
    _eta_atg_eff(declareADProperty<Real>("eta_atg_eff")),
    _eta_s(declareADProperty<Real>("eta_s")),
    _comp_eta(declareADProperty<Real>("comp_eta")),

    // Solid elasticity
    _ol_K_d(declareADProperty<Real>("ol_K_d")),
    _atg_K_d(declareADProperty<Real>("atg_K_d")),
    _K_d(declareADProperty<Real>("K_d")),
    _ol_alpha(declareADProperty<Real>("ol_alpha")),
    _atg_alpha(declareADProperty<Real>("atg_alpha")),
    _alpha(declareADProperty<Real>("alpha")),
    _Skempton(declareADProperty<Real>("Skempton")),
    _K_phi(declareADProperty<Real>("K_phi"))

// // Variable derivatives - un-needed it seems
// _drho_s_dP_f(declareADProperty<Real>("drho_s_dP_f")),
// _drho_f_dP_f(declareADProperty<Real>("drho_f_dP_f")),
// _dX_s_dP_f(declareADProperty<Real>("dX_s_dP_f")),

// Material time derivatives - maybe un-needed
// _drhophi_f_dt(declareADProperty<Real>("drhophi_f_dt")),

{
  // Later here we will set up the functions to determing rho_s rho_f and X_s
}

void
MetaSerpentinite::computeQpProperties()
{

  // Determine densities and volume fractions first

  _phi_atg[_qp] = 1.0 - _phi_f[_qp] - _phi_ol[_qp];

  _rho_T[_qp] = (_rho_f * _phi_f[_qp]) + (_rho_atg * _phi_atg[_qp]) + (_rho_ol * _phi_ol[_qp]);

  // Fluid flow
  _k[_qp] = _k_ref * std::pow(_phi_f[_qp] / _phi_ref, _nk);

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

  // n for atg flow law = 3.8
  _eta_atg_eff[_qp] = _atg_eta_0 * std::pow(eff_strain_rate, (1.0 / 3.8) - 1.0);

  const ADReal eta_s_0 =
      ((_eta_atg_eff[_qp] * _phi_atg[_qp]) + (_ol_eta_0 * _phi_ol[_qp])) / (1.0 - _phi_f[_qp]);

  // quick test if linear viscosity fixes issue

  // cap the solid viscosity

  if (eta_s_0 * std::exp(-_a_eta * (_phi_f[_qp] / _phi_0 - 1.0)) <  // olivine version  _ol_eta_0 * std::exp(-_a_eta * (_phi_f[_qp] / _phi_0 - 1.0))
      _max_eta_s)
  {
    _eta_s[_qp] = eta_s_0 * std::exp(-_a_eta * (_phi_f[_qp] / _phi_0 - 1.0));
  }
  else
  { // Here we use _max_eta_s if the predicted viscosity gets too high
    _eta_s[_qp] = _max_eta_s;
  }

  _comp_eta[_qp] = _zeta * _eta_s[_qp] * (_phi_0 / _phi_f[_qp]);

  // Solid elasticity

  // for Mori-Tanaka
  const ADReal beta_ol = _ol_G * (3.0 * _ol_K + _ol_G) / (3.0 * _ol_K + 4.0 * _ol_G);

  const ADReal beta_atg = _atg_G * (3.0 * _atg_K + _atg_G) / (3.0 * _atg_K + 4.0 * _atg_G);

  // We always calculate this even if there is no olivine at the qp, we accomadate this by simple
  // averaging If we really want to plot this later we can post process it using ol vol frac
  _ol_K_d[_qp] =
      _ol_K -
      _phi_f[_qp] * (std::pow(_ol_K, 2.0) / (pi * _ol_aspect * beta_ol)) *
          std::pow((1 - _phi_f[_qp] + ((_phi_f[_qp] * _ol_K) / (pi * _ol_aspect * beta_ol))), -1.0);

  _atg_K_d[_qp] =
      _atg_K -
      _phi_f[_qp] * (std::pow(_atg_K, 2.0) / (pi * _atg_aspect * beta_atg)) *
          std::pow((1.0 - _phi_f[_qp] + ((_phi_f[_qp] * _atg_K) / (pi * _atg_aspect * beta_atg))),
                   -1.0);

  _K_d[_qp] =
      ((_atg_K_d[_qp] * _phi_atg[_qp]) + (_ol_K_d[_qp] * _phi_ol[_qp])) / (1.0 - _phi_f[_qp]);

  _ol_alpha[_qp] = 1.0 - (_ol_K_d[_qp] / _ol_K);
  _atg_alpha[_qp] = 1.0 - (_atg_K_d[_qp] / _atg_K);

  // Rather than finding the averaged fully solid bulk modulus we just average the individual biots
  // I think it is fine
  _alpha[_qp] =
      ((_atg_alpha[_qp] * _phi_atg[_qp]) + (_ol_alpha[_qp] * _phi_ol[_qp])) / (1.0 - _phi_f[_qp]);

  const ADReal K_s = _K_d[_qp] / (1.0 - _alpha[_qp]);

  _Skempton[_qp] = ((1.0 / _K_d[_qp]) - (1.0 / K_s)) /
                   ((1.0 / _K_d[_qp]) - (1.0 / K_s) + _phi_f[_qp] * ((1.0 / _water_K) - (1.0 / K_s)));

  _K_phi[_qp] = 1.0 / (
    ((1-_phi_f[_qp])/_K_d[_qp]) - (
      1/(((_phi_ol[_qp] * _ol_K)+ (_phi_atg[_qp] * _atg_K))/(1.0 - _phi_f[_qp]))
    )
  );

}
