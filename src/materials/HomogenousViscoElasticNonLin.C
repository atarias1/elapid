#include "HomogenousViscoElasticNonLin.h"
#include "libmesh/utility.h"

registerMooseObject("ElapidApp", HomogenousViscoElasticNonLin);

InputParameters
HomogenousViscoElasticNonLin::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<Real>("rho_f", "Fluid constant density (water).");

  params.addRequiredParam<Real>("rho_s", "Solid constant density.");

  params.addRequiredParam<Real>(
      "k_ref", "A reference permeability value that correlates to a reference porosity.");

  params.addRequiredParam<Real>(
      "phi_ref", "A reference porosity that correlates to a reference permeability.");

  params.addRequiredParam<Real>("nk", "The Kozemy-Carman exponent.");

  params.addRequiredParam<Real>("phi_0", "The initial ambient porosity.");

  params.addRequiredParam<Real>("mu", "Viscosity of water.");

  params.addRequiredParam<Real>("max_eta_s", "Maximum viscosity in case of low strain rate");

  params.addRequiredParam<Real>("eta_0", "Viscosity of solid at strain rate of 1/s.");
  params.addRequiredParam<Real>("nsigma", "Stress Exponent -> 1 for linear.");

  params.addRequiredParam<Real>("a_eta",
                                "Coefficient to modify the porosity effects on solid viscosity.");

  params.addRequiredParam<Real>("zeta", "Coefficient to control the compaction viscosity.");

  params.addRequiredParam<Real>("K", "Solid bulk modulus.");
  params.addRequiredParam<Real>("G", "Solid shear modulus.");

  params.addRequiredParam<Real>("aspect",
                                "Effective aspect ratio of pores in porous aggregates.");

  params.addRequiredParam<Real>("water_K", "Water bulk modulus.");

  params.addRequiredCoupledVar("phi_f", "Fluid pressure used to determine thermodynamic values.");

  params.addRequiredCoupledVar("v_x", "X velocity to find strain rate.");

  params.addRequiredCoupledVar("v_y", "Y velocity to find strain rate.");

  return params;
}

HomogenousViscoElasticNonLin::HomogenousViscoElasticNonLin(const InputParameters & parameters)
  : Material(parameters),

    // Inputs
    _rho_f(getParam<Real>("rho_f")),
    _rho_s(getParam<Real>("rho_s")),
    _k_ref(getParam<Real>("k_ref")),
    _phi_ref(getParam<Real>("phi_ref")),
    _nk(getParam<Real>("nk")),
    _phi_0(getParam<Real>("phi_0")),
    _mu(getParam<Real>("mu")),
    _max_eta_s(getParam<Real>("max_eta_s")),
    _eta_0(getParam<Real>("eta_0")),
    _nsigma(getParam<Real>("nsigma")),
    _a_eta(getParam<Real>("a_eta")),
    _zeta(getParam<Real>("zeta")),
    _K(getParam<Real>("K")),
    _G(getParam<Real>("G")),
    _aspect(getParam<Real>("aspect")),
    _water_K(getParam<Real>("water_K")),

    // Get reference to coupled variables
    _phi_f(adCoupledValue("phi_f")),
    _v_x(adCoupledValue("v_x")),
    _v_y(adCoupledValue("v_y")),
    _grad_v_x(adCoupledGradient("v_x")),
    _grad_v_y(adCoupledGradient("v_y")),

    // To calculate

    // Declared properties
    // Densities & volume fractions
    _phi_s(declareADProperty<Real>("phi_s")),
    _rho_T(declareADProperty<Real>("rho_T")),

    // Explicitly calculated
    _k(declareADProperty<Real>("k")),

    // Viscous solid flow
    _eta_s_eff(declareADProperty<Real>("eta_s_eff")),
    _eta_s(declareADProperty<Real>("eta_s")),
    _comp_eta(declareADProperty<Real>("comp_eta")),

    // Solid elasticity
    _K_d(declareADProperty<Real>("K_d")),
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
HomogenousViscoElasticNonLin::computeQpProperties()
{

  // Determine densities and volume fractions first

  _phi_s[_qp] = 1.0 - _phi_f[_qp];

  _rho_T[_qp] = (_rho_f * _phi_f[_qp]) + (_rho_s * _phi_s[_qp]);

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
  _eta_s_eff[_qp] = _eta_0 * std::pow(eff_strain_rate, (1.0 / _nsigma) - 1.0); // For more general phase need to include an input for stress exponent

  // quick test if linear viscosity fixes issue

  // cap the solid viscosity

  if (_eta_s_eff[_qp] * std::exp(-_a_eta * (_phi_f[_qp] / _phi_0 - 1.0)) <  // olivine version  _ol_eta_0 * std::exp(-_a_eta * (_phi_f[_qp] / _phi_0 - 1.0))
      _max_eta_s)
  {
    _eta_s[_qp] = _eta_s_eff[_qp] * std::exp(-_a_eta * (_phi_f[_qp] / _phi_0 - 1.0));
  }
  else
  { // Here we use _max_eta_s if the predicted viscosity gets too high
    _eta_s[_qp] = _max_eta_s;
  }

  _comp_eta[_qp] = _zeta * _eta_s[_qp] * (_phi_0 / _phi_f[_qp]);

  // Solid elasticity

  // for Mori-Tanaka
  const ADReal beta_s = _G * (3.0 * _K + _G) / (3.0 * _K + 4.0 * _G);
  
  _K_d[_qp] =
      _K -
      _phi_f[_qp] * (std::pow(_K, 2.0) / (pi * _aspect * beta_s)) *
          std::pow((1 - _phi_f[_qp] + ((_phi_f[_qp] * _K) / (pi * _aspect * beta_s))), -1.0);

  _alpha[_qp] = 1.0 - (_K_d[_qp] / _K);

  _Skempton[_qp] = ((1.0 / _K_d[_qp]) - (1.0 / _K)) /
                   ((1.0 / _K_d[_qp]) - (1.0 / _K) + _phi_f[_qp] * ((1.0 / _water_K) - (1.0 / _K)));

  _K_phi[_qp] = 1.0 / (
    ((1-_phi_f[_qp])/_K_d[_qp]) - (
      1/_K
    )
  );

}
