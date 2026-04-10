#include "NDSinglePhaseLinearViscoElastic.h"
#include "libmesh/utility.h"

registerMooseObject("ElapidApp", NDSinglePhaseLinearViscoElastic);

InputParameters
NDSinglePhaseLinearViscoElastic::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<Real>("P_star", "Pressure scale");

  params.addRequiredParam<Real>("rho_f", "Fluid constant density.");

  params.addRequiredParam<Real>("rho_s", "Solid constant density.");

  params.addRequiredParam<Real>(
      "phi_ref", "A reference porosity that correlates with the reference permeability.");

  params.addRequiredParam<Real>("eta_s_0", "Constant nonporous solid viscosity.");

  params.addRequiredParam<Real>(
      "a_eta", "Coefficient to modify the effects of porosity on solid viscosity.");

  params.addRequiredParam<Real>("K", "Solid bulk modulus.");
  params.addRequiredParam<Real>("G", "Solid shear modulus.");

  params.addRequiredParam<Real>("aspect_ratio",
                                "Effective aspect ratio of pores in porous aggregates.");

  params.addRequiredParam<Real>("fluid_K", "Fluid bulk modulus.");

  params.addRequiredParam<Real>("K_dr", "K_d scaling");
  params.addRequiredParam<Real>("alpha_r", "alpha scaling");
  params.addRequiredParam<Real>("K_phir", "K_phi scaling");
  params.addRequiredParam<Real>("B_r", "Skempton scaling");

  params.addRequiredCoupledVar("phi_f", "Fluid pressure used to determine thermodynamic values.");

  return params;
}

NDSinglePhaseLinearViscoElastic::NDSinglePhaseLinearViscoElastic(const InputParameters & parameters)
  : Material(parameters),

    // Inputs
    _P_star(getParam<Real>("P_star")),
    _rho_f(getParam<Real>("rho_f")),
    _rho_s(getParam<Real>("rho_s")),
    _phi_ref(getParam<Real>("phi_ref")),
    _a_eta(getParam<Real>("a_eta")),
    _K(getParam<Real>("K")),
    _G(getParam<Real>("G")),
    _aspect_ratio(getParam<Real>("aspect_ratio")),
    _fluid_K(getParam<Real>("fluid_K")),

    _K_dr(getParam<Real>("K_dr")),
    _alpha_r(getParam<Real>("alpha_r")),
    _K_phir(getParam<Real>("K_phir")),
    _B_r(getParam<Real>("B_r")),

    // Get reference to coupled variables
    _phi_f(adCoupledValue("phi_f")),

    // To calculate

    // Declared properties
    // Densities & volume fractions
    _phi_s(declareADProperty<Real>("phi_s")),
    _rho_T(declareADProperty<Real>("rho_T")),
    _rho_scaled(declareADProperty<Real>("rho_scaled")),

    // Explicitly calculated

    // Viscous solid flow
    _eta_s(declareADProperty<Real>("eta_s")),

    // Solid elasticity
    _K_d(declareADProperty<Real>("K_d")),
    _K_d_hat(declareADProperty<Real>("K_d_hat")),
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
NDSinglePhaseLinearViscoElastic::computeQpProperties()
{

  // Determine densities and volume fractions

  _phi_s[_qp] = 1.0 - _phi_f[_qp];

  _rho_T[_qp] = (_rho_f * _phi_f[_qp]) + (_rho_s * _phi_s[_qp]);

  _rho_scaled[_qp] = _rho_T[_qp] / (_rho_s - _rho_f);

  // Porous solid viscosity

  _eta_s[_qp] = std::exp(-_a_eta * (_phi_f[_qp] / _phi_ref - 1.0));

  // Poro-elasticity

  // for Mori-Tanaka
  const ADReal beta_s = _G * (3.0 * _K + _G) / (3.0 * _K + 4.0 * _G);

  // Probably not the best way to handle poroelastic vars - adjust later
  _K_d[_qp] =
      (_K - _phi_f[_qp] * (std::pow(_K, 2.0) / (pi * _aspect_ratio * beta_s)) *
                std::pow((1 - _phi_f[_qp] + ((_phi_f[_qp] * _K) / (pi * _aspect_ratio * beta_s))),
                         -1.0));

  _K_d_hat[_qp] = _K_d[_qp] / _K_dr;

  _alpha_hat[_qp] = (1.0 - (_K_d[_qp] / _K)) / _alpha_r;

  _Skempton_hat[_qp] =
      (((1.0 / _K_d[_qp]) - (1.0 / _K)) /
       ((1.0 / _K_d[_qp]) - (1.0 / _K) + _phi_f[_qp] * ((1.0 / _fluid_K) - (1.0 / _K)))) /
      _B_r;

  _K_phi_hat[_qp] = (1.0 / (((1 - _phi_f[_qp]) / _K_d[_qp]) - (1 / _K))) / _K_phir;

  _p1[_qp] = _P_star / _K_dr;
  _p2[_qp] = _P_star * _alpha_r / _K_dr;
  _p3[_qp] = (_P_star * _alpha_r) / (_K_dr * _B_r);
  _p4[_qp] = _P_star / _K_phir;
}
