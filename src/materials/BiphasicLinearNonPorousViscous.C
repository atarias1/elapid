#include "BiphasicLinearNonPorousViscous.h"
#include "libmesh/utility.h"

registerMooseObject("ElapidApp", BiphasicLinearNonPorousViscous);

InputParameters
BiphasicLinearNonPorousViscous::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<Real>("rho_x1", "Phase 1 constant density.");

  params.addRequiredParam<Real>("rho_x2", "Phase 2 constant density.");

  params.addRequiredParam<Real>("x1_eta", "Phase 1 viscosity at strain rate of 1/s.");

  params.addRequiredParam<Real>("x2_eta", "Phase 2 viscosity at strain rate 1/s.");

  params.addRequiredCoupledVar("phi_x2", "Phase x2 volume fraction.");

  return params;
}

BiphasicLinearNonPorousViscous::BiphasicLinearNonPorousViscous(const InputParameters & parameters)
  : Material(parameters),

    // Inputs
    _rho_x1(getParam<Real>("rho_x1")),
    _rho_x2(getParam<Real>("rho_x2")),
    _x1_eta(getParam<Real>("x1_eta")),
    _x2_eta(getParam<Real>("x2_eta")),

    // Get reference to coupled variables
    _phi_x2(adCoupledValue("phi_x2")),

    // To calculate

    // Declared properties
    // Densities & volume fractions
    _phi_x1(declareADProperty<Real>("phi_x1")),
    _rho_T(declareADProperty<Real>("rho_T")),

    // Viscous solid flow
    _eta_s(declareADProperty<Real>("eta_s"))

{
}

void
BiphasicLinearNonPorousViscous::computeQpProperties()
{

  // Determine densities and volume fractions first

  _phi_x1[_qp] = 1.0 - _phi_x2[_qp];

  _rho_T[_qp] = (_rho_x1 * _phi_x1[_qp]) + (_rho_x2 * _phi_x2[_qp]);

  // Solid viscosity

  // we just do volume weighted average
  _eta_s[_qp] = ((_x1_eta * _phi_x1[_qp]) + (_x2_eta * _phi_x2[_qp]));
}
