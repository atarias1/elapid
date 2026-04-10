#pragma once
#include "Material.h"

class BiphasicLinearNonPorousViscous : public Material
{
public:
  BiphasicLinearNonPorousViscous(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

private:
  // Inputs

  // Densities

  const Real & _rho_x1; // Phase 1 constant density
  const Real & _rho_x2; // Phase 2 constant density

  // Solid viscous material properties
  const Real & _x1_eta; // Phase 1 viscosity at a strain rate of 1/s
  const Real & _x2_eta; // Phase 2 viscosity at a strain rate of 1/s

  // Coupled variables

  // Phase 2 specified here, _phi_x1 calculated implicitly
  const ADVariableValue & _phi_x2; // vol fraction phase 2

  // Properties to be computed

  // Density and volume fraction properties
  ADMaterialProperty<Real> & _phi_x1; // Phase x1 volume fraction
  ADMaterialProperty<Real> & _rho_T;  // Total density

  // Solid viscous flow
  ADMaterialProperty<Real> & _eta_s; // bulk effective viscosity of porous solid
};
