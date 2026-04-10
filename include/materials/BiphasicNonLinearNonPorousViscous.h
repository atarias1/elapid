#pragma once
#include "Material.h"

class BiphasicNonLinearNonPorousViscous : public Material
{
public:
  BiphasicNonLinearNonPorousViscous(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

private:
  // Inputs

  // Densities

  const Real & _rho_x1; // Phase 1 constant density
  const Real & _rho_x2; // Phase 2 constant density

  const Real & _max_eta_s; // Maximum viscosity of bulk solid (in case of low strain rate)

  // Solid viscous material properties
  const Real & _x1_eta_0; // Phase 1 viscosity at a strain rate of 1/s
  const Real & _x2_eta_0; // Phase 2 viscosity at a strain rate of 1/s

  const Real & _nsigma_x1; // Stress exponent for phase 1
  const Real & _nsigma_x2; // Stress exponent for phase 2

  // Coupled variables

  // Phase 2 specified here, _phi_x1 calculated implicitly
  const ADVariableValue & _phi_x2;      // vol fraction phase 2
  const ADVariableValue & _v_x;         // x velocity
  const ADVariableValue & _v_y;         // y velocity
  const ADVariableGradient & _grad_v_x; // x velocity gradient
  const ADVariableGradient & _grad_v_y; // y velocity gradient

  // Properties to be computed

  // Density and volume fraction properties
  ADMaterialProperty<Real> & _phi_x1; // Phase x1 volume fraction
  ADMaterialProperty<Real> & _rho_T;  // Total density

  // Solid viscous flow
  ADMaterialProperty<Real> & _eta_x1_eff; // effective non-porous phase x1 viscosity
  ADMaterialProperty<Real> & _eta_x2_eff; // effective non-porous phase x2 viscosity
  ADMaterialProperty<Real> & _eta_s;      // bulk effective viscosity of porous solid
};
