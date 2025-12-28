#pragma once
#include "Material.h"

class SinglePhaseNonLinearViscoElastic : public Material
{
public:
  SinglePhaseNonLinearViscoElastic(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

private:
  // Inputs

  // Densities

  const Real & _rho_f; // Fluid constant density
  const Real & _rho_s; // Solid constant density

  // Fluid flow material properties
  const Real & _k_ref;     // Reference Permeability
  const Real & _phi_ref;   // Reference Porosity where Ref k was determined
  const Real & _nk;        // powerlaw exponent for k-phi
  const Real & _phi_0;     // background porosity
  const Real & _mu;        // viscosity of fluid phase
  const Real & _max_eta_s; // Maximum viscosity of solid (in case of low strain rate)

  // Solid viscous material properties
  // Mixing rule for viscosity
  const Real & _eta_0;  // viscosity at a strain rate of 1/s
  const Real & _nsigma; // Stress exponent

  // Porosity effect - we use the exponential rule from Schmalholtz et al., 2023 - can modify
  const Real & _a_eta; // used to modify the exponetial porosity dependent eta drop

  // Compaction viscosity coef.
  const Real & _zeta;

  // Solid elastic material properties

  const Real & _K;            // bulk modulus
  const Real & _G;            // shear modulus
  const Real & _aspect_ratio; // Effective aspect ratio of pores

  const Real & _fluid_K; // fluid bulk modulus

  // Coupled variables
  const ADVariableValue & _phi_f;       // porosity
  const ADVariableValue & _v_x;         // x velocity
  const ADVariableValue & _v_y;         // y velocity
  const ADVariableGradient & _grad_v_x; // x velocity gradient
  const ADVariableGradient & _grad_v_y; // y velocity gradient

  // Properties to be computed

  // Density and volume fraction properties
  ADMaterialProperty<Real> & _phi_s; // solid volume fraction
  ADMaterialProperty<Real> & _rho_T; // Total density

  // Fluid flow
  ADMaterialProperty<Real> & _k; // permeability

  // Solid viscous flow

  ADMaterialProperty<Real> & _eta_s_eff;   // effective non-porous solid viscosity
  ADMaterialProperty<Real> & _eta_s;       // porous solid viscosity
  ADMaterialProperty<Real> & _eta_compact; // volumetric viscosity

  // Solid elasticity
  ADMaterialProperty<Real> & _K_d;   // drained bulk modulus
  ADMaterialProperty<Real> & _alpha; // biot's coef.

  ADMaterialProperty<Real> & _Skempton; // Skempton's coefficient
  ADMaterialProperty<Real> & _K_phi;    // Pore bulk modulus
};
