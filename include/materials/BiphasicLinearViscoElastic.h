#pragma once
#include "Material.h"

class BiphasicLinearViscoElastic : public Material
{
public:
  BiphasicLinearViscoElastic(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

private:
  // Inputs

  // Densities

  const Real & _rho_f;  // Fluid constant density
  const Real & _rho_x1; // Phase 1 constant density
  const Real & _rho_x2; // Phase 2 constant density

  // Fluid flow material properties
  const Real & _k_ref;   // Reference Permeability
  const Real & _phi_ref; // Reference Porosity where Ref k was determined
  const Real & _nk;      // powerlaw exponent for k-phi
  const Real & _phi_0;   // background ambient porosity
  const Real & _mu;      // viscosity of fluid phase

  // Solid viscous material properties
  const Real & _x1_eta_0; // Phase 1 viscosity at a strain rate of 1/s
  const Real & _x2_eta_0; // Phase 2 viscosity at a strain rate of 1/s

  // Porosity effect - we use the exponential rule from Schmalholtz et al., 2023
  const Real & _a_eta; // used to modify the exponetial porosity dependent eta drop

  // Compaction viscosity coef.
  const Real & _zeta;

  // Solid elastic material properties

  const Real & _x1_K;      // Phase 1 bulk modulus
  const Real & _x1_G;      // Phase 1 shear modulus
  const Real & _x2_K;      // atg bulk modulus
  const Real & _x2_G;      // atg shear modulus
  const Real & _x1_aspect; // Effective aspect ratio of pores in porous phase 1 aggregates
  const Real & _x2_aspect; // Effective aspect ratio of pores in porous phase 2 aggregates

  const Real & _fluid_K; // Fluid bulk modulus

  // Coupled variables
  const ADVariableValue & _phi_f; // porosity

  // Phase 2 specified here, _phi_x1 calculated implicitly
  const ADVariableValue & _phi_x2; // vol fraction phase 2

  // Properties to be computed

  // Density and volume fraction properties
  ADMaterialProperty<Real> & _phi_x1; // Phase x1 volume fraction
  ADMaterialProperty<Real> & _rho_T;  // Total density

  // Fluid flow
  ADMaterialProperty<Real> & _k; // Local permeability

  // Solid viscous flow
  ADMaterialProperty<Real> & _eta_s;       // bulk effective viscosity of porous solid
  ADMaterialProperty<Real> & _eta_compact; // volumetric viscosity

  // Solid elasticity
  ADMaterialProperty<Real> & _x1_K_d; // phase x1 drained modulus
  ADMaterialProperty<Real> & _x2_K_d; // phase x2 drained modulus
  ADMaterialProperty<Real> & _K_d;    // bulk drained modulus

  ADMaterialProperty<Real> & _x1_alpha; // x1 biot's coef.
  ADMaterialProperty<Real> & _x2_alpha; // x2 biot's coef.
  ADMaterialProperty<Real> & _alpha;    // bulk biot's coef.

  ADMaterialProperty<Real> & _Skempton; // Skempton's coefficient
  ADMaterialProperty<Real> & _K_phi;    // Pore bulk modulus
};
