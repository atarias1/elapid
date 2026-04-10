#pragma once
#include "Material.h"

class NDBiphasicLinearViscoElastic : public Material
{
public:
  NDBiphasicLinearViscoElastic(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

private:
  // Inputs

  const Real & _P_star; // Pressure Scale
  // Densities

  const Real & _rho_f;  // Fluid constant density
  const Real & _rho_x1; // Phase 1 constant density
  const Real & _rho_x2; // Phase 2 constant density

  // Fluid flow material properties
  const Real & _phi_ref; // Reference Porosity where Ref k was determined

  // Solid viscous material properties
  const Real & _x1_eta_0; // Phase 1 viscosity at a strain rate of 1/s
  const Real & _x2_eta_0; // Phase 2 viscosity at a strain rate of 1/s

  // Porosity effect - we use the exponential rule from Schmalholtz et al., 2023
  const Real & _a_eta; // used to modify the exponetial porosity dependent eta drop

  // Solid elastic material properties

  const Real & _x1_K;      // Phase 1 bulk modulus
  const Real & _x1_G;      // Phase 1 shear modulus
  const Real & _x2_K;      // atg bulk modulus
  const Real & _x2_G;      // atg shear modulus
  const Real & _x1_aspect; // Effective aspect ratio of pores in porous phase 1 aggregates
  const Real & _x2_aspect; // Effective aspect ratio of pores in porous phase 2 aggregates

  const Real & _fluid_K; // Fluid bulk modulus

  const Real & _x1_K_dr;
  const Real & _x2_K_dr;
  const Real & _x1_alpha_r;
  const Real & _x2_alpha_r;
  const Real & _x1_K_phir;
  const Real & _x2_K_phir;
  const Real & _x1_B_r;
  const Real & _x2_B_r;

  // Coupled variables
  const ADVariableValue & _phi_f; // porosity

  // Phase 2 specified here, _phi_x1 calculated implicitly -x1 is matrix
  const ADVariableValue & _phi_x2; // vol fraction phase 2

  // Properties to be computed

  // Density and volume fraction properties
  ADMaterialProperty<Real> & _phi_x1;     // Phase x1 volume fraction
  ADMaterialProperty<Real> & _rho_T;      // Total density
  ADMaterialProperty<Real> & _rho_s;      // Solid density
  ADMaterialProperty<Real> & _rho_scaled; // Density for passing to scaled gravity

  // Fluid flow
  // We keep it simple for now and assume both phase permeability is similar. No need to compute.
  // ADMaterialProperty<Real> & _k; // Local permeability

  // Solid viscous flow
  ADMaterialProperty<Real> & _eta_s; // bulk effective viscosity of porous solid

  // Solid elasticity
  ADMaterialProperty<Real> & _x1_K_d;  // phase x1 drained modulus
  ADMaterialProperty<Real> & _x2_K_d;  // phase x2 drained modulus
  ADMaterialProperty<Real> & _K_d;     // bulk drained modulus
  ADMaterialProperty<Real> & _K_d_hat; // drained bulk modulus

  ADMaterialProperty<Real> & _x1_alpha;  // x1 biot's coef.
  ADMaterialProperty<Real> & _x2_alpha;  // x2 biot's coef.
  ADMaterialProperty<Real> & _alpha;     // biot's coef.
  ADMaterialProperty<Real> & _alpha_hat; // biot's coef. scaled

  ADMaterialProperty<Real> & _Skempton_hat; // Skempton's coefficient scaled
  ADMaterialProperty<Real> & _K_phi_hat;    // Pore bulk modulus scaled

  // pressure scales - may remove later if performance isn't improved
  ADMaterialProperty<Real> & _p1;
  ADMaterialProperty<Real> & _p2;
  ADMaterialProperty<Real> & _p3;
  ADMaterialProperty<Real> & _p4;
};
