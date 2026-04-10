#pragma once
#include "Material.h"

class NDSinglePhaseLinearViscoElastic : public Material
{
public:
  NDSinglePhaseLinearViscoElastic(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

private:
  // Inputs
  const Real & _P_star; // Pressure Scale

  // Densities

  const Real & _rho_f; // Fluid constant density
  const Real & _rho_s; // Solid constant density

  // Fluid flow material properties
  const Real & _phi_ref; // Reference Porosity where Ref k was determined

  // Porosity effect - we use the exponential rule from Schmalholtz et al., 2023 - can modify
  const Real & _a_eta; // used to modify the exponetial porosity dependent eta drop

  // Solid elastic material properties

  const Real & _K;            // bulk modulus
  const Real & _G;            // shear modulus
  const Real & _aspect_ratio; // Effective aspect ratio of pores

  const Real & _fluid_K; // fluid bulk modulus

  // Scales - may modify this later maybe unecessary
  // Values at X(phi_r)
  const Real & _K_dr;
  const Real & _alpha_r;
  const Real & _K_phir;
  const Real & _B_r;

  // Coupled variables
  const ADVariableValue & _phi_f; // porosity

  // Properties to be computed

  // Density and volume fraction properties
  ADMaterialProperty<Real> & _phi_s;      // solid volume fraction
  ADMaterialProperty<Real> & _rho_T;      // Total density
  ADMaterialProperty<Real> & _rho_scaled; // Density for passing to scaled gravity

  // Solid viscous flow

  ADMaterialProperty<Real> & _eta_s; // non-dim porous solid viscosity

  // Solid elasticity - non-dim form - probably need to rethink this
  ADMaterialProperty<Real> & _K_d;       // dim form drained bulk modulus
  ADMaterialProperty<Real> & _K_d_hat;   // drained bulk modulus
  ADMaterialProperty<Real> & _alpha_hat; // biot's coef.

  ADMaterialProperty<Real> & _Skempton_hat; // Skempton's coefficient
  ADMaterialProperty<Real> & _K_phi_hat;    // Pore bulk modulus

  // pressure scales - may remove later if performance isn't improved
  ADMaterialProperty<Real> & _p1;
  ADMaterialProperty<Real> & _p2;
  ADMaterialProperty<Real> & _p3;
  ADMaterialProperty<Real> & _p4;
};
