#pragma once
#include "Material.h"

class MetaSerpentinite : public Material
{
public:
  MetaSerpentinite(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

private:
  // Inputs

  // Densities

  const Real & _rho_f;   // Fluid constant density
  const Real & _rho_atg; // Antigorite constant density
  const Real & _rho_ol;  // Olivine constant density

  // Fluid flow material properties
  const Real & _k_ref;     // Reference Permeability
  const Real & _phi_ref;   // Reference Porosity where Ref k was determined
  const Real & _nk;        // powerlaw exponent for k-phi
  const Real & _phi_0;     // background ambient porosity
  const Real & _mu;        // viscosity of water
  const Real & _max_eta_s; // Maximum viscosity of solid (in case of low strain rate)

  // Solid viscous material properties
  // Mixing rule for viscosity
  const Real & _atg_eta_0; // viscosity at a strain rate of 1/s
  const Real & _ol_eta_0;  // for olivine diffusion creep this is also the effective viscosity
  // Porosity effect - we use the exponential rule from Schmalholtz et al., 2023
  const Real & _a_eta; // used to modify the exponetial porosity dependent eta drop

  // Compaction viscosity coef.
  const Real & _zeta;

  // Solid elastic material properties

  const Real & _ol_K;       // ol bulk modulus
  const Real & _ol_G;       // ol shear modulus
  const Real & _atg_K;      // atg bulk modulus
  const Real & _atg_G;      // atg shear modulus
  const Real & _ol_aspect;  // Effective aspect ratio of pores in porous olivine aggregates
  const Real & _atg_aspect; // Effective aspect ratio of pores in porous antigorite aggregates

  const Real & _water_K; // Water bulk modulus

  // Coupled variables
  const ADVariableValue & _phi_f;     // porosity
  const ADVariableValue & _phi_ol;      // vol fraction olivine (phi_ol)
  const ADVariableValue & _v_x;         // x velocity
  const ADVariableValue & _v_y;         // y velocity
  const ADVariableGradient & _grad_v_x; // x velocity gradient
  const ADVariableGradient & _grad_v_y; // y velocity gradient

  // Properties to be computed

  // Density and volume fraction properties
  ADMaterialProperty<Real> & _phi_atg; // atg volume fraction
  ADMaterialProperty<Real> & _rho_T;

  // Fluid flow
  ADMaterialProperty<Real> & _k;

  // Solid viscous flow
  // ADMaterialProperty<Real> & _eta_ol_eff; // for olivine deforming by diffusion creep this is the
  // same as _ol_eta_0
  ADMaterialProperty<Real> & _eta_atg_eff;
  ADMaterialProperty<Real> & _eta_s;
  ADMaterialProperty<Real> & _comp_eta;

  // Solid elasticity
  ADMaterialProperty<Real> & _ol_K_d;    // ol drained modulus
  ADMaterialProperty<Real> & _atg_K_d;   // atg drained modulus
  ADMaterialProperty<Real> & _K_d;       // bulk drained modulus
  ADMaterialProperty<Real> & _ol_alpha;  // ol biot's coef.
  ADMaterialProperty<Real> & _atg_alpha; // atg biot's coef.
  ADMaterialProperty<Real> & _alpha;     // bulk biot's coef.

  ADMaterialProperty<Real> & _Skempton;  // Skempton's coefficient
  ADMaterialProperty<Real> & _K_phi;     // Pore bulk modulus

  //--- Should be not needed ---//
  // // Material variable derivatives
  // ADMaterialProperty<Real> & _drho_s_dP_f;
  // ADMaterialProperty<Real> & _drho_f_dP_f;
  // ADMaterialProperty<Real> & _dX_s_dP_f;

  // // Material time derivatives

  // ADMaterialProperty<Real> & _drho_T_dt;
  // ADMaterialProperty<Real> & _drho_s_X_s_phi_dt;
};
