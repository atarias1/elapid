# Elapid
Elapid is a MOOSE application for multiphysics finite element modeling of porous media.
Elapid provides custom physics kernels and material models for the (non-)linear deformation of heterogenous porous materials with Maxwell-type viscoelastic rheologies.

## Relevant equations

Currently Elapid provides physics kernels to model Maxwell-type viscoelastic rheology. Derivations of the relevant equations shown here inherit from the work of McKenzie, 1984, regarding the compaction of partially molten rocks. The forms here are implemented according to the derivation in Yarushina and Podlachikov, 2015. Variables and parameters listed at the bottom of the page.

**Solid Deformation ($v_s$)**

$$\nabla  \cdot  \sigma + \rho_T g = 0$$

**Solid Compressibility ($P$)**

$$\nabla \cdot v_s + \frac{1}{K_d}\frac{dP}{dt}-\frac{\alpha}{K_d}\frac{dP_f}{dt} + \frac{P-P_f}{(1-\phi_f)\eta_c} = 0$$

**Fluid Pressure Evolution ($P_f$)**

$$\frac{\alpha}{K_d B}\frac{dP_f}{dt} - \frac{\alpha}{K_d}\frac{dP}{dt} - \nabla \cdot \left(\frac{k(\phi_f)}{\mu}(\nabla P_f + \rho_f g)\right) - \frac{P-P_f}{(1- \phi_f)\eta_c} = 0$$

**Porosity Evolution ($\phi_f$)**

$$\frac{\partial \phi_f}{\partial t} + \frac{1-\phi_f}{K_\phi}\frac{d P}{d t} - \frac{1 - \phi_f}{K_\phi}\frac{d P_f}{d t} - (1 - \phi_f)\frac{P_f - P}{\eta_c} = 0$$

**Solid Conservation Equations ($\phi_i$)**

$$\frac{\partial \phi_i}{\partial t} + \nabla \cdot (v_s \phi_i) = 0$$

## Equivalent Physics Kernels

In some cases, MOOSE default kernels are sufficient, however, we have included additional custom physics kernels to model the above relevant system of equations. Unlisted terms in the above equations are included using default MOOSE kernels.

### Provided Kernels
| Term  | Kernel Name |
| ------------- | ------------- |
| $\nabla \cdot \sigma$                     | ElapidViscousStress2D |
| $\nabla \cdot v_s$                        | ElapidVelocityDiv2D |
| $\frac{1}{K_d}\frac{dP}{dt}$              | ElapidSolidElasticTotalPressure |
| $-\frac{\alpha}{K_d}\frac{dP_f}{dt}$      | ElapidSolidElasticFluidPressure |
| $\frac{P-P_f}{(1-\phi_f)\eta_c}$          | ElapidSolidViscous |
| $\frac{\alpha}{K_d B}\frac{dP_f}{dt}$     | ElapidHydroElasticFluidPressure |
| $-\frac{\alpha}{K_d}\frac{dP}{dt}$        | ElapidHydroElasticTotalPressure |
| $-\nabla \cdot \left(\frac{k(\phi_f)}{\mu}(\nabla P_f + \rho_f g)\right)$ | ElapidDarcy(**) |
| $- \frac{P-P_f}{(1- \phi_f)\eta_c}$       | ElapidHydroViscous |
| $\frac{1-\phi_f}{K_\phi}\frac{d P}{d t}$  | ElapidPoroElasticTotalPressure |
| $-\frac{1-\phi_f}{K_\phi}\frac{dP_f}{dt}$ | ElapidPoroElasticFluidPressure |
| $-(1 - \phi_f)\frac{P_f - P}{\eta_c}$     | ElapidPoroViscous |

(**) Gravity inputs for both ElapidDarcy and the default MOOSE Gravity kernel should be positive. Also, this is formulated for domains where depth is increasingly negative. (Ex. 1 km depth corresponds with a domain xmin = -1000 m)

## Material Models

Elapid currently provides a few material models for general usage. We provide both a monophase solid material and a biphasic solid material. Of these, we also provide materials to model linear solid viscosity and nonlinear power-law solid viscosity. The nonlinear materials can be made linear by assigning a stress exponent of 1, however, the linear models bypass some unnecessary calculations. NOTE: All material models are still nonlinear with respect to the effects of porosity.
These models calculate changing local permeability and changing poroelastic moduli.
Permeability is calcuated using a Kozeny-Carman relation, and the poroelastic moduli are calculated with a Mori-Tanaka mixing rule using a provided effective pore aspect ratio.

Some examples of each material are provided in the examples folder.

### Model descriptions
| Description  | Material Name |
| ------------- | ------------- |
| Single solid phase nonlinear (w.r.t viscosity) porous material | SinglePhaseNonLinearViscoElastic |
| Two solid phase nonlinear (w.r.t viscosity) porous material | BiphasicNonLinearViscoElastic |
| Single solid phase linear (w.r.t viscosity) porous material | SinglePhaseLinearViscoElastic |
| Two solid phase linear (w.r.t viscosity) porous material | BiphasicLinearViscoElastic |



### Variables and Parameters
| Symbol  | Descriptor |
| ------------- | ------------- |
| $\sigma$      | Cauchy stress tensor  |
| $\rho_T$      | Total density  |
| $\rho_f$      | Fluid density |
| $g$           | Gravity  |
| $v_s$         | Solid velocity  |
| $K_d$         | Drained bulk modulus  |
| $B$           | Skempton's coefficient |
| $K_\phi$      | Pore bulk modulus  |
| $\alpha$      | Biot's coefficient  |
| $P$           | Total Pressure  |
| $P_f$         | Fluid Pressure  |
| $\phi_f$      | Fluid volume fraction (Porosity)  |
| $\phi_i$      | Solid i volume fraction  |
| $\eta_c$      | Volumetric (Compaction) viscosity |
| $k(\phi_f)$   | Permeability |
| $\mu$         | Fluid viscosity |


## Possible Future Extensions

### Functionality
N-solid phase compatibility
GPU acceleration
ALE for clast deformation in viscous matrix materials

### Physics
Chemical Reactions
Plasticity


## Moose Information
Fork "Elapid" to create a new MOOSE-based application.

For more information see: [https://mooseframework.inl.gov/getting_started/new_users.html#create-an-app](https://mooseframework.inl.gov/getting_started/new_users.html#create-an-app)

## References

McKenzie, D. A. N. (1984). The generation and compaction of partially molten rock. _Journal of petrology_, _25_(3), 713-765.

Yarushina, V. M., & Podladchikov, Y. Y. (2015). (De) compaction of porous viscoelastoplastic media: Model formulation. _Journal of Geophysical Research: Solid Earth_, _120_(6), 4146-4170.
