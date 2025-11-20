xmax = 0
ymax = 0
xmin = -50 # 5mm
ymin = -80

xres = 20
yres = 40



v_x_BC_top = 1 # Test with higher strain rate 1e-9 is around do_ep = 1e-12 v = edot/L

## Characteristic scales ##
# p_compact = 147035 # Pa
# v_compact = 1.25e-13 # m/s will make res really big not sure if usable for conditioning


# Use this input to equilibriate 

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = ${xres}
  ny = ${yres}

  xmax = ${xmax}
  ymax = ${ymax}
  xmin = ${xmin}
  ymin = ${ymin}

  elem_type = QUAD8
[]

[Variables]

  [p_tot]
    order = FIRST
    family = LAGRANGE
  []

  [p_f]
    order = FIRST
    family = LAGRANGE
  []

  [v_x]
    #scaling = 1e-7 # do the residual analysis later and ensure the scaling is correct
    order = SECOND
    family = LAGRANGE
    initial_condition = 0
  []

  [v_y]
    #scaling = 1e-7
    order = SECOND
    family = LAGRANGE
    initial_condition = 0
  []

  [phi_f]
    order = FIRST
    family = LAGRANGE
  []
[]

[Functions]

  [lithostat]
    type = ParsedFunction
    expression = '1 + 1*2.6*(${ymax} - y)' # Formulated to increase with depth
  []

  # [P_sum]
  #   type = ParsedFunction
  #   value = '5e7*exp(-((x+650)^2 + (y+650)^2)/(2*200^2)) + 1.5e9 + -9.81*2600*(${ymax} - y)'
  # []

  # [dt_cfl]
  #   type = ParsedFunction
  #   value = '100*(abs(${xmin}) / ${xres}) / abs(${v_x_BC_top})'
  # []
[]

[BCs]
  # Shear in +x -y
  [top_v_x]
    type = ADFunctionDirichletBC
    boundary = 'top'
    function = '${v_x_BC_top}'
    variable = v_x
  []

  [top_v_y]
    type = ADDirichletBC
    boundary = 'top'
    value = 0
    variable = v_y
  []

  # No-slip walls
  [bottom_v_y]
    type = ADDirichletBC
    boundary = 'bottom'
    value = 0
    variable = v_y
  []

  [bottom_v_x]
    type = ADFunctionDirichletBC
    boundary = 'bottom'
    function = '-1*${v_x_BC_top}'
    variable = v_x
  []

  [top_p_tot]
    type = ADDirichletBC
    boundary = 'top'
    value = 1
    variable = p_tot
  []

  # [top_p_tot]
  #   type = ADNeumannBC
  #   boundary = 'top'
  #   value = 0
  #   variable = p_tot
  # []

  [bottom_p_tot]
    type = ADFunctionDirichletBC
    boundary = 'bottom'
    function = lithostat
    variable = p_tot
  []

  # [bottom_p_tot]
  #   type = ADNeumannBC
  #   boundary = 'bottom'
  #   value = 0
  #   variable = p_tot
  # []

  [top_p_f]
    type = ADDirichletBC
    boundary = 'top'
    value = 1
    variable = p_f
  []

  # [top_p_f]
  #   type = ADNeumannBC
  #   boundary = 'top'
  #   value = 0
  #   variable = p_f
  # []


  [bottom_p_f]
    type = ADFunctionDirichletBC
    boundary = 'bottom'
    function = lithostat
    variable = p_f
  []

  # [bottom_p_f]
  #   type = ADNeumannBC
  #   boundary = 'bottom'
  #   value = 0
  #   variable = p_f
  # []

  [Periodic]
    [x_V_x]
      variable = v_x
      primary = 'left'
      secondary = 'right'
      translation = '50 0 0'
    []
    [x_V_y]
      variable = v_y
      primary = 'left'
      secondary = 'right'
      translation = '50 0 0'
    []

    [x_p_tot]
      variable = p_tot
      primary = 'left'
      secondary = 'right'
      translation = '50 0 0'
    []

    [x_p_f]
      variable = p_f
      primary = 'left'
      secondary = 'right'
      translation = '50 0 0'
    []

    [x_phi_f]
      variable = phi_f
      primary = 'left'
      secondary = 'right'
      translation = '50 0 0'
    []
  []
[]

[ICs]
  [P_tot_init]
    type = FunctionIC
    function = lithostat
    variable = p_tot
  []

  [P_f_init]
    type = FunctionIC
    function = lithostat
    variable = p_f
  []

  [circle_phi_f]
    type = SmoothCircleIC
    variable = phi_f
    int_width = 0.0
    x1 = -25
    y1 = -60
    radius = 3
    outvalue = 0.001
    invalue = 0.1 # porosity = 0.32
  []
[]

[Materials] # For test reduce all pressure by 4 oom
  [vel] # for advection kernels
    type = ADVectorFromComponentVariablesMaterial
    vector_prop_name = 'velocity'
    u = v_x
    v = v_y
  []

  # [serp]
  #   type = HomogenousViscoElasticNonLin
  #   G = 166
  #   K = 333
  #   a_eta = 0.4
  #   aspect = 0.01
  #   eta_0 = 1.1
  #   k_ref = 1
  #   max_eta_s = 100
  #   mu = 1
  #   nk = 3
  #   nsigma = 3.8
  #   phi_0 = 0.01
  #   phi_f = phi_f
  #   phi_ref = 1
  #   rho_f = 1
  #   rho_s = 2.6
  #   v_x = v_x
  #   v_y = v_y
  #   water_K = 2.2
  #   zeta = 2.0

  #   output_properties = 'k eta_s K_d alpha K_phi phi_s rho_T eta_s K_d alpha Skempton comp_eta'
  #   outputs = exodus
   
  # []

  [serp]
    type = MetaSerpentinite # lowered visc by 2 oom
    a_eta = 0.4         # factor
    atg_G = 39        # Pa
    atg_K = 75        # Pa
    atg_aspect = 0.01   # ratio
    atg_eta_0 = 1.1    # characteristic strain rate 1/s original e10
    k_ref = 1e-2     # m^2 - original estimate 7.9e-21  
    max_eta_s = 1000  # Pas cutoff max viscosity - use olivines - original e22
    mu = 1           # Pas water viscosity
    nk = 3              # K-C exp.
    ol_G = 81         # Pa
    ol_K = 138       # Pa
    ol_aspect = 0.4     # ratio
    ol_eta_0 = 100   # Pas - original e22
    phi_0 = 0.01        # vol frac.
    phi_f = phi_f       # COUPLED VAR
    phi_ol = phi_ol     # COUPLED VAR
    phi_ref = 0.005     # vol frac.
    rho_atg = 2.6      # kg/m^3
    rho_f = 1.0        # kg/m^3
    rho_ol = 3.3       # kg/m^3
    v_x = v_x           # COUPLED VAR
    v_y = v_y           # COUPLED VAR
    water_K = 2.2     # Pa
    zeta = 2.0          # Pas modifies shear viscosity to bulk viscosity
    
    output_properties = 'k eta_s K_d alpha K_phi phi_atg rho_T eta_atg_eff comp_eta ol_K_d atg_K_d ol_alpha atg_alpha Skempton comp_eta'
    outputs = exodus
  []
[]


[Kernels]

  [div_stress_x]
    type = ViscousStress2D
    P_tot = p_tot
    V_x_s = v_x
    V_y_s = v_y
    component = 0
    variable = v_x
  []

  [div_stress_y]
    type = ViscousStress2D
    P_tot = p_tot
    V_x_s = v_x
    V_y_s = v_y
    component = 1
    variable = v_y
  []

  [gravity]
    type = ADGravity
    value = 1 # test with proper kernel and positive gravity to balance pressure with depth (positive seems to be correct according to elastic model?)
    variable = v_y
    density = rho_T
  []

  [div_v]
    type = VelocityDiv2D
    V_x_s = v_x
    V_y_s = v_y
    variable = p_tot
  []

  [elastic1]
    type = MatCoefTimeDerivative
    variable = p_tot
  []

  [hydroelastic]
    type = MatCoefTimeDerivativeCoupled
    P_f = p_f
    variable = p_tot
  []

  [ptot_viscous]
    type = ViscousRheology
    P_f = p_f
    phi_f = phi_f
    variable = p_tot 
  []

  [darcy]
    type = DarcyCoupled

    gravity = '0 1 0' #-9.81 original # needs to be positive?
    mu = 1
    rho_f = 1
    variable = p_f
  []

  [time_hydro]
    type = HydroMatCoefTime
    variable = p_f
  []

  [time_hydro_coupled]
    type = HydroMatCoefTimeCoupled
    P_tot = p_tot
    variable = p_f
  []

  [pf_viscous]
    type = ViscousRheologyNegative
    P_tot = p_tot
    phi_f = phi_f
    variable = p_f
  []

  ## conservation kernels

  [time_phi_f]
    type = MassLumpedTimeDerivative
    variable = phi_f
  []

  [poromechano]
    type = PoroMechanoMatCoefTime
    P_tot = p_tot
    variable = phi_f
  []

  [porohydro]
    type = PoroHydroMatCoefTime
    P_f = p_f
    variable = phi_f
  []

  [porovisco]
    type = PoroViscous
    P_f = p_f
    P_tot = p_tot
    variable = phi_f
  []
[]

# [AuxVariables]
#   [bounds_dummy]
#     order = FIRST
#     family = LAGRANGE
#   []
# []

# [Bounds]

#   [phi_f_upper_bound]
#     type = ConstantBounds
#     variable = bounds_dummy
#     bounded_variable = phi_f
#     bound_type = upper
#     bound_value = 0.9999 # It cannot be one or the denominator goes to 0
#   []
#   [phi_f_lower_bound]
#     type = ConstantBounds
#     variable = bounds_dummy
#     bounded_variable = phi_f
#     bound_type = lower
#     bound_value = 0
#   []
# []

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type -pc_factor_shift' #-pc_factor_mat_solver_package  ' #-snes_type'
    petsc_options_value = 'lu    superlu_dist vinewtonssls NONZERO' #superlu_dist     ' #vinewtonssls'

    # Gives out of memory error --> occurs with mumps apparently?
    # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type -pc_factor_shift' #-pc_factor_mat_solver_package  ' #-snes_type'
    # petsc_options_value = 'lu    mumps vinewtonssls NONZERO' #superlu_dist     ' #vinewtonssls'


    # petsc_options_iname = '-pc_type -pc_hypre_type -snes_type'
    # petsc_options_value = 'hypre boomeramg vinewtonssls'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = bdf2
  automatic_scaling = false
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = true
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'svd'
  # petsc_options = '-pc_svd_monitor'
  #petsc_options = '-snes_test_jacobian -snes_test_jacobian_view'
  #steady_state_detection = true

  # Relax tolerances
  nl_abs_tol = 1e-7 # Accept 2e-6 as converged
  nl_rel_tol = 1e-7
  nl_max_its = 25

  l_max_its = 1000

  num_steps = 1000 # Need 100 to get a checkpoint  
  #end_time = 1e14 # around 3000 years
  [TimeSteppers] # CFL h / |v|
    [RampUpDT]
      type = IterationAdaptiveDT
      optimal_iterations = 7
      dt = 0.02
      growth_factor = 1.5
      cutback_factor = 0.8
    []

    [MaxDT]
      type = FunctionDT
      function = 1e6
    []
  []
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
  interval = 2

  [./my_checkpoint]
    type = Checkpoint
    num_files = 4
    interval = 100
  [../]
[]

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
