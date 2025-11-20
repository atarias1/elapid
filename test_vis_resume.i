
ymax = 0
ymin = -0.005
v_x_BC_top = 1e-10*${ymin} 

[Mesh]
  #Serial number should match corresponding Executioner parameter
  file = test_vis_restruct_my_checkpoint_cp/0100-mesh.cpa.gz
  #This method of restart is only supported on serial meshes
  #distribution = serial
[]

[Problem]
  #Note that the suffix is left off in the parameter below.
  restart_file_base = test_vis_restruct_my_checkpoint_cp/LATEST  # You may also use a specific number here
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
  []

  [v_y]
    #scaling = 1e-7
    order = SECOND
    family = LAGRANGE
  []

  [phi_ol]
    order = CONSTANT
    family = MONOMIAL
  []

  [phi_f]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]

  [lithostat]
    type = ParsedFunction
    expression = '1.5e9 + 9.81*2600*(${ymax} - y)' # Formulated to increase with depth
  []
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
    value = 1.5e9
    variable = p_tot
  []

  [bottom_p_tot]
    type = ADFunctionDirichletBC
    boundary = 'bottom'
    function = lithostat
    variable = p_tot
  []

  [top_p_f]
    type = ADDirichletBC
    boundary = 'top'
    value = 1.5e9
    variable = p_f
  []

  [bottom_p_f]
    type = ADFunctionDirichletBC
    boundary = 'bottom'
    function = lithostat
    variable = p_f
  []

  [Periodic]
    [x_V_x]
      variable = v_x
      primary = 'left'
      secondary = 'right'
      translation = '0.005 0 0'
    []
    [x_V_y]
      variable = v_y
      primary = 'left'
      secondary = 'right'
      translation = '0.005 0 0'
    []

    [x_p_tot]
      variable = p_tot
      primary = 'left'
      secondary = 'right'
      translation = '0.005 0 0'
    []

    [x_p_f]
      variable = p_f
      primary = 'left'
      secondary = 'right'
      translation = '0.005 0 0'
    []
  []
[]

[Materials]
  [vel] # for advection kernels
    type = ADVectorFromComponentVariablesMaterial
    vector_prop_name = 'velocity'
    u = v_x
    v = v_y
  []

  [serp]
    type = MetaSerpentinite # lowered visc by 2 oom
    a_eta = 0.4         # factor
    atg_G = 39e9        # Pa
    atg_K = 75e9        # Pa
    atg_aspect = 0.01   # ratio
    atg_eta_0 = 1.1e10    # characteristic strain rate 1/s original e10
    k_ref = 1e-19     # m^2 - original estimate 7.9e-21  
    max_eta_s = 1.1e22  # Pas cutoff max viscosity - use olivines - original e22
    mu = 1e-3           # Pas water viscosity
    nk = 3              # K-C exp.
    ol_G = 81e9         # Pa
    ol_K = 138e9        # Pa
    ol_aspect = 0.4     # ratio
    ol_eta_0 = 1.1e22   # Pas - original e22
    phi_0 = 0.01        # vol frac.
    phi_f = phi_f       # COUPLED VAR
    phi_ol = phi_ol     # COUPLED VAR
    phi_ref = 0.005     # vol frac.
    rho_atg = 2600      # kg/m^3
    rho_f = 1000        # kg/m^3
    rho_ol = 3300       # kg/m^3
    v_x = v_x           # COUPLED VAR
    v_y = v_y           # COUPLED VAR
    water_K = 2.2e9     # Pa
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
    value = 9.81 # test with proper kernel and positive gravity to balance pressure with depth (positive seems to be correct according to elastic model?)
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

    gravity = '0 9.81 0' #-9.81 original - needs to be positive?
    mu = 1e-3
    rho_f = 1000
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
  [time_phi_ol]
    type = MassLumpedTimeDerivative
    variable = phi_ol
  []

  [advect_phi_ol]
    type = ADConservativeAdvection
    variable = phi_ol
    velocity = 'velocity'
    # upwinding_type = full
  []

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

[DGKernels] # Need to use constant monomial vars I think for advected vars
  [phi_ol_advect]
    type = ADDGAdvection
    variable = phi_ol
    velocity = 'velocity'
  []
[]

[AuxVariables]
  [bounds_dummy]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Bounds]
  [phi_ol_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = phi_ol
    bound_type = upper
    bound_value = 1
  []
  [phi_ol_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = phi_ol
    bound_type = lower
    bound_value = 0
  []

  [phi_f_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = phi_f
    bound_type = upper
    bound_value = 0.9999
  []
  [phi_f_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = phi_f
    bound_type = lower
    bound_value = 0
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type -pc_factor_shift' #-pc_factor_mat_solver_package  ' #-snes_type'
    petsc_options_value = 'lu    superlu_dist vinewtonssls NONZERO' #superlu_dist     ' #vinewtonssls'

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
  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = true
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'svd'
  # petsc_options = '-pc_svd_monitor'
  #petsc_options = '-snes_test_jacobian -snes_test_jacobian_view'
  #steady_state_detection = true

  # Relax tolerances
  nl_abs_tol = 1e-8 # Accept 2e-6 as converged
  nl_rel_tol = 1e-8
  nl_max_its = 15

  l_max_its = 1000

  end_time = 1e14 # around 3000 years
  [TimeSteppers] # CFL h / |v|
    [RampUpDT]
      type = IterationAdaptiveDT
      optimal_iterations = 4
      dt = 100
      growth_factor = 1.25
      cutback_factor = 0.8
    []

    [MaxDT]
      type = FunctionDT
      function = '1e10'
    []
  []
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
  interval = 2

  # [./my_checkpoint]
  #   type = Checkpoint
  #   num_files = 4
  #   interval = 100
  # [../]
[]

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
