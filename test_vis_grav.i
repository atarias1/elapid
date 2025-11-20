xmax = 0
ymax = 0
xmin = -1300
ymin = -1300

xres = 10
yres = 10

v_x_BC_top = 1e-5

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
    order = SECOND
    family = LAGRANGE
    initial_condition = 0
  []

  [v_y]
    order = SECOND
    family = LAGRANGE
    initial_condition = 0
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
    expression = '1.5e9 + 9.81*2600*(${ymax} - y)'
  []

  [P_sum]
    type = ParsedFunction
    value = '5e7*exp(-((x+650)^2 + (y+650)^2)/(2*200^2)) + 1.5e9 + 9.81*2600*(${ymax} - y)'
  []

  [dt_cfl]
    type = ParsedFunction
    value = '100*(abs(${xmin}) / ${xres}) / abs(${v_x_BC_top})'
  []
[]

[BCs]
  # Lid-driven: TOP boundary has prescribed velocity
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

  # [left_v_x]
  #   type = ADDirichletBC
  #   boundary = 'left'
  #   value = 0
  #   variable = v_x
  # []

  # [left_v_y]
  #   type = ADDirichletBC
  #   boundary = 'left'
  #   value = 0
  #   variable = v_y
  # []

  # [right_v_x]
  #   type = ADDirichletBC
  #   boundary = 'right'
  #   value = 0
  #   variable = v_x
  # []

  # [right_v_y]
  #   type = ADDirichletBC
  #   boundary = 'right'
  #   value = 0
  #   variable = v_y
  # []

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
      translation = '1300 0 0'
    []
    [x_V_y]
      variable = v_y
      primary = 'left'
      secondary = 'right'
      translation = '1300 0 0'
    []

    [x_p_tot]
      variable = p_tot
      primary = 'left'
      secondary = 'right'
      translation = '1300 0 0'
    []

    [x_p_f]
      variable = p_f
      primary = 'left'
      secondary = 'right'
      translation = '1300 0 0'
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

  [circle_phi_ol]
    type = SmoothCircleIC
    variable = phi_ol
    int_width = 20
    x1 = -650
    y1 = -650
    radius = 100
    outvalue = 0
    invalue = 0.68 # porosity = 0.32
  []

  [circle_phi_f]
    type = SmoothCircleIC
    variable = phi_f
    int_width = 20
    x1 = -650
    y1 = -650
    radius = 100
    outvalue = 0.01
    invalue = 0.32 # porosity = 0.32
  []
[]

[Materials] # For test reduce all pressure by 4 oom
  [vel] # for advection kernels
    type = ADVectorFromComponentVariablesMaterial
    vector_prop_name = 'velocity'
    u = v_x
    v = v_y
  []

  [eta_s]
    type = ADGenericConstantMaterial
    prop_names = 'eta_s'
    prop_values = '1e12'
  []

  [K_d]
    type = ADGenericConstantMaterial
    prop_names = 'K_d'
    prop_values = '50e9'
  []

  [k]
    type = ADGenericConstantMaterial
    prop_names = 'k'
    prop_values = '1e-15'
  []

  [alpha]
    type = ADGenericConstantMaterial
    prop_names = 'alpha'
    prop_values = '0.5'
  []

  [K_phi]
    type = ADGenericConstantMaterial
    prop_names = 'K_phi'
    prop_values = '50e9'
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
    value = 9.81 # test with proper kernel and positive gravity to balance pressure with depth
    variable = v_y
    density = '2600'
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

  [darcy]
    type = DarcyCoupled

    gravity = '0 -9.81 0' #-9.81 original
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

  # [porovisco]
  #   type = PoroViscous
  #   P_f = p_f
  #   P_tot = p_tot
  #   variable = phi_f
  # []
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
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type -pc_factor_shift' #-pc_factor_mat_solver_package  ' #-snes_type'
    petsc_options_value = 'lu    superlu_dist vinewtonssls NONZERO' #superlu_dist     ' #vinewtonssls'

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

  #petsc_options = '-pc_svd_monitor'
  #petsc_options = '-snes_test_jacobian -snes_test_jacobian_view'
  #steady_state_detection = true

  # Relax tolerances
  nl_abs_tol = 1e-7 # Accept 2e-6 as converged
  nl_rel_tol = 1e-7
  nl_max_its = 15

  end_time = 1e14 # around 3000 years
  [TimeSteppers] # CFL h / |v|
    [RampUpDT]
      type = IterationAdaptiveDT
      optimal_iterations = 5
      dt = 5
      growth_factor = 1.2
      cutback_factor = 0.9
    []

    # [MaxDT]
    #   type = FunctionDT
    #   function = dt_cfl
    # []
  []
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
[]

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
