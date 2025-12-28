xmax = 0
ymax = 0
xmin = -30
ymin = -60

xres = 180
yres = 360

v_x_BC_top = 0.0001
# Since the solid viscosity is constant we can make this low

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

  [phi_f]
    order = FIRST
    family = LAGRANGE
  []
[]

[Functions]

  [lithostat]
    type = ParsedFunction # modified for the conditions of Omlin et al., 2018
    expression = '1*(2*(1-0.002)+ 1*(0.002))*(${ymax} - y)' # Formulated to increase with depth
  []
[]

[BCs]

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
    value = 0
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
    value = 0
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
      translation = '30 0 0'
    []
    [x_V_y]
      variable = v_y
      primary = 'left'
      secondary = 'right'
      translation = '30 0 0'
    []

    [x_p_tot]
      variable = p_tot
      primary = 'left'
      secondary = 'right'
      translation = '30 0 0'
    []

    [x_p_f]
      variable = p_f
      primary = 'left'
      secondary = 'right'
      translation = '30 0 0'
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
    int_width = 0.5
    x1 = -15
    y1 = -48
    radius = 3
    outvalue = 0.002
    invalue = 0.004 #
  []
[]

[Materials]
  [vel] # for advection kernels
    type = ADVectorFromComponentVariablesMaterial
    vector_prop_name = 'velocity'
    u = v_x
    v = v_y
  []

  [single_linear]
    type = SinglePhaseLinearViscoElastic
    G = 166.66
    K = 333.33
    a_eta = 0.4
    aspect_ratio = 0.25
    eta_s_0 = 1
    fluid_K = 10.0
    k_ref = 1.0
    mu = 1
    nk = 3
    phi_0 = 0.002
    phi_f = phi_f
    phi_ref = 0.025
    rho_f = 1
    rho_s = 2
    zeta = 2.0
  []
[]

[Kernels]

  [div_stress_x]
    type = ElapidViscousStress2D
    P_tot = p_tot
    V_x_s = v_x
    V_y_s = v_y
    component = 0
    variable = v_x
  []

  [div_stress_y]
    type = ElapidViscousStress2D
    P_tot = p_tot
    V_x_s = v_x
    V_y_s = v_y
    component = 1
    variable = v_y
  []

  [gravity]
    type = ADGravity
    value = 1 # needs to be positive
    variable = v_y
    density = rho_T
  []

  [div_v]
    type = ElapidVelocityDiv2D
    V_x_s = v_x
    V_y_s = v_y
    variable = p_tot
  []

  [SETP]
    type = ElapidSolidElasticTotalPressure
    variable = p_tot
  []

  [SEFP]
    type = ElapidSolidElasticFluidPressure
    P_f = p_f
    variable = p_tot
  []

  [SV]
    type = ElapidSolidViscous
    P_f = p_f
    phi_f = phi_f
    variable = p_tot
  []

  [darcy]
    type = ElapidDarcy

    gravity = '0 1.0 0' # needs to be positive
    mu = 1.0
    rho_f = 1
    variable = p_f
  []

  [HEFP]
    type = ElapidHydroElasticFluidPressure
    variable = p_f
  []

  [HETP]
    type = ElapidHydroElasticTotalPressure
    P_tot = p_tot
    variable = p_f
  []

  [HV]
    type = ElapidHydroViscous
    P_tot = p_tot
    phi_f = phi_f
    variable = p_f
  []

  [time_phi_f]
    type = MassLumpedTimeDerivative
    variable = phi_f
  []

  [PETP]
    type = ElapidPoroElasticTotalPressure
    P_tot = p_tot
    variable = phi_f
  []

  [PEFP]
    type = ElapidPoroElasticFluidPressure
    P_f = p_f
    variable = phi_f
  []

  [PV]
    type = ElapidPoroViscous
    P_f = p_f
    P_tot = p_tot
    variable = phi_f
  []
[]

[AuxVariables]
  [bounds_dummy]
    order = CONSTANT
    family = MONOMIAL
  []
  [bounds_dummy_f]
    order = FIRST
    family = LAGRANGE
  []
[]

[Bounds]
  [phi_f_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy_f
    bounded_variable = phi_f
    bound_type = upper
    bound_value = 0.9999
  []
  [phi_f_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy_f
    bounded_variable = phi_f
    bound_type = lower
    bound_value = 0
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type -pc_factor_shift'
    petsc_options_value = 'lu    superlu_dist vinewtonssls NONZERO'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = bdf2
  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = true

  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-7
  nl_max_its = 25

  l_max_its = 1000

  num_steps = 100 # Need 100 to get a checkpoint
  #end_time = 1e14 # around 3000 years
  [TimeSteppers] # CFL h / |v|
    [RampUpDT]
      type = IterationAdaptiveDT
      optimal_iterations = 7
      dt = 0.001
      growth_factor = 1.5
      cutback_factor = 0.8
    []

    [MaxDT]
      type = FunctionDT
      function = 0.5
    []
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
