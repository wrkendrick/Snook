[Mesh]
  file = heatflux_box_v4.e
[]

[Variables]
  [./temp]
    initial_condition = '1200.0'
  [../]
[]

[Kernels]
  [./heat]
    type = HeatConduction
	  variable = temp
  [../]
[]

[BCs]
  [./BC_alt]
    type = NeumannBC
	  variable = temp
	  value = 0.0
	  boundary = BC_alt
  [../]
  [./BC_left]
    type = FunctionDirichletBC
	  variable = temp
	  boundary = 'flux_0 flux_1 flux_2 flux_3 flux_4 flux_5 flux_6 flux_7 flux_8 flux_9'
	  function = htpipe_function
    preset = false
  [../]
  [./BC_right]
    type = DirichletBC
	  variable = temp
	  boundary = BC_right
    value = 1200.0
  [../]
[]

[Materials]
  [./blockMat_left]
    type = GenericConstantMaterial
	  prop_names = thermal_conductivity
	  prop_values = 13
	  block = block_1
  [../]
  [./blockMat_right]
    type = GenericConstantMaterial
	  prop_names = thermal_conductivity
	  prop_values = 13
	  block = block_2
  [../]
[]

[Functions]
  [./htpipe_function]
    type = VectorPostprocessorFunction
    component = z
	  argument_column = 'distance'
    value_column = 'flux_aggregate'
    vectorpostprocessor_name = adjusted_LSF
  [../]
  [./htpipe_function_control]
    type = ParsedFunction
    value = '900.0'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  l_max_its = 1000
  nl_max_its = 100
  nl_rel_tol = 1e-8
[]

[Postprocessors]
  [./heatFlux_0]
    type = SideFluxIntegral
	  variable = temp
	  boundary = flux_0
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_1]
    type = SideFluxIntegral
	  variable = temp
	  boundary = flux_1
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_2]
    type = SideFluxIntegral
	  variable = temp
	  boundary = flux_2
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_3]
    type = SideFluxIntegral
	  variable = temp
	  boundary = flux_3
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_4]
    type = SideFluxIntegral
	  variable = temp
	  boundary = flux_4
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_5]
    type = SideFluxIntegral
	  variable = temp
	  boundary = flux_5
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_6]
    type = SideFluxIntegral
	  variable = temp
	  boundary = flux_6
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_7]
    type = SideFluxIntegral
	  variable = temp
	  boundary = flux_7
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_8]
    type = SideFluxIntegral
	  variable = temp
	  boundary = flux_8
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_9]
    type = SideFluxIntegral
	  variable = temp
	  boundary = flux_9
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
[]

[Preconditioning]
[]

[VectorPostprocessors]
  [./flux_aggregate]
    type = VectorOfPostprocessors
    postprocessors = 'heatFlux_0 heatFlux_1 heatFlux_2 heatFlux_3 heatFlux_4 heatFlux_5 heatFlux_6 heatFlux_7 heatFlux_8 heatFlux_9'
    execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./adjusted_LSF]
    type = LeastSquaresFit
    vectorpostprocessor = flux_aggregate
    x_name = 'distance'
    y_name = 'flux_aggregate'
    x_vals = '-0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5'
    num_samples = 11
    execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
[]

[Outputs]
  exodus = true
[]