[Mesh]
  file = geom.e
[]

[Variables]
  [./temp]
    initial_condition = '900.0'
  [../]
[]

[Kernels]
  [./heat]
    type = HeatConduction
	  variable = temp
  [../]
  [./fuelHeat]
    type = HeatSource
	  block = '2'
	  function = fuelVolumetricHeat
	  variable = temp
  [../]
  [./monolithHeat]
    type = HeatSource
	  block = '1'
	  function = monolithVolumetricHeat
	  variable = temp
  [../]
[]

[BCs]
  [./outsideBC]
    type = NeumannBC
	  variable = temp
	  value = 0.0
	  boundary = outside
  [../]
  [./heatpipeBC]
    type = FunctionDirichletBC
	  variable = temp
	  boundary = 'HP-0 HP-1 HP-2 HP-3 HP-4 HP-5 HP-6 HP-7 HP-8 HP-9'
	  function = heatpipe_function
    preset = false
  [../]
[]

[Materials]
  [./fuelMat]
    type = GenericFunctionMaterial
	  prop_names = thermal_conductivity
	  prop_values = fuel_thermal_conductivity
	  block = '2'
  [../]
  [./monolithMat]
    type = GenericFunctionMaterial
	  prop_names = thermal_conductivity
	  prop_values = monolithThermalConductivity
	  block = '1 3'
  [../]
  [./fuelDensity]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '12.870'
    block = '2'
  [../]
  [./monolithDensity]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '7.9'
    block = '1 3'
  [../]
[]

[Functions]
  [./monolithThermalConductivity]
    type = ParsedFunction
    value = '0.162'
  [../]
  [./fuel_thermal_conductivity]
    type = ParsedFunction
    value = 0.0141*((fuelTemp+273.15)^0.39)
    vars = 'fuelTemp'
    vals = 'fuelAverageTemp'
  [../]
  [./fuelVolumetricHeat]
    type = ParsedFunction
    value = '60.0'
  [../]
  [./monolithVolumetricHeat]
    type = ParsedFunction
    value = '0.0'
  [../]
  [./heatpipe_function]
    type = VectorPostprocessorFunction
    component = z
	  argument_column = 'distance'
    value_column = 'flux_aggregate'
    vectorpostprocessor_name = adjusted_LSF
  [../]
  [./heatpipe_function_control]
    type = ParsedFunction
    value = '900.0'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = JFNK
  l_max_its = 1000
  nl_max_its = 100
  line_search = none
[]

[Postprocessors]
  [./fuelAverageTemp]
    type = AverageNodalVariableValue
    variable = temp
    block = '2'
    execute_on = 'INITIAL FINAL'
    allow_duplicate_execution_on_initial = true
  [../]
  [./heatFlux_0]
    type = SideFluxIntegral
	  variable = temp
	  boundary = 'HP-0'
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_1]
    type = SideFluxIntegral
	  variable = temp
	  boundary = 'HP-1'
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_2]
    type = SideFluxIntegral
	  variable = temp
	  boundary = 'HP-2'
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_3]
    type = SideFluxIntegral
	  variable = temp
	  boundary = 'HP-3'
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_4]
    type = SideFluxIntegral
	  variable = temp
	  boundary = 'HP-4'
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_5]
    type = SideFluxIntegral
	  variable = temp
	  boundary = 'HP-5'
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_6]
    type = SideFluxIntegral
	  variable = temp
	  boundary = 'HP-6'
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_7]
    type = SideFluxIntegral
	  variable = temp
	  boundary = 'HP-7'
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_8]
    type = SideFluxIntegral
	  variable = temp
	  boundary = 'HP-8'
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_9]
    type = SideFluxIntegral
	  variable = temp
	  boundary = 'HP-9'
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatpipe_min]
    type = NodalExtremeValue
    variable = temp
    block = '3'
    value_type = min
    execute_on = FINAL
  [../]
  [./heatpipe_max]
    type = NodalExtremeValue
    variable = temp
    block = '3'
    value_type = max
    execute_on = FINAL
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
    x_vals = '-25.0 -20.0 -15.0 -10.0 -5.0 0.0 5.0 10.0 15.0 20.0 25.0'
    num_samples = 11
    execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
[]

[Outputs]
  exodus = true
[]