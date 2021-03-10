[Mesh]
  file = heatflux_box.e
[]

[Variables]
  [./temp]
    initial_condition = '300.0'
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
    type = DirichletBC
	  variable = temp
	  boundary = BC_left
    value = 300.0
  [../]
  [./BC_right]
    type = DirichletBC
	  variable = temp
	  boundary = BC_right
    value = 400.0
  [../]
[]

[Materials]
  [./blockMat]
    type = GenericConstantMaterial
	  prop_names = thermal_conductivity
	  prop_values = 0.05
	  block = 1
  [../]
[]

[Functions]
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  l_max_its = 1000
  nl_max_its = 100
[]

[Postprocessors]
  [./heatFlux_left]
    type = SideFluxIntegral
	  variable = temp
	  boundary = BC_left
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_right]
    type = SideFluxIntegral
	  variable = temp
	  boundary = BC_right
	  diffusivity = thermal_conductivity
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
[]

[Preconditioning]
[]

[VectorPostprocessors]
[]

[Outputs]
  exodus = true
[]