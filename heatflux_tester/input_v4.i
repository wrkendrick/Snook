[Mesh]
  file = heatflux_box_v2.e
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
  [./BC_midleft]
    type = DirichletBC
	  variable = temp
	  boundary = BC_midleft
    value = 350.0
  [../]
  [./BC_midright]
    type = DirichletBC
	  variable = temp
	  boundary = BC_midright
    value = 350.0
  [../]
  [./BC_right]
    type = DirichletBC
	  variable = temp
	  boundary = BC_right
    value = 400.0
  [../]
[]

[Materials]
  [./blockMat_left]
    type = GenericConstantMaterial
	  prop_names = 'thermal_conductivity therm_cond_left'
	  prop_values = '0.05 0.05'
	  block = block_1
  [../]
  [./blockMat_right]
    type = GenericConstantMaterial
	  prop_names = 'thermal_conductivity therm_cond_right'
	  prop_values = '0.1 0.1'
	  block = block_2
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
	  diffusivity = therm_cond_left
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_midleft]
    type = SideFluxIntegral
	  variable = temp
	  boundary = BC_midleft
	  diffusivity = therm_cond_left
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_midright]
    type = SideFluxIntegral
	  variable = temp
	  boundary = BC_midright
	  diffusivity = therm_cond_right
	  execute_on = 'INITIAL NONLINEAR FINAL'
  [../]
  [./heatFlux_right]
    type = SideFluxIntegral
	  variable = temp
	  boundary = BC_right
	  diffusivity = therm_cond_right
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