
Wridge = 3
Wsub = 10
Xo = 0
current_active_region = 6
damp_method = 4
fem_nev = 4
fem_order = 1
geom_x_one_dimension = 0
glob_grid_dxright = 1
glob_grid_showstruc = 1
glob_temp = 300
kpsc_skip_vpiezo = 1
las_length = 450
las_ref1 = 0.18
las_ref2 = 0.95
las_wavegdloss = 12
las_ygridstr = 3.23
linear_solver = 0
mode_solve_lambda = 0.42
n_energy_grid = 1000
newt_rich_method = 0
per_bias_schroedinger = 0
plot_custom_biasX = Current[0]
plot_custom_biasY = ModalGain[0,0]
plot_custom_energy = ModalGain[0,0]
plot_custom_spatial = Mode[0,0]
plot_gain_axis = lambda
plot_num_xpts = 30
plot_num_ypts = 100
plot_region = 14
plot_spatial = 1
plot_xmax = 5
plot_ymax = 5.288
plot_ymin = 0.1
plot_yposition = 0.1
qw_drift = 1
run_auger_rec = 0
run_bpm_numx = 30
run_bpm_numy = 1000
run_bulk_bgr = 1
run_carr_loss = 0
run_gain_maxdens = 2e+013
run_gain_mindens = 4e+012
run_gain_mintemp = 293
run_gain_numpts = 15
run_profile_type = custom
run_residual = 1e-005
run_threshold_bias = 3
sim_llatticetemp = 1


xyrect 1
	element_type = 1
	color = 12
	x = Xo
	y = 0.1
	width = Wsub
	height = 3
	dy_top = 0.001
	dy_bottom = 0.3
	dy_max = 1.5
	use_default_y = 0
	physicaldata
		dope_dendop = 3e+018
		elem_matfile = InGaN_wz.mat
		kpmat_DIELOPT = 2.54^2
		kpmat_MOBN0 = 0.041
		reg_filmat = auto
		reg_material = auto
		reg_regtype = generic
		regmat_thermconduct = 1.3
	end physicaldata
end xyrect

xyrect 2
	element_type = 1
	color = 12
	x = Xo
	y = 3 rel element 1
	width = Wsub
	height = 0.1
	dy_top = 0.001
	dy_bottom = 0.001
	dy_max = 0.05
	use_default_y = 0
	physicaldata
		dope_dendop = 3e+018
		elem_matfile = InGaN_wz.mat
		haensch_tunnel_length = 0.5
		kpmat_DIELOPT = 2.98^2
		kpmat_HLAFF = 5.9443-.2
		kpmat_MOBN0 = 0.039
		reg_conc_i = 0.1
		reg_filmat = auto
		reg_material = auto
		reg_regtype = generic
		regmat_thermconduct = 1
	end physicaldata
end xyrect

xyrect 3
	element_type = 1
	color = 14
	x = Xo
	y = 0.1 rel element 2
	width = Wsub
	height = 1.2
	dy_top = 0.001
	dy_bottom = 0.001
	dy_max = 0.6
	use_default_y = 0
	physicaldata
		dope_dendop = 3e+018
		elem_matfile = AlGaN_wz.mat
		kpmat_DIELOPT = 2.48^2
		kpmat_MOBN0 = 0.001
		reg_conc_i = 0.07
		reg_filmat = auto
		reg_material = auto
		reg_regtype = generic
		regmat_thermconduct = 0.2
	end physicaldata
end xyrect

xyrect 4
	element_type = 1
	color = 12
	x = Xo
	y = 1.2 rel element 3
	width = Wsub
	height = 0.1
	dy_top = 0.0004
	dy_bottom = 0.001
	dy_max = 0.05
	use_default_y = 0
	physicaldata
		dope_dendop = 7e+017
		elem_matfile = InGaN_wz.mat
		kpmat_DIELOPT = 2.54^2
		kpmat_HLAFF = 6.0928-.2
		kpmat_MOBN0 = 0.055
		reg_filmat = auto
		reg_material = auto
		reg_regtype = generic
		regmat_thermconduct = 1.3
	end physicaldata
end xyrect

xyrect 5
	element_type = 2
	color = 10
	color2 = 9
	x = Xo
	y = 0.1 rel element 4
	width = Wsub
	height = 0.038
	dy_top = 0.0004
	dy_bottom = 0.0004
	dy_max = 0.019
	use_default_y = 0
	physicaldata
		dope_dendop = 7e+016
		elem_matfile = InGaN_wz.mat
		elem_numwells = 2
		elem_well_thick = 0.004
		kploop_kanglez_min = 0
		kploop_wavevec_max = 0.08
		kpmat_DIELOPT = 3^2
		kpmat_MOBN0 = 0.03
		kpsc_skip_vcoul = 1
		kpsc_skip_vpiezo = 1
		kpsc_skip_vxc = 0
		kpsys_carr_temp = glob_temp
		kpsys_latt_temp = glob_temp
		kpsys_nconc = run_gain_maxdens
		kpsys_pconc = run_gain_maxdens
		mobel_enhance = 10
		mobhl_enhance = 10
		reg_conc_i = 0.15
		reg_filkp = auto
		reg_filmat = auto
		reg_material = auto
		reg_regtype = qw
		regmat_SRH_taun = 8e-010
		regmat_SRH_taup = 8e-010
		regmat_thermconduct = 0.2
	end physicaldata
	physicaldataB
		dope_dendop = 7e+016
		elem_matfile = InGaN_wz.mat
		elem_sch_thick = 0.01
		kpmat_DIELOPT = 2.61^2
		kpmat_MOBN0 = 0.085
		reg_conc_i = 0.02
		reg_filmat = auto
		reg_material = auto
		reg_regtype = sch
		regmat_thermconduct = 0.2
	end physicaldataB
end xyrect

xyrect 6
	element_type = 1
	color = 14
	x = Xo
	y = 0.038 rel element 5
	width = Wsub
	height = 0.02
	dy_top = 0.002
	dy_bottom = 0.0004
	dy_max = 0.01
	use_default_y = 0
	physicaldata
		dope_dendopa = 1e+019
		elem_matfile = AlGaN_wz.mat
		kpmat_DIELOPT = 2.23^2
		kpmat_HLAFF = 6.202-.1
		kpmat_MOBP0 = 0.001
		reg_conc_i = 0.2
		reg_filmat = auto
		reg_material = auto
		reg_regtype = generic
		regmat_thermconduct = 0.8
	end physicaldata
end xyrect

xyrect 7
	element_type = 1
	color = 12
	x = Xo
	y = 0.02 rel element 6
	width = Wsub
	height = 0.1
	dy_top = 0.001
	dy_bottom = 0.002
	dy_max = 0.05
	use_default_y = 0
	physicaldata
		dope_dendopa = 5e+018
		elem_matfile = InGaN_wz.mat
		kpmat_DIELOPT = 2.54^2
		kpmat_HLAFF = 6.0928-.2
		kpmat_MOBP0 = 0.0015
		reg_filmat = auto
		reg_material = auto
		reg_regtype = generic
		regmat_thermconduct = 1.3
	end physicaldata
end xyrect

xyrect 8
	element_type = 1
	color = 14
	x = Xo
	y = 0.1 rel element 7
	width = Wridge
	height = 0.6
	dx_left = 0.15
	dx_right = 0.15
	dx_max = 0.75
	use_default_x = 0
	dy_top = 0.003
	dy_bottom = 0.001
	dy_max = 0.3
	use_default_y = 0
	physicaldata
		dope_dendopa = 1e+020
		elem_matfile = AlGaN_wz.mat
		kpmat_DIELOPT = 2.48^2
		kpmat_HLAFF = 6.127471
		kpmat_MOBP0 = 5e-005
		reg_conc_i = 0.07
		reg_filmat = auto
		reg_material = auto
		reg_regtype = generic
		regmat_thermconduct = 0.2
	end physicaldata
end xyrect

xyrect 9
	element_type = 1
	color = 12
	x = Xo
	y = 0.6 rel element 8
	width = Wridge
	height = 0.03
	physicaldata
		dope_dendopa = 2e+020
		elem_matfile = InGaN_wz.mat
		kpmat_DIELOPT = 2.54^2
		kpmat_MOBP0 = 0.001
		reg_filmat = auto
		reg_material = auto
		reg_regtype = generic
		regmat_thermconduct = 1.3
	end physicaldata
end xyrect

xyrect 10
	element_type = 3
	color = 7
	x = 3.25 rel element 8
	y = 0.1 rel element 7
	width = (Wsub-Wridge)/2
	height = 0.63
	physicaldata
		elem_matfile = Oxide.mat
		reg_filmat = auto
		reg_material = auto
		reg_regtype = oxide
	end physicaldata
end xyrect

xyrect 11
	element_type = 3
	color = 7
	x = -3.25 rel element 8
	y = 0.1 rel element 7
	width = (Wsub-Wridge)/2
	height = 0.63
	physicaldata
		elem_matfile = Oxide.mat
		reg_filmat = auto
		reg_material = auto
		reg_regtype = oxide
	end physicaldata
end xyrect

xyrect 12
	element_type = 4
	color = 11
	x = Xo
	y = 0
	width = Wsub
	height = 0.1
	physicaldata
		elec_electrode = 1
		elec_temperature = glob_temp
		elec_thermboundary = cauchy
		elec_thermconduct = 0.012
	end physicaldata
end xyrect

xyrect 13
	element_type = 4
	color = 11
	x = 0
	y = 0.03 rel element 9
	width = Wridge
	height = 0.1
	physicaldata
		elec_electrode = 0
		elec_temperature = glob_temp
		elec_thermboundary = none
	end physicaldata
end xyrect


bias 0
	bias_bias = voltage
	bias_biasvalue = 0
	bias_electrode = 0
	bt_active = 500
	bt_biasnumb = 0
	sim_Tamb = glob_temp
	sim_filopt = opt0
	sim_lstead = 1
	sim_resbnd = run_residual
end bias 0


bias 1
	bias_bias = voltage
	bias_biasvalue = 1
	bias_electrode = 0
	bt_active = 502
	bt_biasnumb = 1
	n_steps = 4
	sim_Tamb = glob_temp
	sim_lstead = 1
	sim_resbnd = run_residual
end bias 1


bias 2
	bias_bias = voltage
	bias_biasvalue = 2
	bias_electrode = 0
	bt_active = 502
	bt_biasnumb = 2
	n_steps = 8
	sim_Tamb = glob_temp
	sim_lstead = 1
	sim_resbnd = run_residual
end bias 2


bias 3
	bias_bias = voltage
	bias_biasvalue = 3
	bias_electrode = 0
	bt_active = 502
	bt_biasnumb = 3
	n_steps = 8
	sim_Tamb = glob_temp
	sim_lstead = 1
	sim_resbnd = run_residual
end bias 3


bias 4
	bias_bias = voltage
	bias_biasvalue = 4
	bias_electrode = 0
	bt_active = 502
	bt_biasnumb = 4
	n_steps = 20
	sim_Tamb = glob_temp
	sim_lstead = 1
	sim_resbnd = run_residual
end bias 4


bias 5
	bias_bias = voltage
	bias_biasvalue = 5
	bias_electrode = 0
	bt_active = 502
	bt_biasnumb = 5
	n_steps = 20
	sim_Tamb = glob_temp
	sim_lstead = 1
	sim_resbnd = run_residual
end bias 5


bias 6
	bias_bias = current
	bias_biasvalue = 100
	bias_electrode = 0
	bt_active = 502
	bt_biasnumb = 6
	n_steps = 50
	sim_Tamb = glob_temp
	sim_lim_newt = 50
	sim_lstead = 1
	sim_resbnd = run_residual
end bias 6


bias 7
	bias_bias = current
	bias_biasvalue = 200
	bias_electrode = 0
	bt_active = 500
	bt_biasnumb = 7
	n_steps = 100
	sim_Tamb = glob_temp
	sim_lim_newt = 50
	sim_lstead = 1
	sim_resbnd = run_residual
end bias 7


bias 8
	bias_bias = current
	bias_biasvalue = 300
	bias_electrode = 0
	bt_active = 502
	bt_biasnumb = 8
	n_steps = 50
	sim_Tamb = glob_temp
	sim_lim_newt = 50
	sim_lstead = 1
	sim_resbnd = run_residual
end bias 8


bias 9
	bias_bias = current
	bias_biasvalue = 400
	bias_electrode = 0
	bt_active = 502
	bt_biasnumb = 9
	n_steps = 50
	sim_Tamb = glob_temp
	sim_lim_newt = 50
	sim_lstead = 1
	sim_resbnd = run_residual
end bias 9


bias 10
	bias_bias = current
	bias_biasvalue = 500
	bias_electrode = 0
	bt_active = 502
	bt_biasnumb = 10
	n_steps = 50
	sim_Tamb = glob_temp
	sim_lim_newt = 50
	sim_lstead = 1
	sim_resbnd = run_residual
end bias 10
