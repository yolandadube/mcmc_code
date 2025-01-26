def CLASSPrecisionParameters():

	# precision of CLASS
	k_step_trans_scalars                = 0.4
	selection_sampling_bessel           = 2.0 
	k_scalar_max_tau0_over_l_max        = 2.0 
	l_switch_limber_for_nc_local_over_z = 2000
	#---#
	precisionList = [k_step_trans_scalars, \
					selection_sampling_bessel, \
					k_scalar_max_tau0_over_l_max, \
					l_switch_limber_for_nc_local_over_z]
	
	return precisionList 
#----------------------------------------------------------------------