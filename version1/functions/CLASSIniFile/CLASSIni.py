import os, sys
from configobj import ConfigObj
#--------------------------------------------------------------------------------------------------------------------------

def iniFile(mainPath, cosmologY, gauge, output, nonLinear, z):

	sys.path.append(mainPath+"functions/CLASSInputFormat/")
	sys.path.append(mainPath+"functions/CLASSPrecisionParameters/")
	#---#
	from CLASSPrecision import CLASSPrecisionParameters
	from inputFormat import pythonNumpyArrayToClassInputFormat
	#---#

	# read .ini template file
	templateIniFile = mainPath+"CLASSMultiTracer/multiTracer_explanatory.ini" 
	config          = ConfigObj(templateIniFile, encoding='utf8')

	# fiducial Planck cosmological parameters
	A_s, sigma80, H0, h, omm0, omb0, omcdm0, omega_k0, \
	omega_n0, n_s, gamma, w, fnl = cosmologY

	# CLASS precision parameters
	k_step_trans_scalars, \
	selection_sampling_bessel, \
	k_scalar_max_tau0_over_l_max, \
	l_switch_limber_for_nc_local_over_z = CLASSPrecisionParameters() 

	# create .ini file for CLASS

	# format 'redshifts', 'binWidths', 'b1', 's' to CLASS input
	redshifts = pythonNumpyArrayToClassInputFormat(z)
	
	# configure .ini template file to accommodate changes according to 
	# selected tracer(s)
	config.filename = mainPath+"CLASSDataFile/iniDataFile/CREATED.ini"
	
	# CLASS precision
	config['k_step_trans_scalars']                = k_step_trans_scalars
	config['selection_sampling_bessel']           = selection_sampling_bessel
	config['k_scalar_max_tau0_over_l_max']        = k_scalar_max_tau0_over_l_max
	config['l_switch_limber_for_nc_los_over_z']   = ''
	config['l_switch_limber_for_nc_local_over_z'] = l_switch_limber_for_nc_local_over_z
	# 'mPk' for matter power spectrum, 'mTk' for transfer functions
	config['output']  = output
	# format for output .ini data file
	config['format']  = 'class'
	# headers for output .ini data file
	config['headers'] = 'no'
	# gauge in which calculations are performed
	config['gauge']   = gauge
	# primordial Helium fraction
	config['YHe']     = 'BBN'
	# set cosmological parameters
	config['A_s']          = A_s
	#config['sigma8']       = sigma80
	config['n_s']          = n_s
	config['h']            = h
	#config['H0']           = H0
	config['Omega_b']      = omb0
	#config['omega_b']     = 
	config['Omega_cdm']    = omcdm0
	#config['omega_cdm']   = 
	config['Omega_k']      = omega_k0
	config['w0_fld']       = w
	# primordial non-Gaussianity parameter
	config['f_NL']         = fnl
	# do non-linear P(k) using 'halofit' emulator 
	if (nonLinear==True):	
		config['non linear']   = 'halofit'
	else:
		config['non linear']   = ''
	# scalar modes
	config['modes']            = 's'
	# redshift bin centres
	config['z_pk']             = redshifts
	# linear bias tracer (mainly used for Cls)
	config['selection_bias']   = ''
	#---#
	# leave unselected/empty values as blank
	config.write_empty_values  = True 
	# write updated .ini data file
	config.write()

	# copy .ini file to 'CLASSMultiTracer' folder
	os.system("cp "+config.filename+" "+mainPath+"CLASSMultiTracer/")
	
	return None
#--------------------------------------------------------------------------------------------------------------------------