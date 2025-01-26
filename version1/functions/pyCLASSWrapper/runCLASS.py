import os, sys
#--------------------------------------------------------------------------------------------------------------------------

def pyCLASS(mainPath, cosmologY, gauge, nonLinear, z):
	
	sys.path.append(mainPath+"functions/CLASSIniFile/")
	#---#
	from CLASSIni import iniFile
	#---#

	# empty 'output' folder from CLASSDataFile
	os.system("rm -rf "+mainPath+"CLASSDataFile/outputs/CREATED*")

	# outputs from CLASS
	# 'mPk' for matter power spectrum, 'mTk' for transfer functions
	outputs = ['mPk', 'mTk']

	# create .ini file for CLASS
	for i in range(len(outputs)):
		# .ini file for matter power spectrum at any z
		if (outputs[i]=='mPk'):
			classIniFile = iniFile(mainPath, cosmologY, gauge, outputs[i], nonLinear, z)
		# .ini file for transfer function at z=0
		elif (outputs[i]=='mTk'):
			classIniFile = iniFile(mainPath, cosmologY, gauge, outputs[i], nonLinear, [0.0])
	
		# run CLASS to compute 'outputs' for the created .ini file
		commandLine = './class '+'CREATED.ini'
		os.system('cd '+mainPath+'CLASSMultiTracer/'+'&&'+commandLine)

		# create a directory to store files produced by CLASS
		pathDirectory = mainPath+'CLASSDataFile/outputs/'
		if os.path.exists(pathDirectory):
			# if 'pathDirectory' exists, empty it
			os.system("rm -rf "+pathDirectory+"CREATED*")
		else:
			# if not, then create directory
			os.mkdir(pathDirectory)

		# copy all files produced in 'output' folder from CLASS to created directory
		os.system("cp "+mainPath+"CLASSMultiTracer/output/*.dat"+" "+pathDirectory)
		print("----------------------------------------------------------------------")

	# empty 'output' folder from CLASS
	os.system("rm -rf "+mainPath+"CLASSMultiTracer/output/CREATED*")
	# remove all 'CREATED.ini' from CLASS
	os.system("rm -rf "+mainPath+"CLASSMultiTracer/CREATED*")

	return None
#--------------------------------------------------------------------------------------------------------------------------