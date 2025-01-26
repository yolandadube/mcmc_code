import os, sys
#--------------------------------------------------------------------------------------------------------------------------

def pyCompileCLASS(mainPath):

	# create a '.dat' file to confirm that CLASS is well compiled
	pathDirectory = mainPath+'functions/pyCLASSWrapper/'

	if os.path.exists(pathDirectory+'compileSuccessful.dat'):
		print("----------------------------------------------------------------------")
	else:
		# compile CLASS
		commandLine1 = 'make clean'
		os.system('cd '+mainPath+'CLASSMultiTracer/'+'&&'+commandLine1)
		print("----------------------------------------------------------------------")
		commandLine2 = 'make class'
		os.system('cd '+mainPath+'CLASSMultiTracer/'+'&&'+commandLine2)
		print("----------------------------------------------------------------------")
		filename     = "compileSuccessful.dat"
		compile_     = open(pathDirectory+filename, 'w')
		compile_.write("%s"%('CLASS compilation successful!'))
		compile_.close()

	return None
#--------------------------------------------------------------------------------------------------------------------------