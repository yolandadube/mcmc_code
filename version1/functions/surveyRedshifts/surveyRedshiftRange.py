def redshiftRange(survey, specsSurvey):

	if (survey=="galaxy"):

		zMin = zMax = None

		# H\alpha euclid-like
		if (specsSurvey=="HAlpha"):
			# minimum redshift
			zMin = 0.9
			# maximum redshift
			zMax = 1.8

		# DESI BGS 
		if (specsSurvey=="desiBGS"):
			# minimum redshift
			zMin = 0.0
			# maximum redshift
			zMax = 0.5

		# DESI ELG 
		if (specsSurvey=="desiELG"):
			# minimum redshift
			zMin = 0.6
			# maximum redshift
			zMax = 1.7

		# skaOGalBand2
		if (specsSurvey=="skaOGalBand2"):
			# minimum redshift
			zMin = 0.0
			# maximum redshift
			zMax = 0.5

		# MegaMapper LBG
		if (specsSurvey=="megaMapLBG"):
			# minimum redshift
			zMin = 2.0
			# maximum redshift
			zMax = 5.0

	if (survey=="HI IM"):

		# ska1 Band 1
		if (specsSurvey=="skaOIMBand1"):
			# minimum redshift
			zMin = 0.35
			# maximum redshift
			zMax = 3.05

		# ska1 Band 2
		if (specsSurvey=="skaOIMBand2"):
			# minimum redshift
			zMin = 0.10
			# maximum redshift
			zMax = 0.58

		# meerKAT L Band
		if (specsSurvey=="meerKATLBand"):
			# minimum redshift
			zMin = 0.15
			# maximum redshift
			zMax = 0.53

		# meerKAT UHF Band
		if (specsSurvey=="meerKATUHFBand"):
			# minimum redshift
			zMin = 0.45
			# maximum redshift
			zMax = 1.40
		
		'''
		# chime
		elif (specsSurvey=="chime"):
			# minimum redshift
			zMin = 2.50
			# maximum redshift
			zMax = 6.00
		'''
		
		# hirax256
		if (specsSurvey=="hirax256"):
			# minimum redshift
			zMin = 0.8
			# maximum redshift
			zMax = 2.5

		# hirax1024
		if (specsSurvey=="hirax1024"):
			# minimum redshift
			zMin = 0.8
			# maximum redshift
			zMax = 2.5

		# puma5k
		if (specsSurvey=="puma5k"):
			# minimum redshift
			zMin = 2.0
			# maximum redshift
			zMax = 5.0 #6.0

		# puma32k
		if (specsSurvey=="puma32k"):
			# minimum redshift
			zMin = 2.0
			# maximum redshift
			zMax = 5.0 #6.0

		# DSA2000
		if (specsSurvey=="dsa2000"):
			# minimum redshift
			zMin = 0.0
			# maximum redshift
			zMax = 1.0 
	
	return zMin, zMax
#----------------------------------------------------------------------