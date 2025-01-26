def latexDictionary1(mcmcParameter1):

	# mcmcParameter  ---> parameter name

	# parameter labels
	dictionary = {'A'    :r'$A$',
				  'k0'    :r'$k_{0}$', 
				  'n' :r'$n$', 
				  'm' :r'$m$',
				  'beta'  :r'$\beta$'
				 }

	return dictionary[mcmcParameter1]
#----------------------------------------------------------------------