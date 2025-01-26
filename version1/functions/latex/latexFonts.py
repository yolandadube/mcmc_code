def latexDictionary(mcmcParameter):

	# mcmcParameter  ---> parameter name

	# parameter labels
	dictionary = {'p0'    :r'$P_{0}$',
				  'k0'    :r'$k_{0}$', 
				  'alpha' :r'$\alpha$', 
				  'beta'  :r'$\beta$'
				 }

	return dictionary[mcmcParameter]
#----------------------------------------------------------------------