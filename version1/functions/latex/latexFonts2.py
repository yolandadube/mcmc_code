def latexDictionary2(mcmcParameter2):

	# mcmcParameter  ---> parameter name

	# parameter labels
	dictionary = {'A'    :r'$A$',
				  'k0'    :r'$k_{0}$', 
				  'n' :r'$n$', 
				  'm' :r'$m$',
				  'Delta'  :r'$\Delta$'
				 }

	return dictionary[mcmcParameter2]
#----------------------------------------------------------------------