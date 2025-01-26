def pythonNumpyArrayToClassInputFormat(inputRedshiftArray):

	if (len(inputRedshiftArray)>1):
		# convert input numpy array item to string data type
		inputRedshiftStr = inputRedshiftArray.astype(str)
		# express 'inputRedshiftList' in format required for CLASS
		inputs  = tuple(inputRedshiftStr)
	else:
		inputs = str(inputRedshiftArray[0])

	return inputs
#----------------------------------------------------------------------

def ClassInputFormatToPythonList(inputList):

	# split string item of inputList by delimeter ','
	inputs = [i for item in inputList for i in item.split(',')]

	return inputs
#----------------------------------------------------------------------