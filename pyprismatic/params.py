import pyprismatic.core
class Metadata(object):
	fields=["interpolationFactorX", "interpolationFactorY"]
	def __init__(self, **kwargs):
		self.interpolationFactorX = 4
		self.interpolationFactorY = 4
		for k,v in kwargs:
			if k not in fields:
				print("Invalid metaparameter {} provided".format(k))
			else:
				self.setattr(k, v)
	def toString(self):
		print("interpolationFactorX = {}".format(self.interpolationFactorX))
		print("interpolationFactorY = {}".format(self.interpolationFactorY))
	def go(self):
		self.toString()
		pyprismatic.core.go()
		import numpy as np
		from pyprismatic.fileio import readMRC
		import matplotlib.pyplot as plt
		result = readMRC("/mnt/spareA/clion/PRISM/output_python.mrc")
		plt.figure()
		plt.imshow(np.squeeze(np.sum(result,axis=2)))
		plt.show()
		# result = readMRC(meta.filename_output)
