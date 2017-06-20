import pyprismatic.core
class Metadata(object):
	fields=[
	"interpolationFactorX",
	"interpolationFactorY",
	"filename_atoms",
	"filename_output",
	"realspace_pixelSize",
	"potBound",
	"numFP",
    "fpNum",
	"sliceThickness",
	"cellDim",
	"tileX",
	"tileY",
	"tileZ",
	"E0",
	"alphaBeamMax",
	"NUM_GPUS",
	"NUM_STREAMS_PER_GPU",
	"NUM_THREADS",
	"batch_size_target_CPU",
	"batch_size_target_GPU",
	"batch_size_CPU",
	"batch_size_GPU",
	"gpu_cpu_ratio",
    "probe_stepX",
	"probe_stepY",
	"probeDefocus",
	"C3",
	"C5",
	"probeSemiangle",
	"detector_angle_step",
	"probeXtilt",
	"probeYtilt",
	"scanWindowXMin",
	"scanWindowXMax",
	"scanWindowYMin",
	"scanWindowYMax",
	"random_seed",
	"algorithm",
	"include_thermal_effects",
	"also_do_CPU_work",
	"save2DOutput",
	"save3DOutput",
	"save4DOutput",
	"user_specified_celldims",
	"integration_angle_min",
	"integration_angle_max",
	"transfer_mode"]
	def __init__(self, *args, **kwargs):
		import numpy as np
		self.interpolationFactorX 	  = 4
		self.interpolationFactorY 	  = 4
		self.filename_atoms		      = ""
		self.filename_output	      = "output.mrc"
		self.realspace_pixelSize      = (0.1, 0.1)
		self.potBound				  = 0.1
		self.numFP 				      = 1
		self._fpNum					  = 1
		self.sliceThickness		      = 2.0
		self.cellDim                  = (20.0, 20.0, 20.0)
		self.tileX					  = 3
		self.tileY					  = 3
		self.tileZ					  = 1
		self.E0						  = 80e3
		self.alphaBeamMax			  = 0.024
		self.NUM_GPUS				  = 4
		self.NUM_STREAMS_PER_GPU      = 3	
		self.NUM_THREADS		      = 12
		self.batch_size_target_CPU	  = 1
		self.batch_size_target_GPU    = 1
		self.batch_size_CPU           = 1
		self.batch_size_GPU           = 1
		self.gpu_cpu_ratio            = 100
		self.probe_stepX		      = 0.25
		self.probe_stepY			  = 0.25
		self.probeDefocus			  = 0.0
		self.C3						  = 0.0
		self.C5						  = 0.0
		self.probeSemiangle			  = 0.02
		self.detector_angle_step	  = 0.01
		self.probeXtilt				  = 0.0
		self.probeYtilt				  = 0.0
		self.scanWindowXMin			  = 0.0
		self.scanWindowXMax			  = 1.0
		self.scanWindowYMin  		  = 0.0
		self.scanWindowYMax			  = 1.0
		self.random_seed		      = np.random.randint(0,999999)
		self.algorithm				  = "prism"
		self.include_thermal_effects  = True
		self.also_do_CPU_work		  = True
		self.save2DOutput			  = False
		self.save3DOutput		      = True
		self.save4DOutput			  = False
		self._user_specified_celldims = False
		self.integration_angle_min    = 0
		self.integration_angle_max    = self.detector_angle_step
		self.transfer_mode		      = "auto"
		for k,v in kwargs.items():
			if k not in Metadata.fields:
				print("Invalid metaparameter \"{}\" provided".format(k))
			else:
				setattr(self, k, v)
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
