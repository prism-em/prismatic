import pyprismatic.core
class Metadata(object):
	fields=[
	"interpolationFactorX",
	"interpolationFactorY",
	"filename_atoms",
	"filename_output",
	"realspace_pixelSizeX",
	"realspace_pixelSizeY",
	"potBound",
	"numFP",
	"sliceThickness",
	"cellDimX",
	"cellDimY",
	"cellDimZ",
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
	"integration_angle_min",
	"integration_angle_max",
	"transfer_mode"]
	def __init__(self, *args, **kwargs):
		import numpy as np
		self.interpolationFactorX 	  = 4
		self.interpolationFactorY 	  = 4
		self.filename_atoms		      = ""
		self.filename_output	      = "output.mrc"
		self.realspace_pixelSizeX     = 0.1
		self.realspace_pixelSizeY     = 0.1
		self.potBound				  = 0.1
		self.numFP 				      = 1
		self.sliceThickness		      = 2.0
		self.cellDimX                 = 20.0
		self.cellDimY                 = 20.0
		self.cellDimZ                 = 20.0
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
		# self.toString()
		pyprismatic.core.go(*[4,4,
		 "/home/aj/hdd1/clion/PRISM/SI100.XYZ", "/home/aj/hdd1/clion/PRISM/tt.mrc",
		 0.25, 0.30])

		# pyprismatic.core.go(*[getattr(self, field) for field in Metadata.fields])
		import numpy as np
		from pyprismatic.fileio import readMRC
		import matplotlib.pyplot as plt
		result = readMRC("/mnt/spareA/clion/PRISM/tt.mrc")
		plt.figure()
		plt.imshow(np.squeeze(np.sum(result,axis=2)))
		plt.show()
		# result = readMRC(meta.filename_output)
