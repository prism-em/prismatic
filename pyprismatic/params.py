# Copyright Alan (AJ) Pryor, Jr. 2017
# Transcribed from MATLAB code by Colin Ophus
# Prismatic is distributed under the GNU General Public License (GPL)
# If you use Prismatic, we kindly ask that you cite the following papers:

# 1. Ophus, C.: A fast image simulation algorithm for scanning
#    transmission electron microscopy. Advanced Structural and
#    Chemical Imaging 3(1), 13 (2017)

# 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
#    Implementation of Image Simulation Algorithms for Scanning
#    Transmission Electron Microscopy. arXiv:1706.08563 (2017)

from typing import List, Any
import pyprismatic.core


class Metadata:
    """
    "interpolationFactorX" : PRISM interpolation factor in x-direction
    "interpolationFactorY" : PRISM interpolation factor in y-direction
    "filenameAtoms" : filename containing input atom information in XYZ format (see http://prism-em.com/about/ for more details)
    "filenameOutput" : filename in which to save the 3D output. Also serves as base filename for 2D and 4D outputs if used
    "realspacePixelSizeX" : size of pixel size in X for probe/potential arrays
    "realspacePixelSizeY" : size of pixel size in Y for probe/potential arrays
    "potBound" : limiting radius within which to compute projected potentials from the center of each atom (in Angstroms)
    "numFP" : number of frozen phonon configurations to average over
    "sliceThickness" : thickness of potential slices (in Angstroms)
    "numSlices" : number of slices between intermediate outputs
    "zStart" : depth before intermediate output begins (in Angstroms)
    "cellDimX" : unit cell dimension X (in Angstroms)
    "cellDimY" : unit cell dimension Y (in Angstroms)
    "cellDimZ" : unit cell dimension Z (in Angstroms)
    "tileX" : number of unit cells to tile in X direction
    "tileY" : number of unit cells to tile in Y direction
    "tileZ" : number of unit cells to tile in Z direction
    "E0" : electron beam energy (in eV)
    "alphaBeamMax" : the maximum probe angle to consider (in rad)
    "numGPUs" : number of GPUs to use. A runtime check is performed to check how many are actually available, and the minimum of these two numbers is used.
    "numStreamsPerGPU" : number of CUDA streams to use per GPU
    "numThreads" : number of CPU worker threads to use
    "batchSizeTargetCPU" : desired batch size for CPU FFTs.
    "batchSizeTargetGPU" : desired batch size for GPU FFTs.
    "earlyCPUStopCount" : the WorkDispatcher will cease providing work to CPU workers earlyCPUStopCount jobs from the end. This is to prevent the program waiting for slower CPU workers to complete.
    "probeStepX" : step size of the probe in X direction (in Angstroms)
    "probeStepY" : step size of the probe in Y direction (in Angstroms)
    "probeDefocus" : probe defocus (in Angstroms)
    "C3" : microscope C3 (in Angstroms)
    "C5" : microscope C5 (in Angstroms)
    "probeSemiangle" : probe convergence semi-angle (in rad)
    "detectorAngleStep" : angular step size for detector integration bins (in rad)
    "probeXtilt" : (in Angstroms)
    "probeYtilt" : (in Angstroms)
    "scanWindowXMin" : lower X size of the window to scan the probe (in fractional coordinates)
    "scanWindowXMax" : upper X size of the window to scan the probe (in fractional coordinates)
    "scanWindowYMin" : lower Y size of the window to scan the probe (in fractional coordinates)
    "scanWindowYMax" : upper Y size of the window to scan the probe (in fractional coordinates)
    "scanWindowXMin_r" : lower X size of the window to scan the probe (in Angstroms)
    "scanWindowYMin_r" : lower Y size of the window to scan the probe (in Angstroms)
    "scanWindowXMax_r" : upper X size of the window to scan the probe (in Angstroms)
    "scanWindowYMax_r" : upper Y size of the window to scan the probe (in Angstroms)
    "randomSeed" : number to use for random seeding of thermal effects
    "algorithm" : simulation algorithm to use, "prism" or "multislice"
    "includeThermalEffects" : true/false to apply random thermal displacements (Debye-Waller effect)
    "alsoDoCPUWork" : true/false
    "save2DOutput" : save the 2D STEM image integrated between integrationAngleMin and integrationAngleMax
    "save3DOutput" : true/false Also save the 3D output at the detector for each probe (3D output mode)
    "save4DOutput" : true/false Also save the 4D output at the detector for each probe (4D output mode)
    "saveDPC_CoM"  : true/false Also save the DPC center of mass calculation for each probe
    "savePotentialSlices" : true/false Also save the projected potential array
    "crop4DOutput" : true/false Crop the 4D output smaller than the anti-aliasing boundary (default: False)
    "crop4Damax" : float If crop4D, the maximum angle to which the output is cropped (in mrad) (default: 100)
    "nyquistSampling": set number of probe positions at Nyquist sampling limit
    "integrationAngleMin" : (in rad)
    "integrationAngleMax" : (in rad)
    "transferMode" : memory model to use, either "streaming", "singlexfer", or "auto"
    """

    fields: List[str] = [
        "interpolationFactorX",
        "interpolationFactorY",
        "filenameAtoms",
        "filenameOutput",
        "realspacePixelSizeX",
        "realspacePixelSizeY",
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
        "numGPUs",
        "numStreamsPerGPU",
        "numThreads",
        "batchSizeTargetCPU",
        "batchSizeTargetGPU",
        "earlyCPUStopCount",
        "probeStepX",
        "probeStepY",
        "probeDefocus",
        "C3",
        "C5",
        "probeSemiangle",
        "detectorAngleStep",
        "probeXtilt",
        "probeYtilt",
        "scanWindowXMin",
        "scanWindowXMax",
        "scanWindowYMin",
        "scanWindowYMax",
        "randomSeed",
        "algorithm",
        "includeThermalEffects",
        "alsoDoCPUWork",
        "save2DOutput",
        "save3DOutput",
        "save4DOutput",
        "integrationAngleMin",
        "integrationAngleMax",
        "transferMode",
        "saveDPC_CoM",
        "savePotentialSlices",
        "crop4DOutput",
        "crop4Damax",
        "nyquistSampling",
        "numSlices",
        "zStart",
        "scanWindowXMin_r",
        "scanWindowXMax_r",
        "scanWindowYMin_r",
        "scanWindowYMax_r",
    ]

    str_fields: List[str] = ["algorithm", "transferMode"]

    int_fields: List[str] = [
        "interpolationFactorX",
        "interpolationFactorY",
        "tileX",
        "tileY",
        "tileZ",
        "numFP",
        "numGPUs",
        "numStreamsPerGPU",
        "numThreads",
        "batchSizeTargetCPU",
        "batchSizeTargetGPU",
        "batchSizeCPU",
        "batchSizeGPU",
        "numSlices",
    ]

    float_fields: List[str] = [
        "realspacePixelSizeX",
        "realspacePixelSizeY",
        "potBound",
        "sliceThickness",
        "E0",
        "alphaBeamMax" "earlyCPUStopCount",
        "probeStepX",
        "probeStepY",
        "probeDefocus",
        "C3",
        "C5",
        "probeSemiangle",
        "detectorAngleStep",
        "probeXtilt",
        "probeYtilt",
        "scanWindowXMin",
        "scanWindowXMax",
        "scanWindowYMin",
        "scanWindowYMax",
        "scanWindowXMin_r",
        "scanWindowXMax_r",
        "scanWindowYMin_r",
        "scanWindowYMax_r",
        "integrationAngleMin",
        "integrationAngleMax",
        "zStart",
        "crop4Damax",
    ]

    def __init__(self, *args, **kwargs):
        """
        Fields within Metadata objects can be set either manually or at
        construction time. For example:
        meta = Metadata(interpolationFactorX=8, filenameOutput="test.mrc")
        would initialize a Metadata object with all parameters set to defaults
        except for interpolationFactorX and filenameOutput, which take on the
        values 8 and "test.mrc"
        """
        import numpy as np

        self.interpolationFactorX = 4
        self.interpolationFactorY = 4
        self.filenameAtoms = ""
        self.filenameOutput = "output.h5"
        self.realspacePixelSizeX = 0.1
        self.realspacePixelSizeY = 0.1
        self.potBound = 1.0
        self.numFP = 1
        self.sliceThickness = 2.0
        self.numSlices = 0
        self.zStart = 0.0
        self.cellDimX = 20.0
        self.cellDimY = 20.0
        self.cellDimZ = 20.0
        self.tileX = 3
        self.tileY = 3
        self.tileZ = 1
        self.E0 = 80e3
        self.alphaBeamMax = 0.024
        self.numGPUs = 4
        self.numStreamsPerGPU = 3
        self.numThreads = 12
        self.batchSizeTargetCPU = 1
        self.batchSizeTargetGPU = 2
        self.batchSizeCPU = 1
        self.batchSizeGPU = 1
        self.earlyCPUStopCount = 100.0
        self.probeStepX = 0.25
        self.probeStepY = 0.25
        self.probeDefocus = 0.0
        self.C3 = 0.0
        self.C5 = 0.0
        self.probeSemiangle = 0.02
        self.detectorAngleStep = 0.001
        self.probeXtilt = 0.0
        self.probeYtilt = 0.0
        self.scanWindowXMin = 0.0
        self.scanWindowXMax = 0.99999
        self.scanWindowYMin = 0.0
        self.scanWindowYMax = 0.99999
        self.scanWindowXMin_r = 0.0
        self.scanWindowXMax_r = 1.0
        self.scanWindowYMin_r = 0.0
        self.scanWindowYMax_r = 1.0
        self.randomSeed = np.random.randint(0, 999999)
        self.algorithm = "prism"
        self.includeThermalEffects = False
        self.alsoDoCPUWork = True
        self.save2DOutput = False
        self.save3DOutput = True
        self.save4DOutput = False
        self.saveDPC_CoM = False
        self.crop4DOutput = False
        self.crop4Damax = 0.1
        self.savePotentialSlices = False
        self.nyquistSampling = False
        self.integrationAngleMin = 0
        self.integrationAngleMax = 0.001
        self.transferMode = "auto"
        for k, v in kwargs.items():
            if k not in Metadata.fields:
                print('Invalid metaparameter "{}" provided'.format(k))
            else:
                setattr(self, k, v)

    def _setCellDims(self, filename: str):
        try:
            inf = open(filename, "r")
        except IOError:
            print("Could not set cell dimensions from file {}".format(filename))
            return
        inf.readline()
        self.cellDimX, self.cellDimY, self.cellDimZ = [
            float(i) for i in inf.readline().split()
        ]
        inf.close()

    @property
    def filenameAtoms(self):
        return self._filenameAtoms

    @filenameAtoms.setter
    def filenameAtoms(self, filenameAtoms):
        if filenameAtoms != "":  # do not set cell dimensions for default empty string
            self._filenameAtoms = filenameAtoms
            self._setCellDims(filenameAtoms)

    def readParameters(self, filename: str):
        """Read parameters from ``filename`` previously stored by ``writeParameters()``.
        No input verification is performed.
        """
        try:
            inf = open(filename, "r")
        except IOError:
            print("Could not open parameter file {}".format(filename))

        line = inf.readline()
        while line:
            field, value = line.split(" = ")
            if field in Metadata.str_fields:
                setattr(self, field, value)
            elif field in Metadata.int_fields:
                setattr(self, field, int(value))
            elif field in Metadata.float_fields:
                setattr(self, field, float(value))
            else:
                setattr(self, field, bool(value))
            line = inf.readline()
        inf.close()

    def writeParameters(self, filename: str):
        """Write parameters to ``filename`` but leave out parameters that define
        the input obtained from ``filenameAtoms`` (incuding cell dimensions
        and number of tiles) as well as the output specific settings.
        """
        try:
            outf = open(filename, "w")
        except IOError:
            print("Could not open parameter file {}".format(filename))

        for field in Metadata.fields:
            if (
                ("save" not in field)
                and ("filename" not in field)
                and ("cellDim" not in field)
            ):
                #
                # only save parameters that define a calculation but
                # omit all related to specific input or to what kind of output is
                # generated
                #
                outf.write("{} = {}\n".format(field, getattr(self, field)))
        outf.close()

    def toString(self):
        for field in Metadata.fields:
            print("{} = {}".format(field, getattr(self, field)))

    def go(self):
        self.algorithm: str = self.algorithm.lower()
        self.transferMode: str = self.transferMode.lower()
        l: List[Any] = [getattr(self, field) for field in Metadata.fields]
        pyprismatic.core.go(*l)
