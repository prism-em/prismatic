# This script calculates multiple images using a unit cell tiled varying numbers of times in Z. 
# This type of simulation is useful for determining the true thickness of a sample from experimental images

import pyprismatic as pr
import numpy as np
import matplotlib.pyplot as plt

base_output_name = "thickness_scan"
base_output_ext  = ".mrc"
atom_filename    = "../SI100.XYZ"
tileX=tileY=3
meta = pr.Metadata(filenameAtoms=atom_filename, tileX=tileX, tileY=tileY)
output_filenames=[]
for tileZ in range(1,17):
	meta.tileZ = tileZ
	meta.filenameOutput = base_output_name + str(tileZ) + base_output_ext
	output_filenames.extend([meta.filenameOutput])
	meta.go()

f, ax = plt.subplots(4, 4)
ax = ax.ravel()
for i, filename in enumerate(output_filenames):
	stack = pr.fileio.readMRC(filename)
	img   = np.sum(stack[:, :, :10], axis=2)
	ax[i].imshow(img)
	ax[i].set_title("Thicknesss: {:.2f}$\AA$".format((i + 1) * 5.43))
	ax[i].set_yticklabels([])
	ax[i].set_xticklabels([])
plt.suptitle("Bright Field (0-10 mrad) PRISM images for SI100 of varying thicknesses")
plt.show()
