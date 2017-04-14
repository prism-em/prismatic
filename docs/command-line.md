# Using PRISM from the command line

PRISM contains a command line tool, `prism`, that can be used to run simulations from within a terminal, bash script, etc. Building it requires the CMake variable `PRISM_ENABLE_CLI=1` at compilation time, which is the default behavior.

### Options
The following options are available with `prism`, each documented as **_long form_** **_(short form)_** *parameters* : description

* --**_input-file (-i)_** *filename* : the filename containing the atomic coordinates, which should be a plain text file with comma-separated values in the format x, y, z, Z 
* --**_output-file(-o)_** *filename* : output filename
* --**_interp-factor (-f)_** *number* : PRISM interpolation factor
* --**_num-threads (-j)_** *number* : number of CPU threads to use
* --**_num-streams (-S)_** *number* : number of CUDA streams to create per GPU
* --**_slice-thickness (-s)_** *thickness* : thickness of each slice of projected potential (in Angstroms)
* --**_num-gpus (-g)_** *number* : number of GPUs to use. A runtime check is performed to check how many are actually available, and minimum of these two numbers is used.
* --**_help(-h)_** : print information about the available options
* --**_pixel-size (-p)_** *pixel_size* : size of simulation pixel size
* --**_cell-dimension (-c)_** *x y z* : size of sample in x, y, z directions (in Angstroms)
* --**_algorithm (-a)_** *p/m* : the simulation algorithm to use, either (p)rism or (m)ultislice
* --**_energy (-E)_** *value* : the energy of the electron beam (in keV)
* --**_alpha-max (-A)_** *angle* : the maximum probe angle to consider (in mrad)
* --**_potential-bound (-P)_** *value* : the maximum radius from the center of each atom to compute the potental (in Angstroms)
* --**_also-do-cpu-work (-C)_** *0/1* : boolean value used to determine whether or not to also create CPU workers in addition to GPU ones
* --**_force-streaming-mode_** *0/1* : boolean value to force code to use (true) or not use (false) streaming versions of GPU codes. The default behavior is to estimate the needed memory from input parameters and choose automatically.
* --**_probe-step (-r)_** *step_size* : step size of the probe (in Angstroms)
* --**_num-FP (-F)_** *number* : number of frozen phonon configurations to calculate

