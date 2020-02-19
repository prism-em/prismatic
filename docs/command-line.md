# Using PRISM from the command line

PRISM contains a command line tool, `prism`, that can be used to run simulations from within a terminal, bash script, etc. Building it requires the CMake variable `PRISM_ENABLE_CLI=1` at compilation time, which is the default behavior.

### Options

The following options are available with `prism`, each documented as **_long form_** **_(short form)_** _parameters_ : description

- --**_input-file (-i)_** _filename_ : the filename containing the atomic coordinates, which should be a plain text file with comma-separated values in the format x, y, z, Z
- --**_output-file(-o)_** _filename_ : output filename
- --**_interp-factor (-f)_** _number_ : PRISM interpolation factor
- --**_num-threads (-j)_** _number_ : number of CPU threads to use
- --**_num-streams (-S)_** _number_ : number of CUDA streams to create per GPU
- --**_slice-thickness (-s)_** _thickness_ : thickness of each slice of projected potential (in Angstroms)
- --**_num-gpus (-g)_** _number_ : number of GPUs to use. A runtime check is performed to check how many are actually available, and minimum of these two numbers is used.
- --**_help(-h)_** : print information about the available options
- --**_pixel-size (-p)_** _pixel_size_ : size of simulation pixel size
- --**_cell-dimension (-c)_** _x y z_ : size of sample in x, y, z directions (in Angstroms)
- --**_algorithm (-a)_** _p/m_ : the simulation algorithm to use, either (p)rism or (m)ultislice
- --**_energy (-E)_** _value_ : the energy of the electron beam (in keV)
- --**_alpha-max (-A)_** _angle_ : the maximum probe angle to consider (in mrad)
- --**_potential-bound (-P)_** _value_ : the maximum radius from the center of each atom to compute the potental (in Angstroms)
- --**_also-do-cpu-work (-C)_** _0/1_ : boolean value used to determine whether or not to also create CPU workers in addition to GPU ones
- --**_force-streaming-mode_** _0/1_ : boolean value to force code to use (true) or not use (false) streaming versions of GPU codes. The default behavior is to estimate the needed memory from input parameters and choose automatically.
- --**_probe-step (-r)_** _step_size_ : step size of the probe (in Angstroms)
- --**_num-FP (-F)_** _number_ : number of frozen phonon configurations to calculate
