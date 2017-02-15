# PRISM (MATLAB)

PRISM02_ajp.m is the same as PRISM02.m except that the bottleneck fft2 line
is modified to use a single precision GPU array.

PRISM03_ajp.m is MEX accelerated and requires compilation of PRISM03_mex.cpp. See
`HOWTO_compile_mex.txt` for information about the actual compilation call. You can adjust
the number of threads used with the NUM_THREADS macro defined in `PRISM03.cpp`

It is a little tricky to compile CUDA code and link it with MATLAB. In principle
there is `mexcuda`, but it doesn't seem to be compatible with the latest version
of the CUDA toolkit. A solution is to compile cuda code with `nvcc` and then 
compile/link with MEX. I documented the commands I used to do this in `HOWTO_compile_cuda.txt`
using a text kernel `timesTwo.cu`.

