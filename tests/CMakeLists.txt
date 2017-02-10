add_subdirectory(Array2D_test)

if (PRISM_USE_GPU)
    add_subdirectory(CUDA_test)
endif(PRISM_USE_GPU)
