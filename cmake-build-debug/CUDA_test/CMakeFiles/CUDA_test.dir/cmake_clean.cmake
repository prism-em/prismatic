file(REMOVE_RECURSE
  "CUDA_test.pdb"
  "CUDA_test"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/CUDA_test.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
