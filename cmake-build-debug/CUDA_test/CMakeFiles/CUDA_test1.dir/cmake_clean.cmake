file(REMOVE_RECURSE
  "CUDA_test1.pdb"
  "CUDA_test1"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/CUDA_test1.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
