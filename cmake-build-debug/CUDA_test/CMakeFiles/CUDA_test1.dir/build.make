# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/aj/clion-2016.3.2/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/aj/clion-2016.3.2/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/aj/hdd1/clion/PRISM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/aj/hdd1/clion/PRISM/cmake-build-debug

# Include any dependencies generated for this target.
include CUDA_test/CMakeFiles/CUDA_test1.dir/depend.make

# Include the progress variables for this target.
include CUDA_test/CMakeFiles/CUDA_test1.dir/progress.make

# Include the compile flags for this target's objects.
include CUDA_test/CMakeFiles/CUDA_test1.dir/flags.make

# Object files for target CUDA_test1
CUDA_test1_OBJECTS =

# External object files for target CUDA_test1
CUDA_test1_EXTERNAL_OBJECTS =

CUDA_test/CUDA_test1: CUDA_test/CMakeFiles/CUDA_test1.dir/build.make
CUDA_test/CUDA_test1: CUDA_test/CMakeFiles/CUDA_test1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/aj/hdd1/clion/PRISM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX executable CUDA_test1"
	cd /home/aj/hdd1/clion/PRISM/cmake-build-debug/CUDA_test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CUDA_test1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CUDA_test/CMakeFiles/CUDA_test1.dir/build: CUDA_test/CUDA_test1

.PHONY : CUDA_test/CMakeFiles/CUDA_test1.dir/build

CUDA_test/CMakeFiles/CUDA_test1.dir/requires:

.PHONY : CUDA_test/CMakeFiles/CUDA_test1.dir/requires

CUDA_test/CMakeFiles/CUDA_test1.dir/clean:
	cd /home/aj/hdd1/clion/PRISM/cmake-build-debug/CUDA_test && $(CMAKE_COMMAND) -P CMakeFiles/CUDA_test1.dir/cmake_clean.cmake
.PHONY : CUDA_test/CMakeFiles/CUDA_test1.dir/clean

CUDA_test/CMakeFiles/CUDA_test1.dir/depend:
	cd /home/aj/hdd1/clion/PRISM/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/aj/hdd1/clion/PRISM /home/aj/hdd1/clion/PRISM/CUDA_test /home/aj/hdd1/clion/PRISM/cmake-build-debug /home/aj/hdd1/clion/PRISM/cmake-build-debug/CUDA_test /home/aj/hdd1/clion/PRISM/cmake-build-debug/CUDA_test/CMakeFiles/CUDA_test1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CUDA_test/CMakeFiles/CUDA_test1.dir/depend
