# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/cmake-gui

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/toe/Documents/Recherche/HSP/HSPArborescent

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/toe/Documents/Recherche/HSP/HSPArborescent

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/cmake-gui -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/toe/Documents/Recherche/HSP/HSPArborescent/CMakeFiles /home/toe/Documents/Recherche/HSP/HSPArborescent/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/toe/Documents/Recherche/HSP/HSPArborescent/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named cascadeEA

# Build rule for target.
cascadeEA: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 cascadeEA
.PHONY : cascadeEA

# fast build rule for target.
cascadeEA/fast:
	$(MAKE) -f CMakeFiles/cascadeEA.dir/build.make CMakeFiles/cascadeEA.dir/build
.PHONY : cascadeEA/fast

#=============================================================================
# Target rules for targets named realisable

# Build rule for target.
realisable: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 realisable
.PHONY : realisable

# fast build rule for target.
realisable/fast:
	$(MAKE) -f CMakeFiles/realisable.dir/build.make CMakeFiles/realisable.dir/build
.PHONY : realisable/fast

#=============================================================================
# Target rules for targets named testRead

# Build rule for target.
testRead: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 testRead
.PHONY : testRead

# fast build rule for target.
testRead/fast:
	$(MAKE) -f CMakeFiles/testRead.dir/build.make CMakeFiles/testRead.dir/build
.PHONY : testRead/fast

# target to build an object file
cascadeEA.o:
	$(MAKE) -f CMakeFiles/cascadeEA.dir/build.make CMakeFiles/cascadeEA.dir/cascadeEA.o
.PHONY : cascadeEA.o

# target to preprocess a source file
cascadeEA.i:
	$(MAKE) -f CMakeFiles/cascadeEA.dir/build.make CMakeFiles/cascadeEA.dir/cascadeEA.i
.PHONY : cascadeEA.i

# target to generate assembly for a file
cascadeEA.s:
	$(MAKE) -f CMakeFiles/cascadeEA.dir/build.make CMakeFiles/cascadeEA.dir/cascadeEA.s
.PHONY : cascadeEA.s

# target to build an object file
realisable.o:
	$(MAKE) -f CMakeFiles/realisable.dir/build.make CMakeFiles/realisable.dir/realisable.o
.PHONY : realisable.o

# target to preprocess a source file
realisable.i:
	$(MAKE) -f CMakeFiles/realisable.dir/build.make CMakeFiles/realisable.dir/realisable.i
.PHONY : realisable.i

# target to generate assembly for a file
realisable.s:
	$(MAKE) -f CMakeFiles/realisable.dir/build.make CMakeFiles/realisable.dir/realisable.s
.PHONY : realisable.s

# target to build an object file
testRead.o:
	$(MAKE) -f CMakeFiles/testRead.dir/build.make CMakeFiles/testRead.dir/testRead.o
.PHONY : testRead.o

# target to preprocess a source file
testRead.i:
	$(MAKE) -f CMakeFiles/testRead.dir/build.make CMakeFiles/testRead.dir/testRead.i
.PHONY : testRead.i

# target to generate assembly for a file
testRead.s:
	$(MAKE) -f CMakeFiles/testRead.dir/build.make CMakeFiles/testRead.dir/testRead.s
.PHONY : testRead.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... cascadeEA"
	@echo "... edit_cache"
	@echo "... realisable"
	@echo "... rebuild_cache"
	@echo "... testRead"
	@echo "... cascadeEA.o"
	@echo "... cascadeEA.i"
	@echo "... cascadeEA.s"
	@echo "... realisable.o"
	@echo "... realisable.i"
	@echo "... realisable.s"
	@echo "... testRead.o"
	@echo "... testRead.i"
	@echo "... testRead.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

