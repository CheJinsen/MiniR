# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/a008/coding/c++/MiniR

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/a008/coding/c++/MiniR/build

# Include any dependencies generated for this target.
include src/CMakeFiles/stat.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/stat.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/stat.dir/flags.make

src/CMakeFiles/stat.dir/main.cpp.o: src/CMakeFiles/stat.dir/flags.make
src/CMakeFiles/stat.dir/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/a008/coding/c++/MiniR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/stat.dir/main.cpp.o"
	cd /home/a008/coding/c++/MiniR/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/stat.dir/main.cpp.o -c /home/a008/coding/c++/MiniR/src/main.cpp

src/CMakeFiles/stat.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/stat.dir/main.cpp.i"
	cd /home/a008/coding/c++/MiniR/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/a008/coding/c++/MiniR/src/main.cpp > CMakeFiles/stat.dir/main.cpp.i

src/CMakeFiles/stat.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/stat.dir/main.cpp.s"
	cd /home/a008/coding/c++/MiniR/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/a008/coding/c++/MiniR/src/main.cpp -o CMakeFiles/stat.dir/main.cpp.s

src/CMakeFiles/stat.dir/normal.cpp.o: src/CMakeFiles/stat.dir/flags.make
src/CMakeFiles/stat.dir/normal.cpp.o: ../src/normal.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/a008/coding/c++/MiniR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/stat.dir/normal.cpp.o"
	cd /home/a008/coding/c++/MiniR/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/stat.dir/normal.cpp.o -c /home/a008/coding/c++/MiniR/src/normal.cpp

src/CMakeFiles/stat.dir/normal.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/stat.dir/normal.cpp.i"
	cd /home/a008/coding/c++/MiniR/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/a008/coding/c++/MiniR/src/normal.cpp > CMakeFiles/stat.dir/normal.cpp.i

src/CMakeFiles/stat.dir/normal.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/stat.dir/normal.cpp.s"
	cd /home/a008/coding/c++/MiniR/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/a008/coding/c++/MiniR/src/normal.cpp -o CMakeFiles/stat.dir/normal.cpp.s

# Object files for target stat
stat_OBJECTS = \
"CMakeFiles/stat.dir/main.cpp.o" \
"CMakeFiles/stat.dir/normal.cpp.o"

# External object files for target stat
stat_EXTERNAL_OBJECTS =

../bin/stat: src/CMakeFiles/stat.dir/main.cpp.o
../bin/stat: src/CMakeFiles/stat.dir/normal.cpp.o
../bin/stat: src/CMakeFiles/stat.dir/build.make
../bin/stat: src/CMakeFiles/stat.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/a008/coding/c++/MiniR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ../../bin/stat"
	cd /home/a008/coding/c++/MiniR/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/stat.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/stat.dir/build: ../bin/stat

.PHONY : src/CMakeFiles/stat.dir/build

src/CMakeFiles/stat.dir/clean:
	cd /home/a008/coding/c++/MiniR/build/src && $(CMAKE_COMMAND) -P CMakeFiles/stat.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/stat.dir/clean

src/CMakeFiles/stat.dir/depend:
	cd /home/a008/coding/c++/MiniR/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/a008/coding/c++/MiniR /home/a008/coding/c++/MiniR/src /home/a008/coding/c++/MiniR/build /home/a008/coding/c++/MiniR/build/src /home/a008/coding/c++/MiniR/build/src/CMakeFiles/stat.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/stat.dir/depend

