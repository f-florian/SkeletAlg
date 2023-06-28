###########################################################################
# Copyright (C) 2023 Francesco Florian
# This file is part of SkeletAlg.
#
# SkeletAlg is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SkeletAlg is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SkeletAlg.  If not, see <http://www.gnu.org/licenses/>.
#
# The copyright holders give you permission to combine SkeletAlg
# with code included in the standard release of Netgen (from Joachim
# Schöberl), METIS (from George Karypis at the University of
# Minnesota), OpenCASCADE (from Open CASCADE S.A.S) and ParaView
# (from Kitware, Inc.) under their respective licenses. You may copy
# and distribute such a system following the terms of the GNU GPL for
# Gmsh and the licenses of the other code concerned, provided that
# you include the source code of that other code when and as the GNU
# GPL requires distribution of source code.
#
# Note that people who make modified versions of SkeletAlg are not
# obligated to grant this special exception for their modified
# versions; it is their choice whether to do so. The GNU General
# Public License gives permission to release a modified version
# without this exception; this exception also makes it possible to
# release a modified version which carries forward this exception.
#
# Additional permission under GNU GPL version 3 section 7
# If you modify this Program, or any covered work, by linking or combining it
# with H2Lib (https://github.com/H2Lib/H2Lib), (or a modified version of that
# library), containing parts covered by the terms all right reserved,the
# licensors of this Program grant you additional permission to convey the
# resulting work.
###########################################################################
cmake_minimum_required(VERSION 3.9)
project(DH2Test VERSION 0.0 LANGUAGES CXX)
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0057 NEW)
# see makefile output
set(CMAKE_VERBOSE_MAKEFILE OFF)
# switches
set(debug OFF)
set(openmp ON)
set(logging OFF)
set(timing ON)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)

# set compile option
set(compileOptions -pipe -Ofast -Wall -Wextra -march=native -mavx -mavx2 -mfma -std=gnu++2a -fno-exceptions)

set(utilities "logger" "main")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message("compiler: ${CMAKE_CXX_COMPILER_ID}")
# set complete compile flags depending on the compiler
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  list(APPEND compileOptions -Wdeprecated -Weverything -Wno-zero-as-null-pointer-constant -Wno-vla-extension -Wno-vla -Wno-covered-switch-default -Wno-padded -Wno-float-equal -Wno-old-style-cast -Wno-documentation -Wno-double-promotion -Wno-documentation-unknown-command -Wno-newline-eof -Wno-weak-vtables -Wno-format-nonliteral -Wno-implicit-int-float-conversion -Wno-c99-extensions -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-c++20-compat -Wno-shorten-64-to-32)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  list(APPEND compileOptions -Wall -W -Wextra -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common)
else()
  list(APPEND compileOptions -Wall -W -Wextra -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common)
  message("Warning: unknown compiler, assuming GNU-compatible options!")
endif()

# pass version and name to the program as preprocessor definitions
set(compileDefinitions -Dcmake_MAJOR=${PROJECT_VERSION_MAJOR} -Dcmake_MINOR=${PROJECT_VERSION_MINOR} -Dcmake_NAME="${PROJECT_NAME}" "-DFIELD_COMPLEX" "-DUSE_ZLIB" "-DUSE_AVX" "-DUSE_AVX2" "-DUSE_FMA" "-DUSE_BLAS")
# process switches
if(debug)
  list(APPEND compileDefinitions "-D____debug")
  list(APPEND compileOptions -g -gdwarf-4 -gstrict-dwarf)
else()
  list(APPEND compileOptions -DNDEBUG)  
endif(debug)
if(${logging})
  list(APPEND compileDefinitions "-D____Logging")
endif()
if(${timing})
  list(APPEND compileDefinitions "-D____Timing")
endif()
if(${openmp})
  list(APPEND compileDefinitions "-DUSE_OPENMP")
  list(APPEND compileDefinitions "-DUSE_OPENMP_TASKS")
  list(APPEND compileDefinitions "-DUSE_OPENMP_SIMD")
endif()
  list(APPEND compileOptions -fopenmp)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

add_definitions(${compileDefinitions})
add_compile_options(${compileOptions})

# include module
include("CMakeSystem/buildinstall.cmake")

# add subdirectories to include directories (variable incdir is used in function defBI)
include_directories(SYSTEM "Lyra/include/")
set(incdir ${CMAKE_CURRENT_SOURCE_DIR})
foreach(loopvar IN LISTS utilities)
  list(APPEND incdir ${CMAKE_CURRENT_SOURCE_DIR}/${loopvar})
endforeach(loopvar)

# add subbbdirectories
foreach(loopvar IN LISTS utilities)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${loopvar})
endforeach(loopvar)

# Documentation
message(CHECK_START "Searching for doxygen")
find_package(Doxygen
  OPTIONAL_COMPONENTS dot)
if(${DOXYGEN_FOUND})
  message(CHECK_PASS "found version ${DOXYGEN_VERSION}; optional components found: <todo>")
  set(DOXYGEN_@INCLUDE "Doxyfile")
  doxygen_add_docs(doc
    ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM )
else()
  message(CHECK_FAIL "not found, skipping the target for documentation generation")
endif()

# install license
# install(FILES "LICENSE" DESTINATION "$ENV{PREFIX}/var/lib/${project}" PERMISSIONS WORLD_READ)
