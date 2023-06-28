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
set(sourcesLib "../utils/cutils.cpp" "../utils/data.cpp" "dh2operator.cpp" "compress.cpp" "functions.cpp" "meshreader.cpp" "multisurface.cpp" "options.cpp" "orient.cpp" "surface.cpp" "timing.cpp" "constants.cpp")
set(executableSources "main.cpp")

set(linklibs "h2" "lapack" "blas" "z" "png" "m" "GLU" "GL" "gfortran" "gmshNew")
set(properties " ")
set(targets "logger")

set(executableName "glacier")
set(libraryName "ice")

set(libraryNameStatic "${libraryName}Static")
set(libraryNameObject "${libraryName}Object")
set(targetsExe ${targets})
list(APPEND targetsExe ${libraryNameStatic})

defBI(sourcesLib "OBJECT" linklibs " " properties OFF ${libraryNameObject})
target_compile_options(${libraryNameObject} PUBLIC "-fPIC" "-Dcmake_LIBNAME=\"${libraryNameObject}\"")
set(actualLibSources $<TARGET_OBJECTS:${libraryNameObject}>)
message("objects" ${actualLibSources})
# defBI(${actualLibSources} "STATIC" linklibs targets properties OFF ${libraryNameStatic})
# defBI(${actualLibSources} "SHARED" linklibs targets properties OFF ${libraryName})

message("adding static library \"${progname}\"")
add_library(${libraryNameStatic} STATIC $<TARGET_OBJECTS:${libraryNameObject}>)
target_link_libraries(${libraryNameStatic} ${targets})
message("setting compile def for the objects: -Dcmake_LIBNAME=${libraryName}" )
set(type_v "LIBRARY")

defBI(executableSources "EXECUTABLE" linklibs targetsExe properties OFF ${executableName})


set(testSources "test.cpp")
# defBItestSources "EXECUTABLE" linklibs targetsExe properties OFF "test")

