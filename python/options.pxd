###########################################################################
# Copyright (C) 2023 Francesco Florian
# This file is part of SkeletAlg.

# SkeletAlg is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# SkeletAlg is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with SkeletAlg.  If not, see <http://www.gnu.org/licenses/>.

# The copyright holders give you permission to combine SkeletAlg
# with code included in the standard release of Netgen (from Joachim
# Schöberl), METIS (from George Karypis at the University of
# Minnesota), OpenCASCADE (from Open CASCADE S.A.S) and ParaView
# (from Kitware, Inc.) under their respective licenses. You may copy
# and distribute such a system following the terms of the GNU GPL for
# Gmsh and the licenses of the other code concerned, provided that
# you include the source code of that other code when and as the GNU
# GPL requires distribution of source code.

# Note that people who make modified versions of SkeletAlg are not
# obligated to grant this special exception for their modified
# versions; it is their choice whether to do so. The GNU General
# Public License gives permission to release a modified version
# without this exception; this exception also makes it possible to
# release a modified version which carries forward this exception.

# Additional permission under GNU GPL version 3 section 7
# If you modify this Program, or any covered work, by linking or combining it
# with H2Lib (https://github.com/H2Lib/H2Lib), (or a modified version of that
# library), containing parts covered by the terms all right reserved,the
# licensors of this Program grant you additional permission to convey the
# resulting work.
###########################################################################
cdef extern from "dblock.h":
    ctypedef struct diradmdata:
        double eta1
        double eta2
        double eta3
        
from libcpp.string cimport string

cdef extern from 'main/options.h' namespace 'MainNamespace':
    Options* makeOptions(bool)
    Options* delOptions()
    cdef cppclass Options:
        void fillDefaults()
        bint interactive_
        size_t quadraturePoints_
        size_t interpolationPoints_
        size_t clusterResolution_
        double eps_
        double epsKrylov_
        string meshFile_
        void directionalAdmissibilityData(double, double, double, double, double)
        void compression(string)
        void reference(string)
        void nearfield(string)
        void report()
