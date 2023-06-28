#!env python3
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
import numpy
import config

# if __name__ == "__main__":
opts = config.Config()
opts.interactive = False
opts.quadraturePoints = 4
opts.interpolationPoints  = 6
opts.clusterResolution  = 32
opts.compression(b'implicit')
opts.reference(b'none')
opts.nearfield(b'recompute')
opts.eps  = 5e-7
opts.epsKrylov  = 2e-9
# opts.directionalAdmissibilityData;
# opts.*stopWatch = nullptr;
# opts.*progressBar = nullptr;
opts.meshFile = b"~/Documents/dev/h2test/meshfiles/easycube.msh"

opts.checkInit()
opts.report()
h2lib_wavek=4.0
h2lib_waveki=4.0
h2lib_eta1=10.0
h2lib_eta2=1.0
h2lib_eta3=0.5

print("done")
print("Test 'planewave z ntd'\n");
# Functions::SinusSolution referenceWave({0,0,1}, waveNumber());
# Functions::DtN offset(&referenceWave);
# // Functions::NtD offset(&referenceWave);
# problem.test(&offset);
# // Logger::log(Logger::Color::cyan, Logger::Level::programFlow, "Test 'planewave z dtn'\n");
# // Functions::SinusSolution referenceWave({0,0,1}, waveNumber());
# // Functions::NtD offset(&referenceWave);
# // problem.test(&offset);
