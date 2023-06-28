# distutils: language=c++
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

cimport options
cimport logger
cdef class Config:
    cdef options.Options *data
    def __cinit__(self):
        logger.addLogStdout(logger.min)
        self.data = options.makeOptions(True)

    def __dealloc__(self):
        options.delOptions()

    def checkInit(self):
        self.data.fillDefaults()
    @property
    def interactive(self):
        return self.data.interactive_
    @interactive.setter
    def interactive(self, interactive):
        self.data.interactive_ = interactive

    @property
    def quadraturePoints(self):
        return self.data.quadraturePoints_
    @quadraturePoints.setter
    def quadraturePoints(self, quadraturePoints):
        self.data.quadraturePoints_ = quadraturePoints
    @property
    def interpolationPoints(self):
        return self.data.interpolationPoints_
    @interpolationPoints.setter
    def interpolationPoints(self, interpolationPoints):
        self.data.interpolationPoints_ = interpolationPoints

    @property
    def clusterResolution(self):
        return self.data.clusterResolution_
    @clusterResolution.setter
    def clusterResolution(self, clusterResolution):
        self.data.clusterResolution_ = clusterResolution

    def compression(self, compression):
        self.data.compression(compression)
    def reference(self, reference):
        self.data.reference(reference)
    def nearfield(self, nearfield):
        self.data.nearfield(nearfield)
    def directionalAdmissibilityData(self, eta1, eta2, eta3, wavek, waveki):
        self.data.directionalAdmissibilityData(eta1, eta2, eta3, wavek, waveki)

    @property
    def eps(self):
        return self.data.eps_
    @eps.setter
    def eps(self, eps):
        self.data.eps_ = eps
    @property
    def epsKrylov(self):
        return self.data.epsKrylov_
    @epsKrylov.setter
    def epsKrylov(self, epsKrylov):
        self.data.epsKrylov_ = epsKrylov

    # @property
    # def directionalAdmissibilityData(self):
    #     return self.data.directionalAdmissibilityData_
    # @directionalAdmissibilityData.setter
    # def directionalAdmissibilityData(self, directionalAdmissibilityData):
    #     self.data.directionalAdmissibilityData_ = directionalAdmissibilityData
        
    @property
    def meshFile(self):
        return self.data.meshFile_
    @meshFile.setter
    def meshFile(self, meshFile):
        self.data.meshFile_ = meshFile

    def report(self):
        self.data.report()
    
    
    # @property
    # def stopWatch_(self):
    #     return self.data.stopWatch_
    # @stopWatch_.setter
    # def stopWatch_(self, stopWatch_):
    #     self.data.stopWatch_ = stopWatch_
    # @property
    # def progressBar_(self):
    #     return self.data.progressBar_
    # @progressBar_.setter
    # def progressBar_(self, progressBar_):
    #     self.data.progressBar_ = progressBar_
    # @property
    # def timer_(self):
    #     return self.data.timer_
    # @timer_.setter
    # def timer_(self, timer_):
    #     self.data.timer_ = timer_
