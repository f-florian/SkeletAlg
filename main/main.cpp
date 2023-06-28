/***************************************************************************
 * Copyright (C) 2023 Francesco Florian
 * This file is part of SkeletAlg.
 *
 * SkeletAlg is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SkeletAlg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SkeletAlg.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The copyright holders give you permission to combine SkeletAlg
 * with code included in the standard release of Netgen (from Joachim
 * Schöberl), METIS (from George Karypis at the University of
 * Minnesota), OpenCASCADE (from Open CASCADE S.A.S) and ParaView
 * (from Kitware, Inc.) under their respective licenses. You may copy
 * and distribute such a system following the terms of the GNU GPL for
 * Gmsh and the licenses of the other code concerned, provided that
 * you include the source code of that other code when and as the GNU
 * GPL requires distribution of source code.
 *
 * Note that people who make modified versions of SkeletAlg are not
 * obligated to grant this special exception for their modified
 * versions; it is their choice whether to do so. The GNU General
 * Public License gives permission to release a modified version
 * without this exception; this exception also makes it possible to
 * release a modified version which carries forward this exception.
 *
 * Additional permission under GNU GPL version 3 section 7
 * If you modify this Program, or any covered work, by linking or combining it
 * with H2Lib (https://github.com/H2Lib/H2Lib), (or a modified version of that
 * library), containing parts covered by the terms all right reserved,the
 * licensors of this Program grant you additional permission to convey the
 * resulting work.
 ***************************************************************************/
#include "logger/logger.h"
#include "utils/functionals.h"
#include "utils/data.h"

#include "functions.h"
#include "multisurface.h"
#include "options.h"

using namespace MainNamespace;

SpaceVector<double const> const direction = {1,0,0};

static void testProblem(AdmissibilityParameters const admissibilityParameters,   //!< Wave number s.
                        std::string const meshFile,                              //!< Input mesh file.
                        Functions::MeshOffset &offset,                           //!< Exact solution, for error measurement, and offset values constructed from it.
                        std::string const fileId                                 //!< Filename where the error norm should be written and prefix for the solution file.
                        ) {                                                      //! Compute the solution of the BVP with boundary data given by solution, and record the error on a file.
    FILE *errorFile = outputFile(outputDirectory() + fileId, false);
    if(errorFile == nullptr) abort();

    Logger::log(Logger::Color::cyan, Logger::Level::programFlow, "Constructing geometry...'\n");
    Logger::indent();
    CpuTimer cpuClock;
    BoermTimer wallClock;
    cpuClock.start();
    wallClock.start();
    Surface::MultiSurface geometry(meshFile);       // Setup geometry and metric information
    cpuClock.stop();
    wallClock.stop();
    fprintf(errorFile, "[%6.0lf, %8.0lf,", wallClock.read(), cpuClock.read());
  
    Logger::deindent();
    Logger::log(Logger::Color::cyan, Logger::Level::programFlow, "Solving the problem and computing errors.\n");
    Logger::indent();
    cpuClock.start();
    wallClock.start();
    size_t iterations;
    Surface::ErrorVector error;
    std::tie(iterations, error) = geometry.solve(admissibilityParameters, offset);
    cpuClock.stop();
    wallClock.stop();
    fprintf(errorFile, "%6.0lf,%8.0lf, %4ld, %3.0lf,", wallClock.read(), cpuClock.read(), iterations, cimag(admissibilityParameters.waveNumber));
    Logger::deindent();

    double hmax, ratio; std::tie(hmax,ratio) = geometry.properties();
    fprintf(errorFile, "%7ld, %.3e, %.3e", geometry.problemSize(), hmax, ratio);
    fprintf(errorFile, ", %.5e", error);
    // for(size_t idx = 0; idx < Surface::MultiSurface::errorsCount; ++idx)
    //     fprintf(errorFile, ", %.5e", errors[idx]);
    fprintf(errorFile, "],\n");
    fclose(errorFile);
    Logger::log(Logger::Color::green, Logger::Level::programFlow, "Closed error file.\n");
}

int main(int argc, char **argv) {
    Logger::addLogStdout(Logger::Level::ioBasic);       // Log messages to stdout.
    Options options(argc,argv);                 // Initialize parameters.
    // Logger::log(Logger::Color::magenta, Logger::Level::programFlow, "Done loading options and logger\n");
    Logger::log(Logger::Color::magenta, Logger::Level::programFlow, "Done loading options and logger\n");
    reportParallelization();
    Logger::log(Logger::Color::cyan, Logger::Level::programFlow, "Test 'planewave mixed'. Initializing...\n");
    Functions::PlaneWave referenceWave(direction, options.waveNumber());// Problem definition: exact solution.
    Functions::ZeroExtension offset(&referenceWave);                    // Problem definition: generic zero offset.

    Logger::indent();
    testProblem(options.admissibilityParameters(), options.meshFile(), offset, "mix");
    Logger::deindent();

    Logger::log("Active data (directional cluster bases, matrices, vector): %u, %u, %u\n",
                getactives_dclusterbasis(), getactives_amatrix(), getactives_avector());
    return 0;
}
