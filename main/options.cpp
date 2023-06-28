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
#include "options.h"

#include <cstdlib>
#include <sstream>
#include <cstring>

#include "lyra/lyra.hpp"

#include "logger/logger.h"
#include "utils/data.h"


namespace MainNamespace {
    namespace {
        Options * globalOptions_ = nullptr;
    }
    Options* globalOptions() {return globalOptions_;}
    Options* globalOptions(Options* newOpts) {globalOptions_ = newOpts; return globalOptions_;}
    Options* makeOptions(bool const deferInit) {
        new Options(deferInit); // no need to store it by hand: address automatically stored it the global pointer
        return globalOptions();
    }
    Options* delOptions() {
        delete globalOptions();
        return globalOptions();
    }

    void Options::fillDefaults() {
        if(! interactive_)
            return;

        if(quadraturePoints()  == 0)
            quadraturePoints(static_cast<size_t>(askforint("Quadrature order?", "h2lib_quadorder", 3)));
        if(interpolationPoints() == 0)
            interpolationPoints(static_cast<size_t>(askforint("Interpolation order?", "h2lib_interorder", 4)));
        if(clusterResolution() == 0)
            clusterResolution(static_cast<size_t>(askforint("Cluster resolution?", "h2lib_leafsize", 32)));
        if(eps() == 0)
            eps(askforreal("Compression tolerance?", "h2lib_compeps", 1e-5));
        if(epsKrylov() == 0)
            epsKrylov(askforreal("Krylov solver tolerance?", "h2lib_solvetol", 1e-9));
        if(admissibilityParameters_.eta1 == 0)
            admissibilityParameters_.eta1 = askforreal("Directional admissibility parameter eta1?", "h2lib_eta1", 10.);
        if(admissibilityParameters_.eta2 == 0)
            admissibilityParameters_.eta2 = askforreal("Parabolic  admissibility parameter eta2?", "h2lib_eta2", 2.);
        if(admissibilityParameters_.eta3 == 0)
            admissibilityParameters_.eta3 = askforreal("Decay admissibility parameter eta3?", "h2lib_eta3", 0.5);

        if(admissibilityParameters_.waveNumber == 0) {
            double kappa_r = askforreal("Wave number?", "h2lib_wavek", 1.);
            double kappa_i = askforreal("Imaginary part of the wave number?", "h2lib_waveki", 0.);
            admissibilityParameters_.waveNumber =  kappa_r + kappa_i * M_I;
        }
    }
    Options::~Options() {
        if(globalOptions() == nullptr) {
            Logger::log(Logger::Color::red, Logger::Level::warning,
                        "The global options pointer was already reset to a null pointer!\n"
                        "Skipping h2lib uninitialization\n");
        } else {
            Logger::log(Logger::Level::programFlow, "Uninit h2lib\n");
            uninit_h2lib();
            if(globalOptions(nullptr) != nullptr)
                Logger::log(Logger::Level::warning, "Resetting the global options pointer failed!\n");
        }

        maybeFree(progressBar_);
    }
    Options::Options(int argc, char **argv) {
        // Initializing global objects
        if(globalOptions() != nullptr) {
            Logger::log(Logger::Color::red, Logger::Level::error,
                        "Attepmting to create an 'Options' object, when one exists already!\n"
                        "The global pointer will retain the previous value, but will be reset during this object's destruction!\n");
        } else {
            init_h2lib(&argc, &argv); // init library
            if(globalOptions(this) != this)
                Logger::log(Logger::Level::warning, "Setting the global options pointer failed!");
            else
                Logger::log(Logger::Level::warning, "Setting the global options pointer succeed!");
        }
    
        // Preparing parse
        //! \todo also get environment here.
        bool help = false, version = false;
        admissibilityParameters_.waveNumber = 0;
        admissibilityParameters_.eta1 = admissibilityParameters_.eta2 = admissibilityParameters_.eta3 = 0;

        lyra::cli parser =
            lyra::opt(interactive_)["-i"]["--interactive"]("Interactively handle errors and missing parameters") |
            lyra::opt(meshFile_, "Name of the file contaning the mesh")["-f"]["--file"] |
            lyra::opt(outputDirectory_, "Name of the output directory")["-o"]["--output"] |
            lyra::opt(quadraturePoints_, "Number of quadrature points (per direction?) (order?)")["-q"]["--quadrature"] |
            lyra::opt(interpolationPoints_, "Number of interpoliation points (order?)")["-I"]["--interpolation"] |
            lyra::opt(clusterResolution_, "Cluster resolution (leaf size)")["-r"]["--resolution"] |
            lyra::opt(reference_)["-R"]["--reference"]("Compute reference and compression error") |
            lyra::opt(eps_, "Compression tolerance")["-t"]["--tolerance"] |
            lyra::opt(epsKrylov_, "Tolerance for the Krylov solver")["--epskrylov"] |
            lyra::opt(admissibilityParameters_.eta1, "Directional admissibility parameter")["--eta1"] |
            lyra::opt(admissibilityParameters_.eta2, "Parabolic admissibility parameter")["--eta2"] |
            lyra::opt(admissibilityParameters_.eta3, "Decay admissibility parameter")["--eta3"] |
            lyra::opt([this](double argument){admissibilityParameters_.waveNumber = cimag(admissibilityParameters_.waveNumber) * I + argument;}, "Wave number (real part)")["-k"]["--wavereal"] |
            lyra::opt([this](double argument){admissibilityParameters_.waveNumber = creal(admissibilityParameters_.waveNumber) + I * argument;}, "Wave number (imaginary part)")["-K"]["--waveimag"] |
            lyra::opt(saveMesh_)["-M"]["--saveMesh"] |
            lyra::opt(saveOffset_)["-O"]["--saveOffset"] |
            lyra::opt(saveSolution_)["-S"]["--saveSolution"] |
            lyra::opt(saveRhs_)["-s"]["--saveRhs"] |
            lyra::opt(saveError_)["-e"]["--saveError"] |
            lyra::opt(version)["-v"]["--version"] |
            lyra::help(help);

        lyra::parse_result ok = parser.parse({argc, argv});
        if(!ok)
            Logger::log(Logger::Color::red, Logger::Level::error, "UNKNOWN ERROR WHEN PARSING COMMAND LINE!");

        if(version) {
            printf("%s version %d.%d\nTODO: license info\n", cmake_NAME, cmake_MAJOR, cmake_MINOR);
        }
        if(help) {
            std::stringstream helpStream(cmake_NAME);
            helpStream << parser;
            printf("%s\n", helpStream.str().c_str());
        }
        if(help || version) {
            exit(0);
        }

        BaseTimer *timer=new CpuTimer;
        timer_.addTimer(timer);
        timer=new BoermTimer;
        timer_.addTimer(timer);
        timer_.format("time: ");
        progressBar_ = createProgressBar(); // Create a reporter for progress bars
    
        fillDefaults();
    }
    Options::Options(bool const deferInit) {
        globalOptions(this);
        int argc = 1;
        char *arg = new char[7];
        char **argv = &arg;
        char const tmp[7] = "No-opt";
        memcpy(arg, tmp, 7 * sizeof(char));
        delete[] arg;
        init_h2lib(&argc, &argv);

        BaseTimer *timer=new CpuTimer;
        timer_.addTimer(timer);
        timer=new BoermTimer;
        timer_.addTimer(timer);
        timer_.format("time: ");
        progressBar_ = createProgressBar(); // Create a reporter for progress bars

        if(deferInit)
            return;
        fillDefaults();
    }

    bool Options::compression(bool const compression_p) {
        compression_ = compression_p;
        return compression();
    }

    AdmissibilityParameters Options::admissibilityParameters(double const eta1, double const eta2, double const eta3, double const wavek, double const waveki) {
        admissibilityParameters_.eta1 = eta1;
        admissibilityParameters_.eta2 = eta2;
        admissibilityParameters_.eta3 = eta3;
        admissibilityParameters_.waveNumber = wavek + _Complex_I * waveki;
        return admissibilityParameters_;
    }

    void Options::report() const {
        Logger::log(Logger::Color::green,
                    Logger::Level::developer,
                    "Current options:\n"
                    "Interactive: %d\n"
                    "Quadrature and interpolation points: %ld; %ld\n"
                    "Resolution: %ld\n"
                    "Compression: %d\n"
                    "Mesh file: %s\n",
                    static_cast<int>(interactive_),
                    quadraturePoints_, interpolationPoints_,
                    clusterResolution_, static_cast<int>(compression_),
                    meshFile_.c_str());
    }
}
