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
#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>

#include "h2lib/parameters.h"
#include "h2lib/reporter.h"

#include "logger/logger.h"
#include "utils/cutils.h"
#include "timing.h"

namespace MainNamespace {
    struct AdmissibilityParameters {    //!< Admissibility parameters (see H2Lib).
        double eta1;            //!< Angular admissibility parameter.
        double eta2;            //!< Parabolic admissibility parameter.
        double eta3;            //!< Admissibility parameter for real part.
        field waveNumber;       //!< Wave number.
    };
    /**
     * Wrapper for Lyra, with the correct variables for this project
     *
     * See the parser options for a description of the variables
     */
    struct Options {
        Options(int argc, char** argv);
        Options(bool const deferInit = false);
        ~Options();

        void fillDefaults(); //!< Called by the constructor by default.

        // variables
        bool interactive_ = false;
        bool interactive() const { return interactive_; }
        bool interactive(bool interactive_p) {
            interactive_ = interactive_p;
            return interactive();
        }
        size_t quadraturePoints_ = 0;
        size_t quadraturePoints() const { return quadraturePoints_; }
        size_t quadraturePoints(size_t quadraturePoints_p) {
            quadraturePoints_ = quadraturePoints_p;
            return quadraturePoints();
        }
        size_t interpolationPoints_ = 0;
        size_t interpolationPoints() const { return interpolationPoints_; }
        size_t interpolationPoints(size_t interpolationPoints_p) {
            interpolationPoints_ = interpolationPoints_p;
            return interpolationPoints();
        }
        /** Cluster resolution */
        size_t clusterResolution_ = 0;
        size_t clusterResolution() const { return clusterResolution_; }
        size_t clusterResolution(size_t clusterResolution_p) {
            clusterResolution_ = clusterResolution_p;
            return clusterResolution();
        }
        /** Compression strategy */
        bool compression_ = true;
        bool compression() const { return compression_; }
        bool compression(bool const);
        bool reference_ = false;        //!< Reference matrix
        bool reference() const { return reference_; }
        bool reference(bool const reference_p) {
            reference_ = reference_p;
            return reference();
        }
        // /** Copy or recompute nearfield */
        // NearField nearfield_ = NearField::invalid;
        // NearField nearfield() const { return nearfield_; }
        // NearField nearfield(std::string);
        real eps_ = 0;
        real eps() const { return eps_; }
        real eps(real eps_p) {
            eps_ = eps_p;
            return eps();
        }
        real epsKrylov_ = 0;
        real epsKrylov() const { return epsKrylov_; }
        real epsKrylov(real epsKrylov_p) {
            epsKrylov_ = epsKrylov_p;
            return epsKrylov();
        }
        AdmissibilityParameters admissibilityParameters_ {.eta1 = 0, .eta2 = 0, .eta3 = 0, .waveNumber = 0};
        AdmissibilityParameters admissibilityParameters() const {return admissibilityParameters_;}
        AdmissibilityParameters admissibilityParameters(AdmissibilityParameters const admissibilityParameters_p) {
            admissibilityParameters_ = admissibilityParameters_p;
            return admissibilityParameters();
        }
        AdmissibilityParameters admissibilityParameters(double const eta1, double const eta2, double const eta3, double const wavek, double const waveki);
        field waveNumber() const {return admissibilityParameters().waveNumber;}
        bool saveMesh_ = false;
        bool saveMesh() const { return saveMesh_; }
        bool saveMesh(bool saveMesh_p) {
            saveMesh_ = saveMesh_p;
            return saveMesh();
        }
        bool saveSolution_ = false;
        bool saveSolution() const { return saveSolution_; }
        bool saveSolution(bool saveSolution_p) {
            saveSolution_ = saveSolution_p;
            return saveSolution();
        }
        bool saveOffset_ = false;
        bool saveOffset() const { return saveOffset_; }
        bool saveOffset(bool saveOffset_p) {
            saveOffset_ = saveOffset_p;
            return saveOffset();
        }
        bool saveRhs_ = false;
        bool saveRhs() const { return saveRhs_; }
        bool saveRhs(bool saveRhs_p) {
            saveRhs_ = saveRhs_p;
            return saveRhs();
        }
        bool saveError_ = false;
        bool saveError() const { return saveError_; }
        bool saveError(bool saveError_p) {
            saveError_ = saveError_p;
            return saveError();
        }
        
        reporter *progressBar_ = nullptr;
        std::string meshFile_ = "";
        std::string meshFile() const {return meshFile_;}
        std::string outputDirectory_ = "";
        std::string outputDirectory() const {return outputDirectory_;}
        Clocks timer_;

        void report() const;
    };

    Options* globalOptions();
    Options* globalOptions(Options* newOpts);
    Options* makeOptions(bool const deferInit);
    Options* delOptions();
    
    inline bool interactive() {return globalOptions()->interactive();}
    inline size_t quadraturePoints() {return globalOptions()->quadraturePoints();}
    inline size_t interpolationPoints() {return globalOptions()->interpolationPoints();}
    inline size_t clusterResolution() {return globalOptions()->clusterResolution_;}
    inline bool compression() {return globalOptions()->compression();}
    inline bool reference() {return globalOptions()->reference();}
    // inline Options::NearField nearfield() {return globalOptions()->nearfield();}
    inline real eps() {return globalOptions()->eps();}
    inline real epsKrylov() {return globalOptions()->epsKrylov();}
    inline reporter * progressBar() {return globalOptions()->progressBar_;}
    inline std::string meshFile() {return globalOptions()->meshFile();}
    inline std::string outputDirectory() {return globalOptions()->outputDirectory();}

    constexpr static size_t dimKrylov_K = 20;
    inline size_t dimKrylov() { return dimKrylov_K; }
#ifdef ____Timing
    inline auto startTiming() {return globalOptions()->timer_.start();}
    inline std::string stopTiming() {return globalOptions()->timer_.readAll();}
    constexpr bool timingEnabled = true;
    inline stopwatch* createStopWatch() {return new_stopwatch();}
#else
    constexpr inline auto startTiming() {}
    inline std::string const stopTiming() {return "";}
    constexpr bool timingEnabled = false;
    constexpr inline stopwatch* createStopWatch() {return nullptr;}
#endif
    namespace ReporterHelper {
        inline void fakeUninit(preporter rp [[maybe_unused]]){}
        inline void fakeBegin(uint total [[maybe_unused]], preporter rp [[maybe_unused]]) {}
        inline void fakeStep(preporter rp [[maybe_unused]]) {}
        inline void fakeEnd(preporter rp [[maybe_unused]]) {}
    }
    inline reporter* createFakeProgressBar() {
        reporter * aux = new_reporter();
        aux->uninit = ReporterHelper::fakeUninit;
        aux->begin = ReporterHelper::fakeBegin;
        aux->step = ReporterHelper::fakeStep;
        aux->end = ReporterHelper::fakeEnd;
        return aux;
    }
    inline reporter* createProgressBar() { return NULL; }
    // inline reporter* createProgressBar() { return new_textbar_reporter(); }

    inline bool saveMesh() { return globalOptions()->saveMesh(); }
    inline bool saveSolution() { return globalOptions()->saveSolution(); }
    inline bool saveOffset() { return globalOptions()->saveOffset(); }
    inline bool saveRhs() { return globalOptions()->saveRhs(); }
    inline bool saveError() { return globalOptions()->saveError(); }
}

#endif /* OPTIONS_H */
