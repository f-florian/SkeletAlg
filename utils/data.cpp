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
#include "data.h"

#include <cstdio>

#include "h2lib/basic.h"

#include "logger/logger.h"

namespace MainNamespace {
    void reportParallelization() {
#ifdef USE_OPENMP
        // If OpenMP is enabled, report maximal parallization depth
        Logger::log(Logger::Color::green, Logger::Level::ioBasic, true, "Parallelization: h2lib max depth: %u; openmp max threads: %u; openmp threads limit: %u", max_pardepth, omp_get_max_threads(), omp_get_thread_limit());
#else
        Logger::log(Logger::Color::red, Logger::Level::warning, true, "Parallelization disabled.");
#endif
    }
    FILE* outputFile(std::string const name, bool const replace) {
        Logger::log(Logger::Color::green, Logger::Level::warning, false, "Opening %s in write%s mode...", name.c_str(), replace?"":" (append)");
        if(access(name.c_str(), W_OK)) {
            Logger::log(Logger::Color::red, Logger::Level::warning, true, "Failed!");
            return nullptr;
        }
        Logger::log(Logger::Color::green, Logger::Level::debug, true, "OK");
        return fopen(name.c_str(), replace?"w":"a");
    }
}
