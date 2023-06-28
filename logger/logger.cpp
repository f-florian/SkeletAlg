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
#ifdef ____Logging
#include <cstdio>
#include <cstdarg>
#include <list>
#include <memory>
#include <cassert>

#include "devices.h"
#include "logger.h"

namespace MainNamespace {
    namespace Logger {
        namespace {
            std::list<std::shared_ptr<OutputDevice>> logs_;
        }
        void forwardLog(Color const color, Level const severity, bool const newline, char const *const format, ...) {
            va_list args;
            va_start(args, format);
            va_list argsCopy;
            va_copy(argsCopy, args);
            int bytes = vsnprintf(nullptr, 0, format, argsCopy);
            va_end(argsCopy);
            ++bytes;        // The null character at the end is not in the count.
            char *buffer = static_cast<char*>(malloc(bytes));
            (void) vsnprintf(buffer, bytes, format, args);
            va_end(args);
  
            for(auto &logger: logs_)
                logger->print(severity, color, buffer, newline);
            free(buffer);
        }
        void addLogDevice(std::shared_ptr<OutputDevice> device, Level const logLevel) {
            logs_.push_back(device);
            if(logs_.back()->logLevel(logLevel) != logLevel)
                fprintf(stderr, "Failed to set logLevel on new device");
            log(Color::red, Level::developer, "New log device\n");
        }
        void addLogStdout(Level const logLevel) {
            addLogDevice(std::make_shared<StdStreamDevice>(stdout), logLevel);
        }
        void addLogStderr(Level const logLevel) {
            addLogDevice(std::make_shared<StdStreamDevice>(stdout), logLevel);
        }
        void addLogFile(char const * const filename, Level const logLevel) {
            addLogDevice(std::make_shared<FilePointerDevice>(filename, "a"), logLevel);
        }
        void indent() {
            for(auto &logger: logs_)
                logger->indent();
        }

        void deindent() {
            for(auto &logger: logs_)
                logger->deindent();
        }
    }
}
#endif
