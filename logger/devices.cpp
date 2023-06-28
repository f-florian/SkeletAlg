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
#include <cstdio>

#include "devices.h"

static constexpr size_t indentationSpaces = 2;

namespace MainNamespace {
    namespace Logger {
        Level OutputDevice::logLevel(Level const severity) {
            logLevel_ = severity;
            return logLevel_;
        }
        Level OutputDevice::logLevel() const {
            return logLevel_;
        }
        bool OutputDevice::checkMessageLevel(Level const level) const {
            return (level >= logLevel_);
        }

        StdStreamDevice::StdStreamDevice(std::FILE *stream)
            : stream_(stream) {
            blanks_ = new char[1];
            blanks_[0] = '\0';
        }
        StdStreamDevice::~StdStreamDevice() {
            fflush(stream_);
            delete[] blanks_;
        }
        void StdStreamDevice::print(Level const severity, Color const color, char const * const message, bool const newline) {
            if(checkMessageLevel(severity)) {
                if(newline_)
                    fprintf(stream_, "%s%s%s%s", colorString(color), blanks_, message, colorString(Color::reset));
                else
                    fprintf(stream_, "%s%s%s", colorString(color), message, colorString(Color::reset));
                if(newline)
                    fprintf(stream_, "\n");
                newline_ = newline;
            }
#ifdef ____debug
            fflush(stream_);
#endif
        }
        void StdStreamDevice::indent() {
            indentation_ += indentationSpaces;
            delete[] blanks_;
            blanks_ = new char[indentation_ * indentationSpaces + 1];
            for(size_t idx = 0; idx < indentation_ * indentationSpaces; ++idx)
                blanks_[idx] = ' ';
            blanks_[indentation_ * indentationSpaces] = '\0';
        }
        void StdStreamDevice::deindent() {
            indentation_ -= indentationSpaces;
            blanks_[indentation_] = '\0';
        }

        FilePointerDevice::FilePointerDevice(char const * const __restrict__ filename, char const * const __restrict__ modes) {
            file_ = fopen(filename, modes);
            blanks_ = new char[1];
            blanks_[0] = '\0';
        }
        FilePointerDevice::~FilePointerDevice() {
            fflush(file_);
            fclose(file_);
            delete[] blanks_;
        }
        void FilePointerDevice::print(Level const severity, [[maybe_unused]] Color const color, char const * const message, bool const newline) {
            if(checkMessageLevel(severity)) {
                if(newline_)
                    fprintf(file_, "%s%s", blanks_, message);
                else
                    fprintf(file_, "%s", message);
                if(newline)
                    fprintf(file_, "\n");
                newline_ = newline;
            }
        }
        void FilePointerDevice::indent() {
            indentation_ += indentationSpaces;
            delete[] blanks_;
            blanks_ = new char[indentation_ * indentationSpaces + 1];
            for(size_t idx = 0; idx < indentation_ * indentationSpaces; ++idx)
                blanks_[idx] = ' ';
            blanks_[indentation_ * indentationSpaces] = '\0';
        }
        void FilePointerDevice::deindent() {
            indentation_ -= indentationSpaces;
            blanks_[indentation_] = '\0';
        }
    }
}
