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
#ifndef DEVICES_H
#define DEVICES_H

#include <cstdio>

#include "logger.h"

namespace MainNamespace {
    namespace Logger {
        class OutputDevice {    //! \brief An abstraction for classes which write logs.
        public:
            virtual ~OutputDevice() {}
            virtual void print(Level const severity,            //!< Level of the message
                               Color const color,               //!< Color to use for printing the message (if supported).
                               char const * const message,      //!< Message to write.
                               bool const newline               //!< If true includes a newline at the end and indents next message appropriately.
                               ) = 0;                           //!< \brief Print a message to the device.
            virtual void indent() = 0;  //!< Increase the current indentation level.
            virtual void deindent() = 0;//!< Decrease the current indentation level.
            virtual Level logLevel(Level const severity //!< \brief New severity level.
                                   );                   //!< \brief Set the log level. \return The newly set level. It may differ from the input if setting it is not successful (e.g., if logging is disabled).
            virtual Level logLevel() const;             //!< \brief Getter. \return the current log level.
        protected:
            virtual bool checkMessageLevel(Level const level    //!< Message level.
                                           ) const;             //!< \brief Check whether message should be logged. \return true for messages that should be logged.
            Level logLevel_;    //!< Minimum level of the messages to log.
        };

        /**
         * \brief Logger for a standard stream.
         * Log to a given FILE pointer, writing log messages immediately.
         * Color and lack of seek are assumed.
         * The stream is flushed but the file is not closed at destruction time.
         * No check that the file is still open at any point is performed.
         * It is recommended to use this only for stdout or stderr.
         */
        class StdStreamDevice
            : public OutputDevice {
        public:
            StdStreamDevice(std::FILE *stream   //!< Stream to write to.
                            );                  //!< \brief Constructor
            virtual ~StdStreamDevice() override;
            virtual void print(Level const severity,
                               Color const color,
                               char const * const message,
                               bool const newline
                               ) override;
            virtual void indent() override;
            virtual void deindent() override;
        private:
            bool newline_ = true;       //!< Whether a newline was used and next message should be indented.
            std::FILE *stream_;         //!< Stream used for writing.
            size_t indentation_ = 0;    //!< Current intendation width.
            char *blanks_;              //!< String of spaces to match current intentation, to avoid a for loop for every write.
        };

        /**
         * \brief Logger for a file stored on the filesystem.
         * Log to a given file.
         * Lack of color support is assumed.
         * The file is open during construction and closed at destruction time, after flushing the stream.
         * No check that the file open operation was succesul is performed.
         */
        class FilePointerDevice
            : public OutputDevice {
        public:
            FilePointerDevice(char const * const __restrict__ filename, //!< Filename to open
                              char const * const __restrict__ modes     //!< Opening mode. Only w(+) and a(+) will work, although this is not checked.
                              );                                        //!< \brief Constructor
            virtual ~FilePointerDevice() override;
            virtual void print(Level const severity,
                               Color const color,
                               char const * const message,
                               bool const newline
                               ) override;
            virtual void indent() override;
            virtual void deindent() override;
        private:
            bool newline_ = true;       //!< Whether a newline was used and next message should be indented.
            std::FILE *file_;           //!< File object for writing.
            size_t indentation_ = 0;    //!< Current intendation width.
            char *blanks_;              //!< String of spaces to match current intentation, to avoid a for loop for every write.
        };
    }
}

#endif /* DEVICES_H */
