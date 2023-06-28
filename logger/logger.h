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
#ifndef LOGGER_H
#define LOGGER_H

#include <type_traits>
#include <memory>

#include "utils/enum.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

namespace MainNamespace {
    /**
     * Helper functions to display log messages based on the severity.
     *
     * Unsetting the preprocessor variable _Logging makes all the logging functions inline empty functions, which should then be optimized away at compile time (unless inlining is disabled).
     */
    namespace Logger {
        enum class Color {  //! List of colors accepted by log functions.
            red,              //!< ANSI_COLOR_RED
            green,            //!< ANSI_COLOR_GREEN
            yellow,           //!< ANSI_COLOR_YELLOW
            blue,             //!< ANSI_COLOR_BLUE
            magenta,          //!< ANSI_COLOR_MAGENTA
            cyan,             //!< ANSI_COLOR_CYAN
            reset             //!< ANSI_COLOR_RESET 
        };
        constexpr char const * colorStrings[getUnderlying(Color::reset)+1] = {"\x1b[31m", "\x1b[32m", "\x1b[33m", "\x1b[34m", "\x1b[35m", "\x1b[36m", "\x1b[0m"};
        inline char const * colorString(Color color) {return colorStrings[getUnderlying(color)];}

        namespace {
            /**
             * Auxiliary enum used to count its entries; see Level.
             */
            enum class LevelPrivate : int {
#include "logger.levels"
            };
        }
        /**
         * Possible log levels; negative values are errors, positive are informations, 0 means "output".
         */
        enum class Level : int {
            min = -1 - getUnderlying(LevelPrivate::output),                 //!< don't use it as a message severity level; use it as the log level to see everything
#include "logger.levels"
        };
        static_assert(getUnderlying(Level::output) == 0, "The intended value for \"output\" should be 0!");
        class OutputDevice;
#ifdef ____Logging
        void addLogDevice(std::shared_ptr<OutputDevice> device,     //!< New device to add.
                          Level const logLevel                      //!< Logging level of the device.
                          );                                        //!< \brief Add a logging device.
        void addLogStdout(Level const logLevel      //!< Logging level of the device.
                          );                        //!< \brief Add a logger for stdout.
        void addLogStderr(Level const logLevel      //!< Logging level of the device.
                          );                        //!< \brief Add a logger for stderr.
        void addLogFile(char const * const filename,//!< Name of the file to open.
                        Level const logLevel        //!< Logging level of the device.
                        );                          //!< \brief Add a logger for a given file. The logs are appended to it.
        void indent();      //!< \brief Increase indentation level in the logger.
        void deindent();    //!< \brief Decrease indentation level in the logger.

        template <class... Args>
        inline void log(Color const color,              //!< Color of the log, if the device supports it.
                        Level const severity,           //!< Log severity.
                        bool const newline,             //!< Whether to add a new line. Intended as a mean to counteract the spurious indentation spaces when continuing on the same line. Might have more advantages as well.
                        char const * const format,      //!< Format, as in printf
                        Args... args                    //!< Additional arguments, according to the rules of printf.
                        ) {                             //! Log a message.
            forwardLog(color, severity, newline, format, args...);
        }
        template <class... Args> inline void log(Color const color, Level const severity, char const * const format, Args... args) {forwardLog(color, severity, false, format, args...);}       //!< Equivalent to log(color, severity, format, false, args)
        template <class... Args> inline void log(Level const severity, char const * const format, Args... args) {forwardLog(Color::reset, severity, false, format, args...);}                   //!< Equivalent to log(Color:reset, severity, format, args)
        template <class... Args> inline void log(char const * const format, Args... args) {forwardLog(Color::reset, Level::warning, false, format, args...);}                                   //!< Equivalent to log(Level::Warning, format, args)

        void forwardLog(Color const color,              //!< As in OutputDevice::print.
                        Level const severity,           //!< As in OutputDevice::print.
                        bool const newline,             //!< Whether to add a new line. Intended as a mean to counteract the spurious indentation spaces when continuing on the same line. Might have more advantages as well.
                        char const *const format,       //!< As in printf.
                        ...                             //!< As in printf.
                        );                              //! Build the message to print and forward it to the loggers.
        constexpr bool logEnabled = true;   //!< Compile-time constant indicating whether logging support was enabled.
#else
        constexpr bool logEnabled = false;
        inline void addLogDevice([[maybe_unused]] std::shared_ptr<OutputDevice>, [[maybe_unused]] Level const) {}
        constexpr inline void addLogStdout([[maybe_unused]] Level const) {}
        constexpr inline void addLogStderr([[maybe_unused]] Level const) {}
        constexpr inline void addLogFile([[maybe_unused]]char const * const, [[maybe_unused]] Level const) {}
        constexpr inline void indent() {}
        constexpr inline void deindent() {}

        constexpr inline void log([[maybe_unused]] Color const, [[maybe_unused]] Level const, [[maybe_unused]] bool const, [[maybe_unused]] char const * const, ...) {}
        constexpr inline void log([[maybe_unused]] Color const, [[maybe_unused]] Level const, [[maybe_unused]] char const * const, ...) {}
        constexpr inline void log([[maybe_unused]] Level const, [[maybe_unused]] char const * const, ...) {}
        constexpr inline void log([[maybe_unused]] char const * const, ...) {}
        constexpr inline void forwardLog([[maybe_unused]] Color const, [[maybe_unused]] Level const, [[maybe_unused]] bool, [[maybe_unused]] char const *const, ...) {}
#endif
    }
}

#endif /* LOGGER_H */
