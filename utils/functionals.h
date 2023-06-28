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
#ifndef FUNCTIONALS_H
#define FUNCTIONALS_H

#include <type_traits>
#include <functional>

#include "main/options.h"

namespace MainNamespace {
    // template<class T1, class T2> using maybePair = std::conditional<std::is_void_v<T2>, T1, std::pair<T1,T2>>;
    template<class Type>
    using pointerOrReference = std::conditional_t<std::is_pointer_v<Type>, Type,
                                                  std::conditional<std::is_reference_v<Type>, Type,
                                                                   std::add_lvalue_reference_t<Type>>>;
    template<class Type>
    using constPointer = std::conditional_t<std::is_pointer_v<Type>,
                                            std::add_const_t<std::add_pointer_t<std::add_const_t<std::remove_pointer_t<Type>>>>,
                                            std::add_const_t<Type>>;
    // template <class Type>
    // std::enable_if<std::conditional<std::is_pointer_v<Type>, std::true, std::is_reference_v<Type>>> makePointerOrReference(Type arg) {return arg};

    // template <class Type>
    // std::enable_if<std::conditional<std::is_pointer_v<Type>, std::false, std::is_reference_v<Type>>> makePointerOrReference(Type arg) {return arg};

    template<class Return, class... Args>
    typename std::enable_if_t<std::negation_v<std::is_void<Return>>, Return
                              > timeExecution(std::function<Return(Args...)> function,                          //!< Function whose running time should be measured.
                                              std::function<void(std::string const,                             //!< Elapsed time as a string.
                                                                 pointerOrReference<constPointer<Return>>       //!< The value returned by \p function
                                                                 )> logger,                                     //!< Function which is called with the elapsed time and the return value as arguments. The intended use is logging messages.
                                              Args... args                                                      //!< Arguments to pass to \p function.
                                              ) {                                                               //! Execute \p function with \p args... as arguments, log the elapsed time and possibly some other message, according to \p logger. This version is called when \p function returns someting other than void. \return The value returned by \p function.
        startTiming();
        Return returnValue = function(args...);
        logger(stopTiming(), returnValue);
        return returnValue;
    }
    template<class Return, class... Args>
    typename std::enable_if_t<std::is_void_v<Return>, void
                              > timeExecution(std::function<void(Args...)> function,            //!< Function whose running time should be measured.
                                              std::function<void(std::string const)> logger,    //!< Function which is called with the elapsed time. The intended use is logging messages.
                                              Args... args                                      //!< Arguments to pass to \p function.
                                              ) {                                               //! Execute \p function with \p args... as arguments, log the elapsed time and possibly some other message, according to \p logger. This version is called when \p function returns void.
        startTiming();
        function(args...);
        logger(stopTiming());
    }
    template<class Return, class... Args>
    typename std::enable_if_t<std::negation_v<std::is_void<Return>>, Return
                              > timeExecution(Return (*function)(Args...),
                                              std::function<void(std::string const,
                                                                 pointerOrReference<constPointer<Return>>
                                                                 )> logger,
                                              Args... args
                                              ) {
        return timeExecution<Return, Args...>(static_cast<std::function<Return(Args...)>>(function), logger, args...);
  }
  template<class Return, class... Args>
  typename std::enable_if_t<std::is_void_v<Return>, void
                            > timeExecution(void (*function)(Args...),
                                            std::function<void(std::string const)> logger,
                                            Args... args
                                            ) {
      return timeExecution<void, Args...>(static_cast<std::function<void(Args...)>>(function), logger, args...);
  }
}
#endif /* FUNCTIONALS_H */
