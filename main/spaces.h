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
#ifndef SPACES_H
#define SPACES_H

#include <array>
#include <vector>

#include "logger/logger.h"

#include "utils/arrays.h"
#include "utils/enum.h"

namespace MainNamespace {
    enum class SpaceSlices {    //! Enum to manage arrays of differend boundary operator or associated auxiliary elements.
        DirichletSolution = 0,  //!< Elements with given boundary data, in the Dirichlet trace.
        DirichletData,          //!< Element with unknown (solution) data, in the Dirichlet trace.
        NeumannSolution,        //!< Elements with given boundary data, in the Neumann trace.
        NeumannData,            //!< Element with unknown (solution) data, in the Neumann trace.
        COUNT                   //!< Counting the number of elements in the enumerator.
    };
    enum class OperatorSlices {         //! Enum to manage arrays of differend boundary operator or associated auxiliary elements.
        SolutionSolutionSingle = 0,     //!< Object associated to the sigle layer operator, with input and output elements with given boundary data.
        SolutionSolutionDouble,         //!< Object associated to the double layer operator, with input and output elements with given boundary data.
        SolutionSolutionHypersingular,  //!< Object associated to the hypersingular operator, with input and output elements with given boundary data.
        DataSolutionSingle,             //!< Object associated to the single layer operator, with input an element with given boundary data, output an element with unknown (solution) data.
        COUNTDof = DataSolutionSingle,  //!< Count of the slices that only operate on degrees of freedom.
        DataSolutionDouble,             //!< Object associated to the double layer operator, with input an element with given boundary data, output an element with unknown (solution) data.
        SolutionDataDouble,             //!< Object associated to the double layer operator, with input output an element with unknown (solution) data, an element with given boundary data (for K').
        DataSolutionHypersingular,      //!< Object associated to the hypersingular operator, with input an element with given boundary data, output an element with unknown (solution) data.
        COUNT                           //!< Counting the number of elements in the enumerator.
    };
    enum class DataSpaces {     //! Enum to manage arrays of differend boundary operator or associated auxiliary elements.
        Dirichlet = 0,          //!< Element with unknown (solution) data, in the Dirichlet trace.
        Neumann,                //!< Element with unknown (solution) data, in the Neumann trace.
        COUNT                   //!< Counting the number of elements in the enumerator.
    };
    enum class DataOperator {           //! Enum to manage arrays of differend boundary operator or associated auxiliary elements.
        DataSolutionSingle = 0,         //!< Object associated to the single layer operator, with input an element with given boundary data, output an element with unknown (solution) data.
        DataSolutionDouble,             //!< Object associated to the double layer operator, with input an element with given boundary data, output an element with unknown (solution) data.
        SolutionDataDouble,             //!< Object associated to the double layer operator, with input output an element with unknown (solution) data, an element with given boundary data (for K').
        DataSolutionHypersingular,      //!< Object associated to the hypersingular operator, with input an element with given boundary data, output an element with unknown (solution) data.
        COUNT                           //!< Counting the number of elements in the enumerator.
    };
  
    template<class Objects> using IndexSpace = std::array<std::vector<Objects>, size<SpaceSlices>()>;   //!< Type for an array with indices information.
    template<class Objects> using SlicedSpace = PointerArray<Objects, size<SpaceSlices>()>;             //!< Type for an array with pointers to all slices of some space-related object.
    template<class Objects> using SlicedOperator = PointerArray<Objects, size<OperatorSlices>()>;       //!< Type for an array with pointers to all slices of some operator-related object.
    template<class Objects> using Operators = std::array<Objects, size<OperatorSlices>()>;              //!< Type for an array with all slices of some operator-related object.
    constexpr inline std::pair<SpaceSlices, SpaceSlices> operatorSpaces(OperatorSlices const slice) {
        switch (slice) {
        case OperatorSlices::SolutionSolutionSingle:
            return {SpaceSlices::NeumannSolution, SpaceSlices::NeumannSolution};
        case OperatorSlices::DataSolutionSingle:
            return {SpaceSlices::NeumannSolution, SpaceSlices::NeumannData};
        case OperatorSlices::SolutionSolutionDouble:
            return {SpaceSlices::NeumannSolution, SpaceSlices::DirichletSolution};
        case OperatorSlices::DataSolutionDouble:
            return {SpaceSlices::NeumannSolution, SpaceSlices::DirichletData};
        case OperatorSlices::SolutionDataDouble:
            return {SpaceSlices::NeumannData, SpaceSlices::DirichletSolution};
        case OperatorSlices::SolutionSolutionHypersingular:
            return {SpaceSlices::DirichletSolution, SpaceSlices::DirichletSolution};
        case OperatorSlices::DataSolutionHypersingular:
            return {SpaceSlices::DirichletSolution, SpaceSlices::DirichletData};
        case OperatorSlices::COUNT:
            Logger::log(Logger::Color::red, Logger::Level::error, "Invalid operator slice selected (COUNT)\n");
            break;
        default:
            Logger::log(Logger::Color::red, Logger::Level::error, "Unknown operator slice selected, check your source!\n");
            break;
        }
        return {SpaceSlices::COUNT, SpaceSlices::COUNT};
    }
    constexpr inline std::pair<SpaceSlices, SpaceSlices> operatorSpaces(size_t const slice) {
        return operatorSpaces(static_cast<OperatorSlices>(slice));
    }
}

#endif /* SPACES_H */

