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
#ifndef DH2OPERATOR_H
#define DH2OPERATOR_H

#include <array>

#include "utils/cutils.h"

namespace MainNamespace {
    /**
     * Wrapper for a Directional H2 matrix and its associated data.
     */
    class DH2MatrixWrapper {
    public:
        DH2MatrixWrapper() = default;
        DH2MatrixWrapper(DH2MatrixWrapper const &other) = delete;
        DH2MatrixWrapper(DH2MatrixWrapper &&other) noexcept = default;
        size_t rows() const;
        size_t columns() const;
        void build(helmholtz3d *const description,      //!< Data structure describing some parameters and the operator to approximate.
                   dcluster const *const rowCluster,    //!< Row directional cluster.
                   dcluster const *const columnCluster, //!< Column directional cluster.
                   diradmdata const &admissibilityData  //!< Admissibility and direction parameters.
                   );                                   //!< Take ownership of, the pointer to the bem description and build the matrix and all needed data structures.
        void apply(field const scale,           //!< Scale.
                   avector const *const source, //!< Input.
                   avector *const destination   //!< [out] Vector to add to.
                   ) const;                     //!< Apply the operator to \p source * \p scale, and add the result to \p destination.
        void applyTransposed(field const scale,         //!< Scale.
                             avector *const source,     //!< Input. Not const because it needs to be conjugated, but this is reverted before the function returns.
                             avector *const destination //!< [out] Vector to add to.
                             ) const;                   //!< Apply the transposed of the operator to \p source * \p scale, and add the result to \p destination.
    
    private:
        DumbPointer<dh2matrix> matrix_; //!< The matrix
        DumbPointer<helmholtz3d> bemDescription_;       //!< BEM descriptions: information on geometry and how to construct the operators depending on the paratemers and the wave number.
        DumbPointer<dblock> directionalBlock_;          //!< Root of directional block tree.
        // std::array<DumbPointer<cluster>, 2> clusters_ = {nullptr, nullptr};             //!< Roots of nondirectional cluster trees.
        std::array<DumbPointer<dcluster>, 2> directionalClusters_ = {nullptr, nullptr}; //!< Roots of directional cluster trees.

        // diradmdata directionalAdmissibilityData_ = {0, 0, 0, 0, nullptr, nullptr};      //!< Data for clustering: admissibility parameters and wave number.
    };
}

#endif /* DH2OPERATOR_H */
