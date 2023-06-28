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
#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <array>
#include <complex>
#include <numeric>

#include "logger/logger.h"

#include "constants.h"

namespace MainNamespace {
    template <typename Field> using SpaceVector = std::array<Field, spaceDimension>;      //!< Vectors on spaceDimension dimensional vector space on the field Field.
    template <typename Field> constexpr inline SpaceVector<Field> pointer2constVector(std::add_const_t<Field> const coordinate[3]) {
        return {coordinate[0],coordinate[1],coordinate[2]};
    }
    template <typename Field> constexpr inline SpaceVector<std::remove_const_t<Field>> pointer2vector(std::add_const_t<Field> const coordinate[3]) {
        return {coordinate[0],coordinate[1],coordinate[2]};
    }
    template <typename Field> constexpr inline SpaceVector<std::add_const_t<Field>> makeConst(SpaceVector<Field> const vector) {
        return {vector[0],vector[1],vector[2]};
    }
    template <typename Field> inline SpaceVector<Field>& operator-= (SpaceVector<Field> &v1, SpaceVector<Field> const v2) {
        #pragma omp simd
        for(size_t idx = 0; idx < spaceDimension; ++idx)
            v1[idx] -= v2[idx];
        return v1;
    }
    template <typename Field> inline SpaceVector<std::remove_const_t<Field>> operator- (SpaceVector<Field> const v1, SpaceVector<Field> const v2) {
        SpaceVector<std::remove_const_t<Field>> aux;
        #pragma omp simd
        for(size_t idx = 0; idx < spaceDimension; ++idx)
            aux[idx] = v1[idx] - v2[idx];
        return aux;
    }
    template <typename Field> inline SpaceVector<Field>& operator+= (SpaceVector<Field> &v1, SpaceVector<Field> const v2) {
        #pragma omp simd
        for(size_t idx = 0; idx < spaceDimension; ++idx)
            v1[idx] += v2[idx];
        return v1;
    }
    template <typename Field> inline SpaceVector<std::remove_const_t<Field>> operator+ (SpaceVector<Field> const v1, SpaceVector<Field> const v2) {
        SpaceVector<std::remove_const_t<Field>> aux;
        #pragma omp simd
        for(size_t idx = 0; idx < spaceDimension; ++idx)
            aux[idx] = v1[idx] + v2[idx];
        return aux;
    }
    template <typename Field> inline SpaceVector<std::remove_const_t<Field>> operator/ (SpaceVector<Field> const vector, Field const scale) {
        SpaceVector<std::remove_const_t<Field>> aux;
        #pragma omp simd
        for(size_t idx = 0; idx < spaceDimension; ++idx)
            aux[idx] = vector[idx] / scale;
        return aux;
    }
    template <typename Field1, typename Field2> inline SpaceVector<std::remove_const_t<Field2>> operator* (SpaceVector<Field1> const vector, Field2 const scale) {
        SpaceVector<std::remove_const_t<Field2>> aux;
        #pragma omp simd
        for(size_t idx = 0; idx < spaceDimension; ++idx)
            aux[idx] = vector[idx] * scale;
        return aux;
    }
}

#endif /* LINEARALGEBRA_H */
