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
#include "orient.h"

#include <cassert>
#include <cstring>

#include "h2lib/surface3d.h"

namespace MainNamespace {
    namespace Surface {
        static surface3d* copySurface(surface3d const *const surface    //!< Surface to copy.
                                      ) {                               //! Copy a surface3d object. \return A copy of \p surface.
            surface3d* copy = new_surface3d(surface->vertices, surface->edges, surface->triangles);
            #ifdef USE_OPENMP
            #pragma omp parallel
            #endif
            {
                std::memcpy(copy->x, surface->x, 3*sizeof(real)*static_cast<size_t>(surface->vertices));
                std::memcpy(copy->e, surface->e, 2*sizeof(uint)*static_cast<size_t>(surface->edges));
                std::memcpy(copy->t, surface->t, 3*sizeof(uint)*static_cast<size_t>(surface->triangles));
                std::memcpy(copy->s, surface->s, 3*sizeof(uint)*static_cast<size_t>(surface->triangles));
                copy->m = nullptr;
            }
            return copy;
        }

        surface3d* orientCopy(surface3d const *const surface, OrientingInfo const &correct) {
            assert(correct.size() == static_cast<size_t>(surface->triangles));
            surface3d *copy = copySurface(surface);
            #ifdef USE_OPENMP
            #pragma omp parallel for simd
            #endif
            for(uint triangle = 0; triangle < copy->triangles; ++triangle)
                if(!correct[triangle]) {
                    std::swap(copy->t[triangle][0], copy->t[triangle][1]);
                    std::swap(copy->s[triangle][0], copy->s[triangle][1]);
                }
            prepare_surface3d(copy);
            return copy;
        }
    }
}
