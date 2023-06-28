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
#include "functions.h"

#include <cstddef>
#include <execution>
#include <functional>
#include <numeric>

#include "h2lib/avector.h"
#include "h2lib/helmholtz3d.h"
#include "h2lib/krylov.h"
#include "h2lib/sparsematrix.h"

#include "linearalgebra.h"
#include "logger/logger.h"
#include "options.h"

namespace MainNamespace {
    namespace Functions {
        PlaneWave::PlaneWave(SpaceVector<double const> const frequency, field const waveNumber) {
            waveNumber_ = waveNumber;
            real normalizationReal = 0;
            for(size_t j = 0; j < spaceDimension; ++j)
                normalizationReal += frequency[j] * frequency[j];

            field normalization = waveNumber_ / sqrt(normalizationReal);
            for(size_t j = 0; j < spaceDimension; ++j)
                frequency_[j] = frequency[j] * normalization;
            Logger::log(Logger::Level::debug, "Normalizing frequencies (scaling factor: %lf%+lfi; frequencies: %lf%+lfi, %lf%+lfi, %lf%+lfi)\n", creal(normalization), cimag(normalization), creal(frequency_[0]), cimag(frequency_[0]), creal(frequency_[1]), cimag(frequency_[1]), creal(frequency_[2]), cimag(frequency_[2]));
        }
        field PlaneWave::dirichlet(SpaceVector<real const> const point, real const scale) const {
            field exponent = 0;
            for(size_t j = 0; j < spaceDimension; ++j){
                exponent += frequency_[j]*point[j];
            }
            exponent *= scale;
            return cexp(exponent) - cexp(-exponent);
        }
        field PlaneWave::neumann(SpaceVector<real const> const point, SpaceVector<real const> const normal, real const scale) const {
            field exponent = 0;
            for(size_t j = 0; j < spaceDimension; ++j){
                exponent += frequency_[j]*point[j];
            }
            exponent *= scale;
            field direction = 0;
            for(size_t j = 0; j < spaceDimension; ++j){
                direction += frequency_[j]*normal[j];
            }
            return (cexp(exponent) + cexp(-exponent)) * direction * scale;
        }

        MeshOffset::MeshOffset(Reference const *const solution)
            : solution_(solution),
          traceScaling_(csqrt(solution_->waveNumber())) {}
        void MeshOffset::c0Dirichlet(Mesh const &mesh, std::vector<uint> const &indices, bool const dof, real const waveScale, avector *const coefficients) const {
            #ifdef USE_OPENMP
            #pragma omp parallel for simd
            #endif
            for(size_t vertex = 0; vertex < indices.size(); ++vertex) {
                uint index = indices[vertex];
                uint shift = index < mesh.vertexDofs_ ? 0 : mesh.triangleDofs_;
                setentry_avector(coefficients, index + shift, dirichlet(dof, mesh.vertex(index), waveScale) * traceScaling_);
            }
        }
        void MeshOffset::l2Neumann(helmholtz3d const *const parameters, Mesh const &mesh, std::vector<uint> const &indices, bool const dof, real const waveScale, real const stiffness, avector *const coefficients) const {
            field finalScaling = 2 / traceScaling_ * stiffness;
            #ifdef USE_OPENMP
            #pragma omp parallel for simd
            #endif
            for(size_t triangle = 0; triangle < indices.size(); ++triangle) {
                uint index = indices[triangle];
                uint shift = index < mesh.triangleDofs_ ? mesh.vertexDofs_ : mesh.vertices();
                SpaceVector<real const> normal = mesh.normal(index);
                std::array<size_t,3> vertexIdxs = mesh.triangle(index);
                        
                field sum = 0.0;
                for(size_t nodeIdx = 0; nodeIdx < parameters->nq; ++nodeIdx) {
                    real barycentricWeights[3] = {1-parameters->xq[0][nodeIdx]-parameters->xq[0][nodeIdx], parameters->xq[0][nodeIdx], parameters->xq[0][nodeIdx]};
                    SpaceVector<real> node = mesh.vertex(vertexIdxs[0])*barycentricWeights[0];
                    node += mesh.vertex(vertexIdxs[1])*barycentricWeights[1];
                    node += mesh.vertex(vertexIdxs[2])*barycentricWeights[2];

                    sum += neumann(dof, {node[0], node[1], node[2]} , normal, waveScale) * parameters->wq[nodeIdx];
                }
                setentry_avector(coefficients, indices[triangle] + shift, sum * finalScaling);
            }
        }
        field MeshOffset::dirichlet(SpaceVector<real const> const point, real const scale) const {
            return solution_->dirichlet(point, scale);
        }
        field MeshOffset::neumann(SpaceVector<real const> const point, SpaceVector<real const> const normal, real const scale) const {
            return solution_->neumann(point, normal, scale);
        }

        field ZeroExtension::dirichlet(bool const dof, SpaceVector<real const> const point, real const scale) const {
            if(dof)
                return 0;
            return solution_->dirichlet(point, scale);
        }
        field ZeroExtension::neumann(bool const dof, SpaceVector<real const> const point, SpaceVector<real const> const normal, real const scale) const {
            if(dof)
                return 0;
            return solution_->neumann(point, normal, scale);
        }

    }
}
