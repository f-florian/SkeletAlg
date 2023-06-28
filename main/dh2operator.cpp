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
#include "dh2operator.h"
#include "compress.h"

#include "h2lib/dcluster.h"
#include "h2lib/dh2compression.h"
#include "h2lib/h2matrix.h"
#include "h2lib/helmholtz3d.h"
#include "h2lib/truncation.h"

#include "logger/logger.h"
#include "utils/cutils.h"

#include "options.h"

namespace MainNamespace {
    static inline dcluster* cloneDcluster(dcluster const * const other  //!< Directional cluster to copy.
                                          ) {                           //! Clone a directional cluster. \return Copy of \p other.
        dcluster *data = new_dcluster(other->size, other->idx, other->sons, other->dim);
        for(size_t dimension = 0; dimension < other->dim; ++dimension) {    // Copy bounding box.
            data->bmin[dimension] = other->bmin[dimension];
            data->bmax[dimension] = other->bmax[dimension];
        }
        for(size_t son = 0; son < other->sons; ++son)       // Clone son, recursively.
            #ifdef USE_OPENMP
            #pragma omp task
            #endif
            data->son[son] = cloneDcluster(other->son[son]);
        #ifdef USE_OPENMP
        #pragma omp taskwait
        #endif
        update_dcluster(data);
        return data;
    }
    void DH2MatrixWrapper::build(helmholtz3d *const description, dcluster const *const rowCluster, dcluster const *const columnCluster, diradmdata const &admissibilityData) {
        if((rowCluster == nullptr) || (columnCluster == nullptr))
            return;
        bemDescription_ = description;
        bemDescription_->rp = progressBar();
        Logger::log(Logger::Color::yellow, Logger::Level::programFlow, true, "Building operator:");
        directionalClusters_[0] = cloneDcluster(rowCluster);
        directionalClusters_[1] = cloneDcluster(columnCluster);
        directionalBlock_ = nullptr;    // Free any existing data.
        matrix_ = buildOperator(directionalClusters_, directionalBlock_, admissibilityData, *bemDescription_);

        #ifdef ____debug
        size_t maxColumn = 0, maxRow = 0;
        for(size_t idx = 0; idx < rowCluster->size; ++idx) {
            if(rowCluster->idx[idx] > maxRow)
                maxRow = rowCluster->idx[idx];
        }
        for(size_t idx = 0; idx < columnCluster->size; ++idx) {
            if(columnCluster->idx[idx] > maxColumn)
                maxColumn = columnCluster->idx[idx];
        }
        Logger::log(Logger::Color::yellow, Logger::Level::developer, "Matrix dimension: %d x %d (%ld, %ld)\n", matrix_->rb->t->size, matrix_->cb->t->size, maxRow, maxColumn);
        #endif
    }
    void DH2MatrixWrapper::apply(field const scale, avector const *const source, avector *const destination) const {
        if(matrix_ != nullptr)
            addeval_parallel_dh2matrix(scale, *matrix_, source, destination);
    }
    void DH2MatrixWrapper::applyTransposed(field const scale, avector *const source, avector *const destination) const {
        if(matrix_ == nullptr)
            return;
        conj_avector(source);       // addevaltrans actually means adjoint in h2lib
        conj_avector(destination);  // conjugate so that the second time it is reverted for what was already there.
        addevaltrans_parallel_dh2matrix(scale, *matrix_, source, destination);
        conj_avector(source);       // revert conjugation
        conj_avector(destination);  // addevaltrans actually means adjoint in h2lib
    }
    size_t DH2MatrixWrapper::rows() const {
        if(matrix_ == nullptr)
            return 0;
        else
            return directionalClusters_[1]->size;
    }
    size_t DH2MatrixWrapper::columns() const {
        if(matrix_ == nullptr)
            return 0;
        else
            return directionalClusters_[0]->size;
    }
}
