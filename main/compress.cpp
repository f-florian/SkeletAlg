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
#include <functional>

#include "h2lib/dcluster.h"
#include "h2lib/dh2compression.h"
#include "h2lib/h2matrix.h"
#include "h2lib/helmholtz3d.h"
#include "h2lib/truncation.h"

#include "logger/logger.h"
#include "utils/cutils.h"
#include "utils/functionals.h"

#include "options.h"
#include "surfacehelpers.h"

#include "compress.h"

namespace MainNamespace {
    static inline void reportClusterCompression(dclusterbasis const * const basis,              //!< Cluster basis. \todo Why? 
                                                dclusteroperator const * const clOperator,      //!< Cluster operator. \todo Why?  
                                                size_t const dof,                               //!< Degrees of freedom, to report the memory used per degree of freedom.
                                                std::string const timingReport                  //!< List of the times to build the ??, as reported by varius timers. \todo build what?
                                                ) {                                             //! \brief Log a message reporting the compression of the given cluster.
      size_t sz = getsize_dclusterbasis(basis);
      size_t szo = getsize_dclusteroperator(clOperator);
      Logger::log("%.1f MB (%.2f KB/DoF); for basis change: %.1f MB (%.2f KB/DoF); %s seconds;"
                  "Maximal rank %u, Rank sum %u, Rank branch sum %u\n",
                  bytes2MB(sz), bytes2KB(sz) / static_cast<double>(dof), bytes2MB(szo), bytes2KB(szo) / static_cast<double>(dof), timingReport.c_str(),
                  getmaxrank_dclusterbasis(basis), basis->ktree, basis->kbranch);
      
      uint rankLevels, *rank;
      getlevelranks_dclusterbasis(basis, &rankLevels, &rank);
      Logger::log("  Ranks: [");
      for(size_t idx=0; idx<=rankLevels; idx++)
        Logger::log(",%u", rank[idx]);
      Logger::log("]\n");
      freemem(rank);
    }

    static inline dclusteroperator *buildImplicitClusterBasis(helmholtz3d *const bm,            //!< Pointer to the helmoltz object. \todo why? what?
                                                              dcluster const *const root,       //!< Root of the directional cluster.
                                                              leaf_func_t leaf,                 //!< idk \todo idk.
                                                              size_t const dofs                 //!< Number of degrees of freedom (to report the memory per degee of freedom).
                                                              ) {                               //! \brief Build a cluster basis. \todo why implicit? \return directional cluster operater.
        return timeExecution(buildweight_implicit_dclusterbasis,
                             [dofs](std::string const time, dclusteroperator const * const clusterBasis) {
                                 uint sz = getsize_dclusteroperator(clusterBasis);
                                 Logger::log(" %.1f MB (%.2f KB/DoF), %s;", bytes2MB(sz), kBRatio(sz, dofs), time.c_str());
                             },
                             root,
                             (rank_func_t) bm->rank,
                             leaf,
                             (transfer_func_t) bm->transfer,
                             static_cast<void*>(bm), progressBar());
    }
    
    static inline void copy(diradmdata const &data,     //!< Directional admissibility data source.
                            diradmdata &destination     //!< [out] Destination.
                            ) {                         //! The field xy is (re)allocated, ld is cleared and left as a null pointer. Parameters are copied. \brief Copy a diradmdata stuct.
      destination.eta1 = data.eta1;
      destination.eta2 = data.eta2;
      destination.eta3 = data.eta3;
      destination.kappa = data.kappa;

      maybeFree(destination.xy);
      destination.xy = allocreal(3);    // This must be so for some reason.
      maybeFree(destination.ld);
      destination.ld = nullptr;
    }
    static inline void init(diradmdata &destination     //!< [out] Object to initialize.
                            ) {                         //! Set pointers to null to prevent deallocation of non-allocated memory.
        destination.xy = nullptr;
        destination.ld = nullptr;
    }
    static inline dh2matrix* compressRowColumnFar(dblock const * const root,
                                                  helmholtz3d * const bemDescription,
                                                  dclusteroperator const * const rowCluster,
                                                  dclusteroperator const * const columnCluster,
                                                  truncmode const * const trMode
                                                  ) {
      startTiming();
      Logger::log(Logger::Color::green, Logger::Level::programFlow, "Compress row column far");
      dclusteroperator *rowoperator = buildfromcluster_dclusteroperator(root->rc);
      dclusterbasis *rowbasis = buildrowbasis_implicit_dh2matrix(root, (rank_func_t) bemDescription->rank, (leaf_func_t) bemDescription->row_leaf, (transfer_func_t) bemDescription->transfer, (coupling_func_t) bemDescription->coupling, bemDescription, rowCluster, columnCluster, trMode, eps(), rowoperator, progressBar());
      Logger::log(Logger::Color::green, Logger::Level::programFlow, "Compress row cluster basis, eps=%.2e: ", eps());
      reportClusterCompression(rowbasis, rowoperator, bemDescription->gr->triangles, stopTiming());
      
      startTiming();
      dclusteroperator *columnoperator = buildfromcluster_dclusteroperator(root->cc);
      dh2matrix *dh2Matrix = buildcolmatrix_implicit_dh2matrix(root, rowbasis, rowoperator, (rank_func_t) bemDescription->rank, (leaf_func_t) bemDescription->col_leaf, (transfer_func_t) bemDescription->transfer, (coupling_func_t) bemDescription->coupling, bemDescription, columnCluster, trMode, eps(), columnoperator, progressBar());
      Logger::log(Logger::Color::green, Logger::Level::programFlow, "Compress column cluster basis and filling farfield: ");
      reportClusterCompression(dh2Matrix->cb, columnoperator, bemDescription->gr->triangles, stopTiming());

      del_dclusteroperator(columnoperator);
      del_dclusteroperator(rowoperator);
      return dh2Matrix;
    }
    static void makeDirections(std::array<DumbPointer<dcluster>,2> &directionalClusters,
                               diradmdata const &directions,
                               DumbPointer<dblock> &directionalBlock
                               ) {
        diradmdata directionsLocal; init(directionsLocal); copy(directions, directionsLocal);
        Logger::log(Logger::Color::yellow, Logger::Level::debug, "Creating directions (η₁/𝕽𝔢 k: %lf)", directionsLocal.eta1 / abs(creal(directionsLocal.kappa)));
        directionsLocal.ld = builddirections_box2_dcluster(*directionalClusters[0], *directionalClusters[1], directionsLocal.eta1 / abs(creal(directionsLocal.kappa)));
        uint directionCount = getalldirections_dcluster(*directionalClusters[0]);
        Logger::log(Logger::Color::green, Logger::Level::programFlow, "%u directions (%.1f per cluster);\n", directionCount, static_cast<real>(directionCount) / static_cast<real>(directionalClusters[0]->desc));
        for(size_t indexCounter = 0; indexCounter <= directionsLocal.ld->depth; indexCounter++)
          Logger::log(Logger::Color::green, Logger::Level::debug, "Level %2zd: maxdiam %.3e,%3zu directions\n", indexCounter, directionsLocal.ld->maxdiam[indexCounter], directionsLocal.ld->directions[indexCounter]);
        
        directionalBlock = build_dblock(*directionalClusters[0], *directionalClusters[1], 0, parabolic_admissibility, &directionsLocal);
        Logger::log(Logger::Color::green, Logger::Level::programFlow, "Directional block trees: (eta1=%.1f, eta2=%.1f, eta3=%.1f): %u blocks (%.1f per cluster). Used directions row, column (per cluster):",
                    directionsLocal.eta1,
                    directionsLocal.eta2,
                    directionsLocal.eta3,
                    directionalBlock->desc,
                    static_cast<real>(directionalBlock->desc) / static_cast<real>(directionalClusters[0]->desc));
        markused_dblock(*directionalBlock);
        for(size_t idx = 0; idx < 2; ++idx) {
            propagateused_dcluster(*directionalClusters[idx]);
        }
        for(size_t idx = 0; idx < 2; ++idx) {
            uint directions = getuseddirections_dcluster(*directionalClusters[idx]);
            Logger::log("%u (%.1f);", directions, static_cast<real>(directions) / static_cast<real>(directionalClusters[idx]->desc));
        }
        // maybeFree(directionsLocal.ld);
        // maybeFree(directionsLocal.xy);
    }
    dh2matrix* buildOperator(std::array<DumbPointer<dcluster>, 2> &directionalClusters, DumbPointer<dblock> &directionalBlock, diradmdata const &directions, helmholtz3d *const bemDescription) {
        makeDirections(directionalClusters, directions, directionalBlock);    // Build directions and related data.
        truncmode *truncationParameters = new_blockreleucl_truncmode();   // Truncation strategy, don't change.
        Logger::log(Logger::Color::green, Logger::Level::programFlow, "\nCompressing (non directional) cluster bases. Implicit (row,column) basis weights:");
        std::array<DumbPointer<dclusteroperator>, 2> clusterOperators;
        clusterOperators[0] = buildImplicitClusterBasis(bemDescription, *directionalClusters[0], (leaf_func_t) bemDescription->row_leaf, bemDescription->gr->triangles);
        clusterOperators[1] = buildImplicitClusterBasis(bemDescription, *directionalClusters[1], (leaf_func_t) bemDescription->col_leaf, bemDescription->gr->triangles);

        Logger::log(Logger::Color::green, Logger::Level::programFlow, "Building (compressed) boundary integral operators...\n");
        dh2matrix* matrix = compressRowColumnFar(*directionalBlock, bemDescription, *clusterOperators[0], *clusterOperators[1], truncationParameters);
        timeExecution<void, helmholtz3d const*, dh2matrix*
                      >(std::function<void(helmholtz3d const *, dh2matrix *)>(fill_near_dh2matrix_helmholtz3d),
                        [](std::string const time) { Logger::log("Nearfield: %s; ", time.c_str());},
                        bemDescription,
                        matrix);

        del_truncmode(truncationParameters);
        size_t fullSize = getsize_dh2matrix(matrix);
        size_t nearSize = getnearsize_dh2matrix(matrix);
        size_t dof = bemDescription->gr->triangles + bemDescription->gr->vertices;
        Logger::log("%.1f MB (%.2f KB/DoF); %.1f MB nearfield (%.2f KB/DoF, %.1f%%).\n", bytes2MB(fullSize), kBRatio(fullSize, dof), bytes2MB(nearSize), kBRatio(nearSize, dof), percentRatio(nearSize, fullSize));

        if(reference()) {
            amatrix *reference = new_amatrix(matrix->rb->t->size, matrix->cb->t->size);
            fill_block_amatrix_helmholtz3d(bemDescription, matrix->rb->t, matrix->cb->t, reference);
            Logger::log("%.3e\n", norm2diff_amatrix_nopermute_dh2matrix(matrix, reference)/norm2_parallel_amatrix(reference));
            del_amatrix(reference);
        }
        return matrix;
    }
}
