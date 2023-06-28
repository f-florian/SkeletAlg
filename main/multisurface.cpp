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
#include "multisurface.h"

#include <unistd.h>

#include <array>
#include <cstdio>
#include <cstring>
#include <execution>
#include <functional>
#include <numeric>
#include <tuple>
#include <utility>

#include "h2lib/amatrix.h"
#include "h2lib/helmholtz3d.h"
#include "h2lib/krylov.h"
#include "h2lib/surface3d.h"

#include "logger/logger.h"
#include "utils/cutils.h"
#include "utils/data.h"
#include "utils/functionals.h"
#include "linearalgebra.h"
#include "meshreader.h"
#include "options.h"
#include "spaces.h"
#include "surface.h"
#include "surfacehelpers.h"

namespace MainNamespace {
    namespace Surface {
        typedef void (*StepFunction)(field, void const *, avector const *const, avector *const);
        using ChildFunction = std::function<void(SingleSurface &)>;
        using ChildConstFunction = std::function<void(SingleSurface const &)>;
        static constexpr size_t gmresMaxIterations = 20000;
        static constexpr real kstiffness[] = {1,1}, kdensity[] = {1,5};

        static void getFromFile(std::string const filename,     //!< Name of the mesh file.
                                MeshReader &meshReader          //!< Gmsh interface object. Is modified, as it saves the elements permutations.
                                ) {                             //! Load the mesh from a given file. \return A Mesh.
            Logger::log(Logger::Color::green, Logger::Level::debugData, false, "Reading file %s: ", filename.c_str());
            if(access(filename.c_str(), R_OK)) {
                Logger::log(Logger::Color::red, Logger::Level::warning, true, "File does not exist, or unreadable.");
                abort();
            } else {
                Logger::log(Logger::Color::green, Logger::Level::developer, true, "File readable.");
            }
            Logger::indent();
            meshReader.parseMesh(filename);
            prepare_surface3d(meshReader.mesh().mesh_);
            Logger::deindent();

            Logger::log(Logger::Color::green, Logger::Level::developer, true, "Done loading mesh with (v, e, t): (%ld, %ld, %ld); sizes: D: %ld, N: %ld.", meshReader.mesh().vertices(), meshReader.mesh().mesh_->edges, meshReader.mesh().triangles(), meshReader.mesh().vertexDofs_, meshReader.mesh().triangleDofs_);
            meshReader.mesh().mesh_->m=NULL;
        }

        static std::pair<size_t, avector*> doGmres(avector const *const rhs,            //!< Right-hand side of the linear system.
                                                   StepFunction const stepFunction,     //!< Pointer to a function performing the matrix-vector product.
                                                   void const *const data               //!< Additional input to \p stepFunction holding data needed for the evaluation (currently a pointer to MultiSurface).
                                                   ) {                                  //! Solve the linear system with right-hand side \p rhs and matrix-vector product given by \p stepFunction. \return Number of iterations and solution.
            uint k;
            avector *solution = new_avector(rhs->dim); clear_avector(solution);
            DumbPointer<avector> rhat = new_avector(solution->dim);
            DumbPointer<avector> q = new_avector(solution->dim);
            DumbPointer<amatrix> qr = new_amatrix(solution->dim, dimKrylov());
            DumbPointer<avector> tau = new_avector(dimKrylov());

            double tolerance = norm2_avector(rhs) * epsKrylov();
            Logger::log("Running GMRES, tolerance %.2e (%.2e)\n", tolerance, epsKrylov());

            size_t gmresIteration=1;    // One iteration always happens, the increment happens after the error log. 
            startTiming();
            init_gmres(stepFunction, data, rhs, solution, *rhat, *q, &k, *qr, *tau);
            for(; gmresIteration < gmresMaxIterations; ++gmresIteration) {
                step_gmres(stepFunction, data, rhs, solution, *rhat, *q, &k, *qr, *tau);
                double errorrr = residualnorm_gmres(*rhat, k);
                if(errorrr < tolerance)
                    break;
                Logger::log(Logger::Level::dataOperation, "\r  Step %4zd (%2u): Residual norm %.2e/%.2e", gmresIteration, k, errorrr, tolerance);
            }
            finish_gmres(stepFunction, data, rhs, solution, *rhat, *q, &k, *qr, *tau);
            if(gmresIteration < gmresMaxIterations)     // If the break is never hit, then the for condition is violated. 
                Logger::log(Logger::Color::green, Logger::Level::programFlow, "Residual norm: %.2e (%s seconds, %ld iterations).\n", residualnorm_gmres(*rhat, k), stopTiming().c_str(), gmresIteration);
            else
                Logger::log(Logger::Color::red, Logger::Level::warning, "Residual norm: %.2e (%s seconds), did not converge in %ld iterations!\n", residualnorm_gmres(*rhat, k), stopTiming().c_str(), gmresIteration - 1);

            return std::make_pair(gmresIteration, solution);
        }
        
        static void forwardChildren(MultiSurface::Children &children,   //!< Children list.
                                    ChildFunction function              //!< Function to apply to all \p children.
                                    ) {                                 //! Apply the \p function to all \p children.
            for(SingleSurface &child : children)
                function(child);
        }
        static void forwardChildrenC(MultiSurface::Children const &children,     //!< Children (constant) list.
                                    ChildConstFunction function                 //!< Function to apply to all (constant) \p children.
                                    ) {                                         //! Apply the (constant) \p function to all \p children.
            for(SingleSurface const &child : children)
                function(child);
        }

        MultiSurface::MultiSurface(std::string const meshFileName) {
            Logger::log(Logger::Color::green, Logger::Level::programFlow, true, "Loading the mesh...");
            Logger::indent();
            getFromFile(meshFileName, gmshInterface_);                                                                          // Load mesh
            metric_ = new_helmholtz3d(geometry(), 0, static_cast<uint>(quadraturePoints()), static_cast<uint>(interpolationPoints()));  // Metric data
            Logger::deindent();

            Logger::log(Logger::Color::green, Logger::Level::programFlow, true, "Creating %ld children", mesh().membership_.size());
            children_.reserve(mesh().membership_.size());
            for(size_t surfaceIdx = 0; surfaceIdx < mesh().membership_.size(); ++surfaceIdx) {
                Logger::log(Logger::Color::green, Logger::Level::objectManagement, "Building child...");
                Logger::indent();
                children_.emplace_back(mesh(), surfaceIdx, kstiffness[surfaceIdx], kdensity[surfaceIdx], &gmshInterface_);
                Logger::deindent();
            }
        }
        MultiSurface::~MultiSurface() {
            maybeFree(metric_);
        }
        Result MultiSurface::solve(AdmissibilityParameters const admissibilityParameters, Functions::MeshOffset &offset) {
            waveNumber(admissibilityParameters.waveNumber);
            metric_->kappa = waveNumber();
            diradmdata admissibilityData = {admissibilityParameters.eta1, admissibilityParameters.eta2, admissibilityParameters.eta3, waveNumber(), nullptr, nullptr};

            DumbPointer<avector> rhs = new_avector(static_cast<uint>(problemSize())); // allocate rhs
            clear_avector(*rhs);

            forwardChildren(children_, [this, &offset](SingleSurface &child){child.makeOffset(offset, metric_);});
            forwardChildren(children_, [admissibilityData, &rhs](SingleSurface &child){child.makeMultiRhs(admissibilityData, *rhs);});
            gmshInterface_.saveView(MeshReader::View::rhs, *rhs);
            forwardChildren(children_, [admissibilityData](SingleSurface &child){child.buildOperator(admissibilityData);});
            forwardChildrenC(children_, [&rhs](SingleSurface const &child){child.makeSingleRhs(*rhs);});
            gmshInterface_.saveView(MeshReader::View::rhs, *rhs);

            auto stepFunction = [](field scaling, void const * data, avector const * source, avector * destination) {
                MultiSurface const * const surface = static_cast<MultiSurface const*>(data);
                surface->apply(scaling, source, destination);
            };
            size_t iterations; DumbPointer<avector> solution; std::tie(iterations, solution) = doGmres(*rhs, stepFunction, static_cast<void const*>(this));     // Don't clear solution, it's done inside gmrs. Besides, it's nullptr before.
            return std::make_pair(iterations, this->error(*solution, offset));
        }
        ErrorVector MultiSurface::error(avector *const solution, Functions::MeshOffset const &reference) const {
            gmshInterface_.saveView(MeshReader::View::solution, solution);
            Logger::log(Logger::Color::magenta, Logger::Level::debug, true, "Computing error");
            real solutionNorm = 0, errorNorm = 0;

            for(decltype(children_.begin()) child = children_.begin(); child != children_.end(); ++child) {
                real childSolution, childError;
                std::tie(childSolution, childError) = child->computeDifferenceNorm(solution, reference, metric_);
                solutionNorm += childSolution;
                errorNorm += childError;
            }
            return sqrt(errorNorm/solutionNorm);
        }
        void MultiSurface::apply(field const scaling, avector const * const source, avector * const destination) const {
            for(SingleSurface const &child : children_)
                child.apply(scaling, source, destination);
        }
        std::pair<double, double> MultiSurface::properties() const {
            double hmin, hmax, anglemin, angleedge;
            getproperties_surface3d(geometry(), &hmin, &hmax, &anglemin, &angleedge);
            Logger::log("Global mesh properties: %u vertices, %u edges, %u triangles; hmin %.2e, hmax %.2e, anglemin %g, angleedge %g\n",
                        geometry()->vertices, geometry()->edges, geometry()->triangles, hmin, hmax, anglemin, angleedge);
            for(SingleSurface const &child: children_)
                child.properties();
            return std::make_pair(hmax, hmax/hmin);
        }
    }
}
