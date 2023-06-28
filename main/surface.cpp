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
#include <bits/iterator_concepts.h>
#include <cassert>
#include <cstddef>

#include "h2lib/clustergeometry.h"
#include "h2lib/dblock.h"
#include "h2lib/dclusterbasis.h"
#include "h2lib/dh2compression.h"
#include "h2lib/helmholtz3d.h"
#include "h2lib/surface3d.h"

#include "logger/logger.h"
#include "utils/cutils.h"
#include "utils/enum.h"
#include "utils/functionals.h"

#include "dh2operator.h"
#include "options.h"
#include "meshreader.h"
#include "orient.h"
#include "spaces.h"
#include "surfacehelpers.h"

#include "surface.h"

namespace MainNamespace {
    namespace Surface {
        typedef helmholtz3d (*helmholtzConstructor) (surface3d const *gr, field kappa, uint nq, uint ni);
        static real energyNorm(avector const *const vector,             //!< Solution vector. Only the first component, corresponding to single trace degrees of freedom, are considered.
                               SingleSurface const *const operators,    //!< Object capable of applying the Calderón operator. 
                               field const waveNumber                   //!< Wave number.
                               ) {                                      //! Approximate the energy norm of \p vector.
            avector* tmpCalderon = new_avector(vector->dim); clear_avector(tmpCalderon);// Auxiliary vector to store the result of the Calderón operator.
            field coercivityScaling = 1;              // conj(s)/abs(s) is implicit in the choice of the test function. Compare thesis with Banjai Lubic Sayas (doi  10.1007/s00211-014-0650-0)
            operators->apply(coercivityScaling, vector, tmpCalderon);
            real norm = creal(dotprod_avector(vector, tmpCalderon));
            Logger::log(Logger::Color::blue, Logger::Level::debug, true, "Energy norm %lf, scaling %.2lf%+.2lfi", norm, creal(coercivityScaling), cimag(coercivityScaling));
            del_avector(tmpCalderon);
            return norm;
        }
        static void makeTmpApply(field const scaling,                   //!< Scaling for the operator.
                                 helmholtz3d *const description,        //!< Data structure describing some parameters and the operator to approximate.
                                 dcluster const *const rowCluster,      //!< Row directional cluster.
                                 dcluster const *const columnCluster,   //!< Column directional cluster.
                                 diradmdata const &admissibilityData,   //!< Admissibility and direction parameters.
                                 avector *const source,                 //!< Input.
                                 avector *const destination,            //!<[in,out] output, but add to it instead of replacing.
                                 bool const transpose = false           //!< Whether to transpose the operator before applying it.
                                 ) {                                    //! Costruct a disposable directional matrix and apply it.
            DH2MatrixWrapper matrix;
            if(transpose) {
                matrix.build(description, columnCluster, rowCluster, admissibilityData);
                Logger::log(Logger::Color::green, Logger::Level::debug, "(trans:) source %d, destination %d, row %d column %d\n", source->dim, destination->dim, matrix.rows(), matrix.columns());
                matrix.applyTransposed(scaling, source, destination);
            } else {
                matrix.build(description, rowCluster, columnCluster, admissibilityData);
                Logger::log(Logger::Color::green, Logger::Level::debug, "source %d, destination %d row %d column %d\n", source->dim, destination->dim, matrix.columns(), matrix.rows());
                matrix.apply(scaling, source, destination);
            }
        }
        /***********
         * Indices *
         ***********/
        static void flipForOrientation(std::vector<bool> const orientation,     //!< Orientation of all triangles; true if correct.
                                       avector *const vector,                   //!< Vector of coefficients to reorient.
                                       size_t const offset = 0                  //!< Only start to flip from the index \p offset.
                                       ) {                                      //! Flip the coefficients of some "Neumann data" to change the orientation between global and local.
            #ifdef USE_OPENMP
            #pragma omp parallel for simd
            #endif
            for(size_t idx = offset; idx < vector->dim; ++idx)
                if(!orientation[idx])
                    (void) scaleentry_avector(vector, idx, -1);
        }

        static inline avector const * new_const_sub_avector(avector const * const vector,       //!< Owner vector.
                                                            uint const dim,                     //!< Dimension of the view.
                                                            uint const off                      //!< Offset to start the view.
                                                            ) {                                 //! Take a view of a constant avector. \return Constant view.
            return const_cast<avector const *>(new_sub_avector(const_cast<avector *>(vector), dim, off));
        }
        static inline void del_const_sub_avector(avector const * vector //!< Consant view to delete.
                                                 ) {                    //! Delete a constant vector, using const_cast. \warning This is intended for views only, not for vectors owning their elements. \todo Check that the vector actually doesn't have ownership.
            del_avector(const_cast<avector *>(vector));
        }
        /*********
         * Other *
         *********/
        static inline cluster* buildGeometricCluster(clustergeometry * const geometry,  //!< Cluster geometry to use for the current cluster.
                                                     std::vector<uint> &indices         //!< Indices to use for all slices.
                                                     ) {                                //! \brief Construct a cluster from options and the given parameters. \return A pointer to the new cluster.
            if(indices.size() == 0)
                return nullptr;
            return buildgeometric_clustergeometry(geometry, 64, static_cast<uint>(clusterResolution()), indices.size(), indices.data());
        }

        /*****************
         * SingleSurface *
         *****************/
        SingleSurface::SingleSurface(Mesh const &mesh, size_t const surfaceIndex, real const stiffness_p, real const density_p, MeshReader *gmshInterface)
            : mesh_(mesh),
          density_(density_p),
          stiffness_(stiffness_p),
          gmshInterface_(gmshInterface) {
            orient_.resize(mesh_.triangles());
            Logger::log(Logger::Color::green, Logger::Level::debug, "Surface %ld: membership vector length: %ld\n", surfaceIndex, mesh_.membership_[surfaceIndex].size());
            assert(orient_.size() + mesh_.vertices() == mesh_.membership_[surfaceIndex].size());
            for(size_t idx = 0; idx < mesh.triangleDofs_; ++idx)
                orient_[idx] = (mesh_.membership_[surfaceIndex][idx + mesh_.vertexDofs_] >= 0);
            for(size_t idx = mesh.triangleDofs_; idx < orient_.size(); ++idx)
                orient_[idx] = (mesh_.membership_[surfaceIndex][idx + mesh_.vertices()] >= 0);
            geometry_ = orientCopy(mesh_.mesh_, orient_);
            // assert(check()); //normally false for multi-domain.

            // Add indices to the right vector (if any).
            at(indices_, SpaceSlices::DirichletSolution).reserve(mesh.vertexDofs_);
            at(indices_, SpaceSlices::DirichletData).reserve(mesh.vertices() - mesh.vertexDofs_);
            at(indices_, SpaceSlices::NeumannSolution).reserve(mesh.triangleDofs_);
            at(indices_, SpaceSlices::NeumannData).reserve(mesh.triangles() - mesh.triangleDofs_);

            #ifdef USE_OPENMP
            #pragma omp parallel sections
            #endif
            {
                #ifdef USE_OPENMP
                #pragma omp section
                #endif
                for(size_t vertex = 0; vertex < mesh.vertexDofs_; ++vertex)
                    if(mesh_.membership_[surfaceIndex][vertex] != 0)
                        at(indices_, SpaceSlices::DirichletSolution).push_back(vertex);
                #ifdef USE_OPENMP
                #pragma omp section
                #endif
                for(size_t vertex = mesh_.vertexDofs_; vertex < mesh.vertices(); ++vertex)
                    if(mesh_.membership_[surfaceIndex][vertex + mesh.triangleDofs_] != 0)
                        at(indices_, SpaceSlices::DirichletData).push_back(vertex);
                #ifdef USE_OPENMP
                #pragma omp section
                #endif
                for(size_t triangle = 0; triangle < mesh.triangleDofs_; ++triangle)
                    if(mesh_.membership_[surfaceIndex][triangle + mesh.vertexDofs_] != 0)
                        at(indices_, SpaceSlices::NeumannSolution).push_back(triangle);
                #ifdef USE_OPENMP
                #pragma omp section
                #endif
                for(size_t triangle = mesh_.triangleDofs_; triangle < mesh_.triangles(); ++triangle)
                    if(mesh_.membership_[surfaceIndex][triangle + mesh_.vertices()] != 0)
                        at(indices_, SpaceSlices::NeumannData).push_back(triangle);
            }
            for(std::vector<uint> & indexVector : indices_) {
                indexVector.shrink_to_fit();
                Logger::log(Logger::Color::green, Logger::Level::debug, "%ld indices.\n", indexVector.size());
            }
            setupCluster();
        }
        void SingleSurface::buildOperator(diradmdata const &admissibilityData) {
            assert(static_cast<size_t>(OperatorSlices::COUNTDof) == 3);
            std::array<helmholtz3d*, static_cast<size_t>(OperatorSlices::COUNTDof)> bemDescriptions {
                new_slp_helmholtz3d(   geometry(), waveNumber(),      static_cast<uint>(quadraturePoints()), static_cast<uint>(interpolationPoints())),
                new_dlp_cl_helmholtz3d(geometry(), waveNumber(), 0, static_cast<uint>(quadraturePoints()), static_cast<uint>(interpolationPoints())),
                new_hyp_ll_helmholtz3d(geometry(), waveNumber(),      static_cast<uint>(quadraturePoints()), static_cast<uint>(interpolationPoints()))
                };
            
            Logger::indent();
            for(size_t slice = 0; slice < static_cast<size_t>(OperatorSlices::COUNTDof); ++slice) {
                std::pair<SpaceSlices, SpaceSlices> clusterIndices = operatorSpaces(slice);
                operators_[slice].build(bemDescriptions[slice], directionalClusters_[clusterIndices.first], directionalClusters_[clusterIndices.second], admissibilityData);
            }
        }
        void SingleSurface::setupCluster() {
            Logger::log(Logger::Color::yellow, Logger::Level::programFlow, "Creating cluster trees (resolution %u). Normal: (descendents, size)", clusterResolution());
            clustergeometry *clusterTriangleGeometry, *clusterVertexGeometry;
                clusterTriangleGeometry = buildgeometry_surface3d(geometry());     // auxiliary data only needed here
                clusterVertexGeometry = buildgeometry_vertex_surface3d(geometry());// auxiliary data only needed here
            // #ifdef USE_OPENMP
            // #pragma omp parallel
            // #pragma omp single nowait
            // #endif
            // {
                // #ifdef USE_OPENMP
                // #pragma omp task depend(in: clusterTriangleGeometry) depend(out: this->clusters_[SpaceSlices::NeumannSolution])
                // #endif
                clusters_[SpaceSlices::NeumannSolution]   = buildGeometricCluster(clusterTriangleGeometry, at(indices_, SpaceSlices::NeumannSolution));
                // #ifdef USE_OPENMP
                // #pragma omp task depend(in: clusterTriangleGeometry) depend(out: this->clusters_[SpaceSlices::NeumannData])
                // #endif
                clusters_[SpaceSlices::NeumannData]       = buildGeometricCluster(clusterTriangleGeometry, at(indices_, SpaceSlices::NeumannData));
                // #ifdef USE_OPENMP
                // #pragma omp task depend(in: clusterVertexGeometry) depend(out: this->clusters_[SpaceSlices::DirichletSolution])
                // #endif
                clusters_[SpaceSlices::DirichletSolution] = buildGeometricCluster(clusterVertexGeometry,   at(indices_, SpaceSlices::DirichletSolution));
                // #ifdef USE_OPENMP
                // #pragma omp task depend(in: clusterVertexGeometry) depend(out: this->clusters_[SpaceSlices::DirichletData])
                // #endif
                clusters_[SpaceSlices::DirichletData]     = buildGeometricCluster(clusterVertexGeometry,   at(indices_, SpaceSlices::DirichletData));
            // }
            #ifdef USE_OPENMP
            #pragma omp parallel
            #pragma omp single nowait
            #endif
            {
                for(size_t slice = 0; slice < size<SpaceSlices>(); ++slice)
                    #ifdef USE_OPENMP
                    #pragma omp task depend(in: this->clusters_[slice])
                    #endif
                    {
                        if(clusters_[slice] == nullptr)
                            directionalClusters_[slice] = nullptr;
                        else
                            directionalClusters_[slice] = buildfromcluster_dcluster(clusters_[slice]);
                    }
                #ifdef USE_OPENMP
                #pragma omp task depend(in: this->clusters_)
                #endif
                CFunctions::freeCPointer(clusterVertexGeometry);
                #ifdef USE_OPENMP
                #pragma omp task depend(in: this->clusters_)
                #endif
                CFunctions::freeCPointer(clusterTriangleGeometry);
                #ifdef USE_OPENMP
                #pragma omp task depend(in: this->clusters_)
                #endif
                for(size_t slice = 0; slice < size<SpaceSlices>(); ++slice)
                    if(clusters_[slice] == nullptr)
                        Logger::log(Logger::Color::green, Logger::Level::programFlow, "(empty) ");
                    else
                        Logger::log(Logger::Color::green, Logger::Level::programFlow, "%u, %u; ", clusters_[slice]->desc, clusters_[slice]->size);
            }
            Logger::log(Logger::Color::green, Logger::Level::programFlow, ". Done\n");
        }
        void SingleSurface::apply(field const scaling, avector const * const source, avector * const destination) const {
            // take views
            avector const *sourceD = new_const_sub_avector(source, vertexDofs(), 0), *sourceTrueN = new_const_sub_avector(source, triangleDofs(), vertexDofs());
            avector *destinationD = new_sub_avector(destination, vertexDofs(), 0), *destinationN = new_sub_avector(destination, triangleDofs(), vertexDofs());
            avector *sourceN = new_avector(sourceTrueN->dim);
            copy_avector(sourceTrueN, sourceN);
            flipForOrientation(orient_, sourceN);
            flipForOrientation(orient_, destinationN);     // Do it as well so that the result of what is there already is intact.

            // Apply operators
            constAt(operators_, OperatorSlices::SolutionSolutionSingle).apply(scaling * stabilizer(), sourceN, destinationN);
            constAt(operators_, OperatorSlices::SolutionSolutionDouble).apply(-scaling, sourceD, destinationN);
            constAt(operators_, OperatorSlices::SolutionSolutionDouble).applyTransposed(scaling, sourceN, destinationD);
            constAt(operators_, OperatorSlices::SolutionSolutionHypersingular).apply(scaling / stabilizer(), sourceD, destinationD);
            flipForOrientation(orient_, destinationN);

            del_avector(sourceN);
            // delete views
            del_avector(destinationD);
            del_avector(destinationN);
            del_const_sub_avector(sourceD);
            del_const_sub_avector(sourceTrueN);
        }
        field SingleSurface::waveNumber() const {
            return waveNumber_;
        }
        field SingleSurface::waveNumber(field const waveNumber_p) {
            waveNumber_ = waveNumber_p;
            stabilizer_ = - waveNumber_ * _Complex_I * sqrt(stiffness_/density_);
            return waveNumber();
        }
        field SingleSurface::stabilizer() const {
            return stabilizer_;
        }
        field SingleSurface::stabilizer(field const stabilizer_p) {
            stabilizer_ = stabilizer_p;
            waveNumber_ = stabilizer_ * _Complex_I * sqrt(density_/stiffness_);
            return stabilizer();
        }
        void SingleSurface::makeOffset(Functions::MeshOffset const &solution, helmholtz3d const *const parameters) {
            (void)stabilizer(parameters->kappa);      // This sets s and κ.
            offset_ = new_avector(fullBasisSize());
            clear_avector(*offset_);
            solution.c0Dirichlet(mesh_, at(indices_, SpaceSlices::DirichletSolution), true, sqrt(density_/stiffness_), *offset_);
            solution.c0Dirichlet(mesh_, at(indices_, SpaceSlices::DirichletData),    false, sqrt(density_/stiffness_), *offset_);
            solution.l2Neumann(parameters, mesh_, at(indices_, SpaceSlices::NeumannSolution), true, sqrt(density_/stiffness_), stiffness_, *offset_);
            solution.l2Neumann(parameters, mesh_, at(indices_, SpaceSlices::NeumannData),    false, sqrt(density_/stiffness_), stiffness_, *offset_);
            gmshInterface_->saveView(MeshReader::View::offset, *offset_);
        }
        void SingleSurface::makeMultiRhs(diradmdata const &directionalAdmissibilityData, avector *const rhs) {
            Logger::log(Logger::Color::cyan, Logger::Level::dataOperation, "Computing disposable part of Calderón and applying to offset for rhs... ");
            std::array<helmholtz3d*, 4> bemDescriptions = {
                new_slp_helmholtz3d(   geometry(), waveNumber(),      static_cast<uint>(quadraturePoints()), static_cast<uint>(interpolationPoints())),
                new_dlp_cl_helmholtz3d(geometry(), waveNumber(), 0.5, static_cast<uint>(quadraturePoints()), static_cast<uint>(interpolationPoints())),
                new_dlp_cl_helmholtz3d(geometry(), waveNumber(),-0.5, static_cast<uint>(quadraturePoints()), static_cast<uint>(interpolationPoints())),
                new_hyp_ll_helmholtz3d(geometry(), waveNumber(),      static_cast<uint>(quadraturePoints()), static_cast<uint>(interpolationPoints()))
            };
            field scaling[4] {-stabilizer(), 1, -1, -1/stabilizer()};
            avector *source[2] {        // Takes the correct view of offset_ (source)
                new_sub_avector(*offset_, triangles(), vertices()),     // The vector needs to have full vertex size for the operator to work correctly. Shift only by "extraneous" vertices.
                new_sub_avector(*offset_, vertices(), triangleDofs())   // The vector needs to have full vertex size for the operator to work correctly. Shift only by "extraneous" triangleDofs, not vertexDofs.
            };
            avector *destination[2] {   // Takes the correct view of rhs (destination)
                new_sub_avector(rhs, triangleDofs(), vertexDofs()),     // For the destination only the dof part is used (since rhs is in the homogeneous single trace space)
                new_sub_avector(rhs, vertexDofs(), 0)                   // For the destination only the dof part is used (since rhs is in the homogeneous single trace space)
            };

            Logger::log(Logger::Color::cyan, Logger::Level::debug, false, "Flipping... ");
            flipForOrientation(orient_, destination[0]);
            flipForOrientation(orient_, source[0], triangleDofs());

            Logger::log(Logger::Color::cyan, Logger::Level::debug, false, "Applying... ");
            fflush(stdout);
            for(size_t idx = 0; idx < 4; ++idx) {
                makeTmpApply(scaling[idx],                      // Correct scaling.
                             bemDescriptions[idx],              // Auxiliary data from h2lib.
                             directionalClusters_[2*(1 - idx/2)],// NeumannSolution twice, then DirictletSolution twice.
                             directionalClusters_[3 - 2*(idx%2)],// Alternate NeumannData and DirictletData.
                             directionalAdmissibilityData,      // Misc parameters.
                             source[idx%2],                     // Input.
                             destination[idx/2],                // Output.
                             idx == 2                           // The only one to transpose.
                             );
            }
            Logger::log(Logger::Color::cyan, Logger::Level::debug, false, "Flipping... ");
            flipForOrientation(orient_, destination[0]);
            flipForOrientation(orient_, source[0], triangleDofs());

            for(size_t idx = 0; idx < 2; ++idx) {
                del_avector(source[idx]);
                del_avector(destination[idx]);
            }
            Logger::log(Logger::Color::cyan, Logger::Level::debug, true, "Done.");
        }
        void SingleSurface::makeSingleRhs(avector *const rhs) const {
            Logger::log(Logger::Color::cyan, Logger::Level::dataOperation, "Applying regular part of Calderón to offset for rhs...");

            apply(-1, *offset_, rhs);      // the rhs gat a - sign
        }

        std::pair<real, real> SingleSurface::computeDifferenceNorm(avector *const solution, Functions::MeshOffset const &reference, helmholtz3d const *const parameters) const {
            Logger::log(Logger::Color::green, Logger::Level::programFlow, true, "Computing and plotting comlpete solution on subdomain.");
            avector *completeSolution = clone_avector(*offset_);
            avector *completeSolutionRestrict = new_sub_avector(completeSolution, solution->dim, 0);
            add_avector(1, solution, completeSolutionRestrict);
            gmshInterface_->saveView(MeshReader::View::solution, completeSolution);
            del_avector(completeSolution);
            del_avector(completeSolutionRestrict);

            Logger::log(Logger::Color::green, Logger::Level::programFlow, true, "Computing and plotting restricted reference.");
            DumbPointer<avector> discreteReference = new_avector(solution->dim);
            clear_avector(*discreteReference);
            reference.c0Dirichlet(mesh_, at(indices_, SpaceSlices::DirichletSolution), false, sqrt(density_/stiffness_), *discreteReference);
            reference.l2Neumann(parameters, mesh_, at(indices_, SpaceSlices::NeumannSolution), false, sqrt(density_/stiffness_), stiffness_, *discreteReference);

            Logger::log(Logger::Color::green, Logger::Level::debug, false, "Subtracting offset to reference to make it single-trace... ");
            avector const *offsetRelevant = new_const_sub_avector(*offset_, solution->dim, 0);
            add_avector(-1, offsetRelevant, *discreteReference);      // Get to the single trace function
            del_const_sub_avector(offsetRelevant);
            Logger::log(Logger::Color::green, Logger::Level::debug, false, "Computing reference norm... ");
            double solutionNorm = energyNorm(*discreteReference, this, stabilizer());
            Logger::log(Logger::Color::green, Logger::Level::debug, true, "%lf.", solutionNorm);
            gmshInterface_->saveView(MeshReader::View::solution, *discreteReference);

            Logger::log(Logger::Color::green, Logger::Level::debug, false, "Subtracting solution to reference to compute error... ");
            #ifdef USE_OPENMP
            #pragma omp parallel for simd
            #endif
            for(size_t idx = 0; idx < at(indices_, SpaceSlices::DirichletSolution).size(); ++idx) {
                addentry_avector(*discreteReference, at(indices_, SpaceSlices::DirichletSolution)[idx], -getentry_avector(solution, at(indices_, SpaceSlices::DirichletSolution)[idx]));
            }
            #ifdef USE_OPENMP
            #pragma omp parallel for simd
            #endif
            for(size_t idx = 0; idx < at(indices_, SpaceSlices::NeumannSolution).size(); ++idx) {
                addentry_avector(*discreteReference, at(indices_, SpaceSlices::NeumannSolution)[idx]+vertexDofs(), -getentry_avector(solution, at(indices_, SpaceSlices::NeumannSolution)[idx]+ vertexDofs()));
            }
            // add_avector(-1, solution, *discreteReference);      // Now discreteReference is the dof-wise error
            gmshInterface_->saveView(MeshReader::View::error, *discreteReference);
            Logger::log(Logger::Color::green, Logger::Level::debug, true, "Computing and returning error norm.");
            return std::make_pair(solutionNorm, energyNorm(*discreteReference, this, stabilizer()));
        }

        size_t SingleSurface::triangles() const {
            return mesh_.triangles();
        }
        size_t SingleSurface::vertices() const {
            return mesh_.vertices();
        }
        size_t SingleSurface::neumannBasisSize() const {
             return triangles(); 
        }
        size_t SingleSurface::dirichletBasisSize() const {
             return vertices(); 
        }
        size_t SingleSurface::fullBasisSize() const {
             return dirichletBasisSize() + neumannBasisSize(); 
        }
        size_t SingleSurface::singleSpaceSize() const {
             return dirichletBasisSize() + neumannBasisSize();
        }

        size_t SingleSurface::vertexDofs() const {
            return mesh_.vertexDofs_;
        }
        size_t SingleSurface::triangleDofs() const {
            return mesh_.triangleDofs_;
        }
        surface3d const * SingleSurface::geometry() const {
            return geometry_;
        }

        bool SingleSurface::check() const {
            bool gr_closed = isclosed_surface3d(geometry());
            bool gr_oriented = isoriented_surface3d(geometry());
            if(gr_closed && gr_oriented) {
                Logger::log(Logger::Color::green, Logger::Level::developer, true, "Surface is closed and oriented");
                return true;
            }
            Logger::log(Logger::Color::red, Logger::Level::warning, true, "Surface is:%s closed,%s oriented", (gr_closed ? "" : " NOT"), (gr_oriented ? "" : " NOT"));
            return false;
        }
        void SingleSurface::properties() const {
            double hmin, hmax, anglemin, angleedge;
            getproperties_surface3d(geometry(), &hmin, &hmax, &anglemin, &angleedge);
            Logger::log("  %u vertices, %u edges, %u triangles\n  hmin %.2e, hmax %.2e, anglemin %g, angleedge %g\n",
                        geometry_->vertices, geometry_->edges, geometry_->triangles, hmin, hmax, anglemin, angleedge);
        }
    }
}
