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
#ifndef MULTISURFACE_H
#define MULTISURFACE_H

#include <complex.h>

#include <vector>

#include "h2lib/helmholtz3d.h"
#include "h2lib/surface3d.h"

#include "logger/logger.h"
#include "utils/cutils.h"
#include "functions.h"
#include "surface.h"
#include "meshreader.h"

namespace MainNamespace {
    struct AdmissibilityParameters;
    namespace Surface {
        // static constexpr size_t errorsCount = 7;    //!< Compute the error in this many ways.
        using ErrorVector = real;      //!< A list of errors, computed in different ways.
        using Result = std::pair<size_t, ErrorVector>;          //!< The numerically relevant result: number of GMRES iterations, and list of errors.

        /**
         * \brief Representation of a boundary.
         *
         * Holds the boundary mesh and the associated mass matrices, as well as a list of SingleSurface "children",
         * i.e., subsets of the boundary mesh which are the border of subdomains, and therefore closed manifolds.
         *
         * Through its children it holds (a discretization of) the calderon operator \b C and therefore the calderon bilinear form c.
         * This object and its children hold all the information needed to solve a boundary problem in variational form "Find y such that c(x,y) = rhs_(x) for all x in the chosen basis".
         * It also knows which basis elements are unkonwns and which are given boundary data.
         * 
         * When solving the problem boundary data are treated with an offset function (offset_).
         * This generate a right hand side (rhs_) and leaves homogeneous boundary conditions.
         */
        class MultiSurface {
        public:
            using Children = std::vector<SingleSurface>;
            MultiSurface(std::string const meshFileName //!< Filename of the mesh to load.
                         );                             //!< Contructor. Parse a mesh file in gmsh format and build the global mesh.
            ~MultiSurface();                            //!< \brief Destructor.

            Result solve(AdmissibilityParameters const admissibilityParameters, //!< Admissibility parameters  to build the operators.
                         Functions::MeshOffset &offset                          //!< Offset function, and information about the reference solution.
                         );                                                     //!< Compute the right-hand side, the operators, then solve the linear system. \return The number of iterations, and the error with respect to the reference, computed in different ways.
            
            std::pair<double, double> properties() const;               //!< \brief Print surface properties. \return Mesh size (widest element) and ratio between minimum and maximum size.
            inline size_t triangles() const { return static_cast<size_t>(geometry()->triangles); }      //!< \brief Getter. \return Number of triangles in the mesh.
            inline size_t vertices() const { return static_cast<size_t>(geometry()->vertices); }        //!< \brief Getter. \return Number of vertices in the mesh.
            inline size_t dirichletBasisSize() const { return mesh().vertexDofs_; }                      //!< \brief Getter. \return The size of the Dirichlet trial space.
            inline size_t neumannBasisSize() const { return mesh().triangleDofs_; }                      //!< \brief Getter. \return The size of the Neumann trial space.
            inline size_t problemSize() const { return dirichletBasisSize() + neumannBasisSize(); }     //!< \brief Get the dimension of the space of Cauchy space (i.e., Dirichlet and Neumann together). \return The size of the Cauchy trial space.
        private:
            inline field waveNumber() const { return metric_->kappa; }                                  //!< \brief Getter. \return The wave number.
            inline size_t quadratureNodes() const { return metric_->nq; }                               //!< \brief Getter. \return The number of quadrature nodes. \todo per dimensor or total?.
            inline size_t interpolationNodes() const { return metric_->ni; }                            //!< \brief Getter. \return The number of interpolation nodes. \todo per dimensor or total?.

            /**
             * Fill destination with the elements of source or 0.
             *
             * Fill the elements of destination corresponding to a base element with the appropriate element in source, set the others to 0.
             * Use nonZero_ to determine if an index corresponds to a base element.
             * @param source vector with the components not corresponding to a basis element removed.
             * @param destination vector with the components not corresponding to a basis element set to 0.
             */
            inline void fillZeros(avector const * const source, avector * const destination) const {
                assert(destination->dim == triangles() + vertices());
                assert(source->dim == problemSize());
                clear_avector(destination);
                for(size_t idxDestination = 0; idxDestination < problemSize(); ++idxDestination)
                    setentry_avector(destination, idxDestination, getentry_avector(source, idxDestination));
            }
            inline void filterBasis(avector const * const source,       //!< Source vector of length \ref triangles() + \ref vertices().
                                    avector * const destination         //!<[out] Filtered vector of length \ref problemSize(). Must be already allocated.
                                    ) const {                           //! Copy the rest to \p destination. \brief Remove the components of \p source not part of the basis of the trial space.
                assert(source->dim == triangles() + vertices());
                assert(destination->dim == problemSize());
                for(size_t idxDestination = 0; idxDestination < problemSize(); ++idxDestination)
                    setentry_avector(destination, idxDestination, getentry_avector(source, idxDestination));
            }
            void apply(field const scaling,             //!< Multiplicative factor to apply.
                       avector const * const source,    //!< Input for the operator.
                       avector * const destination      //!<[out] \p destination + \p scaling * \b C (\p source).
                       ) const;                         //!< Apply the Calderon operator \b C to source, but only on components for the trial space.
            ErrorVector error(avector *const solution,                  //!< Discrete solution vector.
                              Functions::MeshOffset const &reference    //!< Exact solution (reference).
                              ) const;                                  //!< \todo fix doc. The (D) and (N) versions are computed also including the components that should be 0, the others only on the unknown part. \brief Compute the error in several ways. \return The following error norms (in order): energy; max; max(D); max(N); l2; l2(D); l2(N). \todo relative error.

            // private getters and setters
            Mesh& mesh() {return gmshInterface_.mesh();}
            Mesh const & mesh() const {return gmshInterface_.mesh();}
            void waveNumber(field const kappa) {
                if(metric_ != nullptr)
                    metric_->kappa = kappa;
            }
            surface3d* geometry() {
                return mesh().mesh_;
            }
            surface3d const * geometry() const {
                return mesh().mesh_;
            }
   
            Children children_; //!< List of subdomains.
            mutable MeshReader gmshInterface_;  //!< Interface to gmsh.
            helmholtz3d *metric_ = nullptr;     //!< Data structure needed to generate metric data (e.g. mass matrix).
        };
    }
}

#endif /* MULTISURFACE_H */
