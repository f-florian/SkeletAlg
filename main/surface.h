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
#ifndef SURFACE_H
#define SURFACE_H

#include <complex.h>
#include <vector>

struct _dcluster;
typedef struct _dcluster dcluster;
struct _diradmdata;
typedef struct _diradmdata diradmdata;
struct _dh2matrix;
typedef struct _dh2matrix dh2matrix;
struct _truncmode;
typedef struct _truncmode truncmode;

#include "h2lib/surface3d.h"

#include "utils/arrays.h"
#include "utils/cutils.h"
#include "utils/enum.h"

#include "dh2operator.h"
#include "functions.h"
#include "spaces.h"

namespace MainNamespace {
    struct AdmissibilityParameters;
    class MeshReader;
    namespace Surface {
        // control whether dirichlet or neumann is first in the combined vector
        constexpr inline uint dirichletOffset(size_t const dim [[maybe_unused]]){return 0;}
        constexpr inline uint neumannOffset(const size_t dim){return static_cast<uint>(dim);}

        /**
         * Manage a closed surface, with its associated operators.
         * \todo expand.
         */
        class SingleSurface {
        public:
            SingleSurface(Mesh const &mesh,             //!< Global mesh.
                          size_t const surfaceIndex,    //!< Indices of the current mesh.
                          real const stiffness_p,       //!< Domain stiffness parameter.
                          real const density_p,         //!< Domain density parameter.
                          MeshReader *gmshInterface     //!< Interface to gmsh, e.g., to print offset or solution.
                          );                            //!< \brief Construct an object to manage data associated to a closed boundary.
            SingleSurface(const SingleSurface &other) = delete; //!< Copy constructor would break a few things. But a copy should never be needed anyway.
            SingleSurface(SingleSurface &&other) = default;     //!< Move constructor.
            void buildOperator(diradmdata const &admissibilityData      //!< Parameters for constructing clusters.
                               );                                       //!< Build the operators.
            void makeOffset(Functions::MeshOffset const &solution,      //!< Offset function.
                            helmholtz3d const *const parameters         //!< Paratemers, like quadrature.
                            );                                          //!< Project the offset on the current mesh.
            void makeMultiRhs(diradmdata const &admissibilityData,      //!< Admissibility data
                              avector *const rhs                        //!<[in, out] Rhs vector to add to. It must be allocated and may contain data.
                              );                                        //!< Project the offset on the current mesh, then use it to compute the right-hand side. \todo Make const (fix orientation flip).
            void makeSingleRhs(avector *const rhs                       //!<[in, out] Rhs vector to add to. It must be allocated and may contain data.
                               ) const;                                 //!< Project the offset on the current mesh, then use it to compute the right-hand side.
            /**
             * \brief Apply the (scaled) Calderón operator (previously constructed) to the source vector, and add it to destination.
             *
             * If \f$A\f$ is the current discretization of the Calderón operator, compute \f[
             * y \leftarrow y + \alpha A x
             * \f]
             */
            void apply(field const scaling,             //!< \f$\alpha\f$.
                       avector const * const source,    //!< \f$x\f$.
                       avector * const destination      //!< [in,out] \f$y\f$.
                       ) const;                         //! \brief Apply the Calderón operator to the "solution" part of \p source.
            void properties() const;    //!< Log the value of some properties in a message. \brief Compute and report mesh properties.
            void scaleNeumann(avector *const traces     //!<[in,out] The full vector of traces.
                              ) const;                  //!< Scale the portion of the vector corresponding to the Neumann trace using the stiffness parameter.
            void unscaleNeumann(avector *const traces   //!<[in,out] The full vector of traces.
                                ) const;                //!< Undo the scaling of the portion of the vector corresponding to the Neumann trace using the stiffness parameter.
            field waveNumber() const;                   //!< Getter. \return The (Helmholtz-type) wave number κ, as in κ²u + Δu = 0.
            field waveNumber(field const waveNumber_p   //!< New wave number κ.
                             );                         //!< Setter for waveNumber_ and stabilizer_. \return The new value of κ.
            field stabilizer() const;                   //!< Getter. \return The (wave-type) wave number s, as in ps²u - AΔu = 0.
            field stabilizer(field const stabilizer_p   //!< New stabilization scaling s.
                             );                         //!< Setter for waveNumber_ and stabilizer_. \return The new value of s.
            std::pair<real,real> computeDifferenceNorm(avector *const solution,                 //!< The approximate solution vector.
                                                       Functions::MeshOffset const &reference,  //!< The reference solution.
                                                       helmholtz3d const *const parameters      //!< Object with quadrature parameters.
                                                       ) const;                                 //!< Compute the energy error norm. \return The error norm of the difference between the \p reference solution and the approximate \p solution.
        private:
            size_t singleSpaceSize() const;     //!< Getter. \return the size of the single trace space, including all degrees of freedom on the boundary.
            size_t triangles() const;           //!< Getter. \return the number of triangles in the full mesh.
            size_t vertices() const;            //!< Getter. \return the number of vertices in the full mesh.
            size_t neumannBasisSize() const;    //!< Getter. \return ?
            size_t dirichletBasisSize() const;  //!< Getter. \return ?
            size_t fullBasisSize() const;       //!< Getter. \deprecated
            size_t vertexDofs() const;          //!< Getter. \return The number of vertices that are degrees of freedom in the current mesh.
            size_t triangleDofs() const;        //!< Getter. \return The number of triangles that are degrees of freedom in the current mesh.
            surface3d const * geometry() const; //!< Getter. \return A pointer to the H2Lib mesh structure.
            
            inline void cleanData() {   //!< Clear all operators before constructing new ones (e.g, for a different wavenumber).
                clusters_.clear();
                directionalClusters_.clear();
            }
            bool check() const; //!< \brief Check that the surface is a boundary. \return true if the geometry is a closed, oriented manifold.

            void setupCluster();        //!< \brief Setup cluster trees. \todo better description

            Mesh const& mesh_;          //!< Reference to the whole mesh.
            surface3d * geometry_;      //!< Reoriented surface3d object. \todo delete. \todo solve orientation issues first.
            std::vector<bool> orient_;  //!< Indicates where orientation is the same in parent and current surface.

            IndexSpace<uint> indices_ = {};                     //!< Indices used in the current surface.
            SlicedSpace<cluster> clusters_ = {};                //!< Roots of nondirectional cluster trees.
            SlicedSpace<dcluster> directionalClusters_ = {};    //!< Roots of directional cluster trees.
            Operators<DH2MatrixWrapper> operators_ = {};        //!< Boundary operator, sliced to minimize the part to fill.
            DumbPointer<avector> offset_ = nullptr;             //!< Offset (which is multi trace hence belongs here).
            field waveNumber_;  //!< Wave number κ, as in the Helmholtz-type equation Δu+κ²u=0.
            field stabilizer_;  //!< Stabilization scaling for the single layer and hypersingular operators: s:=-iκ as in ps²u - AΔu=0.
            real density_;      //!< Density (mass per volume) of the enclosed domain.
            real stiffness_;    //!< Stiffness of the enclosed domain. \todo Should be a matrix.
            MeshReader *gmshInterface_ = nullptr;       //!< Interface to gmsh, mainly to write solutions.
        };
    }
}

#endif /* SURFACE_H */
