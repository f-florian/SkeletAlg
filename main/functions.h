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
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "h2lib/settings.h"

#include "utils/cutils.h"

#include "meshreader.h"

#include <array>
#include <vector>

namespace MainNamespace {
    template <class spaceFunction>
    constexpr field (*dirichletPointer())(pcreal, pcreal, void *) {
        return [](pcreal x, [[maybe_unused]] pcreal n, void *data) -> field {
            return static_cast<spaceFunction*>(data)->dirichlet({x[0], x[1], x[2]});
        };
    }
    template <class spaceFunction>
    constexpr field (*neumannPointer())(pcreal, pcreal, void *) {
        return [](pcreal x, pcreal n, void *data) {
            return static_cast<spaceFunction*>(data)->neumann({x[0], x[1], x[2]}, {n[0], n[1], n[2]});
        };
    }

    namespace Functions {
        //! \brief Virtual class that specifies the reference traces for testing. \todo spacetime?
       class Reference {
        public:
            virtual ~Reference(){}
            virtual field dirichlet(SpaceVector<real const> const point,//!< Point where the trace is computed.
                                    real const scale                    //!< Scaling for the frequency, due to material parameters.
                                    ) const = 0;                        //!< Compute the (discrete) (Dirichlet) trace of the function at \p point.
            virtual field neumann(SpaceVector<real const> const point,  //!< Point where the trace is computed.
                                  SpaceVector<real const> const normal, //!< Normal vector at \p point.
                                  real const scale                      //!< Scaling for the frequency, due to material parameters.
                                  ) const = 0;                          //!< Compute the (discrete) Neumann trace of the function at \p point.
           virtual field waveNumber() const = 0;       //!< Getter. \return The (modulus of the) wave number.
        };
  
        //! Plane wave reference.
        class PlaneWave
            :public Reference {
        public:
            PlaneWave() = delete;
            PlaneWave(SpaceVector<double const> const frequency,//!< Wave direction (will be normalized).
                      field const waveNumber                    //!< Wave number (modulus).
                      );                                        //!< Construct a reference plane wave.
            virtual field dirichlet(SpaceVector<real const> const point, real const scale) const override;
            virtual field neumann(SpaceVector<real const> const point, SpaceVector<real const> const normal, real const scale) const override;
            virtual field waveNumber() const override {return waveNumber_;}
        private:
            SpaceVector<field> frequency_;      //!< Normalized wave direction. Its norm is waveNumber_.
            field waveNumber_;
        };
  
        //! \brief Offset that takes its values from a mesh data.
        class MeshOffset
            : public Reference {
        public:
            MeshOffset(Reference const *const solution  //!< The exact solution, needed at least to compute the traces on the non-dof elements.
                       );                               //!< Constuctor. Set the trace scaling and the actual exact solution. This class is pure virtual, but these data have to be initialized for every derived class.
            virtual void c0Dirichlet(Mesh const &mesh,                  //!< Mesh data.
                                     std::vector<uint> const &indices,  //!< Indices of the dofs to fill.
                                     bool const dof,                    //!< Whether the indices correspond to degrees of freedom or data.
                                     real const waveScale,              //!< Scaling of the wave number induced by parameters.
                                     avector *const coefficients        //!<[out] The offset vector to fill.
                                     ) const;                           //!< Interpolate the Dirichlet trace with continuous, picewise linear functions. The trace is scaled by s^(1/2).
            virtual void l2Neumann(helmholtz3d const *const parameters, //!< Helmholtz object with the quadrature parameters.
                                   Mesh const &mesh,                    //!< Mesh description.
                                   std::vector<uint> const &indices,    //!< Indices of the dofs to fill.
                                   bool const dof,                      //!< Whether the indices correspond to degrees of freedom or data.
                                   real const waveScale,                //!< Scaling of the wave number induced by parameters.
                                   real const stiffness,                //!< Stiffness (for scaling the Neumann trace with respect to the normal trace).
                                   avector *const coefficients          //!<[out] The offset vector to fill.
                                   ) const;                             //!< Compute the L² projection of the Nemuann trace into piecewise constant functions. The trace is scaled by s^(-1/2).
            virtual field dirichlet(bool const dof,                     //!< Whether the current point is a dof. The offset on dofs may differ from the trace of the actual solution.
                                    SpaceVector<real const> const point,//!< Point where to compute the offset.
                                    real const scale                    //!< Scaling for the frequency, due to material parameters.
                                    ) const = 0;                        //!< Compute an offset for the Dirictlet trace of the current solution.
            virtual field neumann(bool const dof,                       //!< Whether the current point is a dof. The offset on dofs may differ from the trace of the actual solution.
                                  SpaceVector<real const> const point,  //!< Point where to compute the offset.
                                  SpaceVector<real const> const normal, //!< Normal vector at the given point.
                                  real const scale                      //!< Scaling for the frequency, due to material parameters.
                                  ) const = 0;                          //!< Compute an offset for the Neumann trace of the current solution.
            virtual field dirichlet(SpaceVector<real const> const point, real const scale) const override;
            virtual field neumann(SpaceVector<real const> const point, SpaceVector<real const> const normal, real const scale) const override;
            virtual field waveNumber() const override {return solution_->waveNumber();}
        protected:
            Reference const *solution_; //!< Pointer to the actual solution generating the offset.
            field traceScaling_;        //!< sqrt(-i wavenumber) for trace scaling.
        };

        //! This only works decently for some solutions
        class ZeroExtension
            : public MeshOffset {
        public:
            ZeroExtension(Reference const *const solution       //!< The exact solution, needed to compute the traces on the non-dof elements.
                          ) : MeshOffset(solution) {}           //!< Constuct a zero extension of a given function offset.
            virtual field dirichlet(bool const dof, SpaceVector<real const> const point, real const scale) const override;
            virtual field neumann(bool const dof, SpaceVector<real const> const point, SpaceVector<real const> const normal, real const scale) const override;
        };
    }
}

#endif /* FUNCTIONS_H */
