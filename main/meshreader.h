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
#ifndef MESHREADER_H
#define MESHREADER_H

#include <array>
#include <map>
#include <vector>
#include <string>

#include "utils/enum.h"
#include "linearalgebra.h"

struct _surface3d;
typedef struct _surface3d surface3d;
struct _avector;
typedef struct _avector avector;

namespace MainNamespace {
    struct Mesh {
        using Membership = short int;                                   //!< Data type for membership information. Currently int, as it can indicate orientation.
        using MembershipVector = std::vector<Membership>;               //!< Vector of membership information.
        using SurfaceList = std::map<std::string, MembershipVector>;    //!< List of surfaces (as pairs (name, membership)).
        using SurfaceMembership = std::vector<MembershipVector>;        //!< Vector of surfaces, each constisting of membership (and orientation) information.
        Mesh();
        ~Mesh();
        Mesh(Mesh const &) = delete;
        Mesh(Mesh &&);
        Mesh& operator=(Mesh &&);
        void writeMesh();       //!< Write a mesh in a custom format, at a custom location. Format: one line with number of vertices, triangles, vertexDofs, triangleDofs; then vertices (coordinates), then triangles (vertex indices triplets).
        SpaceVector<double const> vertex(size_t const vertexIdx //!< Vertex index.
                                         ) const;               //!< Get a vertex. \return The coordinates as an array of double.
        SpaceVector<double const> normal(size_t const triangleIdx       //!< Triangle index.
                                         ) const;                       //!< Get the normal to a triangle. \return The coefficients in the standard basis as an array of double.
        std::array<size_t,3> triangle(size_t const triangleIdx  //!< Triangle index.
                                      ) const;                  //!< Get a triangle. \return The indices of the vertices as an array of size_t.

        surface3d* mesh_ = nullptr;     //!< H2Lib mesh data structure.
        SurfaceMembership membership_;  //!< Vectors of surfaces with dof membership.
        size_t vertexDofs_ = 0;         //!< Number of vertices which are actually degrees of freedom and not known data.
        size_t triangleDofs_ = 0;       //!< Number of triangles which are actually degrees of freedom and not known data.
        size_t triangles() const;       //!< Getter. \return The number of triangles.
        size_t vertices() const;        //!< Getter. \return The number of vertices.
    };

    class MeshReader {
    public:
        enum class View : size_t {
            solution,
            rhs,
            offset,
            error,
            COUNT
        };
        MeshReader();
        ~MeshReader();
        void parseMesh(std::string const filename       //!< Name of the meshfile to parse.
                       );                               //!< Parse a gmsh mesh file into a surface3d structure. The mesh_ object is filled in the appropriate way.
        void saveView(View const type,          //!< Name for the view.
                      avector const *const data //!< Data to represent.
                      );                        //!< Save the Cauchy data associated to the mesh, with the given name.
        Mesh &mesh() {return mesh_;}
        Mesh const &mesh() const {return mesh_;}
    private:
        void toGmshIndices();   //!< Increase all indices by 1 to make them 1-based, as per GMSH conventions.
        void fromGmshIndices(); //!< Decrease all indices by 1 to make them 0-based.
        void readData();
        std::string fileName_ = "";     //!< File name for the mesh.
        std::string modelName_ = "";    //!< Model name (as seen by gmsh).
        Mesh mesh_;    //!< Mesh object.
        std::vector<size_t> reordering_;//!< Reordering of the triangles and vertices tags (reverse map: index j stores the position of the element with tag j+1). It belongs here since it is only useful when interacting with gmsh.
        std::vector<size_t> nodeTag_;   //!< Reordering of the vertex (node) tags (forward map: index j stores the tag of the element in that position). It belongs here since it is only useful when interacting with gmsh.
        std::vector<size_t> elementTag_;//!< Reordering of the triangle () tags (forward map: index j stores the tag of the element in that position). It belongs here since it is only useful when interacting with gmsh.
        std::array<int,2*size<View>()> viewTags_;//!< Mapping view tags from type to tags.
        std::array<bool,size<View>()> enabled_;//!< Enabled views.
        std::array<size_t, size<View>()> step_ = {0,0,0,0};
    };
}

#endif /* MESHREADER_H */
