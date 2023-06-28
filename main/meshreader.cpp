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
#include "meshreader.h"

#include <array>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <list>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gmsh.h>

#include "h2lib/surface3d.h"

#include "logger/logger.h"
#include "utils/cutils.h"
#include "options.h"
#include "utils/enum.h"

namespace MainNamespace {
    using SurfaceList = std::map<std::string, Mesh::MembershipVector>;    //!< List of surfaces (as pairs (name, membership)).
    static std::string const dirichletName = "d";       //!< Dirichlet boundary conditions identifier.
    static std::string const neumannName = "n";         //!< Neumann boundary conditions identifier.
    static std::string const names[size<MeshReader::View>()] = {"solution", "rhs", "offset","error"};
    
    enum Type {         //! Ordering of lines, triangle, nodes in an array containing definitions for more than one dimension.
        edgesIdx = 0,
        trianglesIdx,
        nodesIdx
    };

    static inline void parsePhysicalName(std::string const name,//!< Name to process.
                                         int const tag,         //!< Tag of the physical group.
                                         int const dimension,   //!< Dimension of the physical group.
                                         SurfaceList& surfaces  //!< [out] For each surface, for every node, 0 if the node does not belong to the surface. The sign specify the orientation.
                                         ) {
        /** \brief Parse which nodes belong to which surface.
         *
         * Parse a string into a list of strings, separated by any number of whitespaces.
         * These are interpreted as the names of the surfaces which share the nodes in the physical group, or as boundary conditions in the given point.
         *
         * For all surface names (or boundary condition):
         * * if this surface was not parsed before, add it to the \p surfaces structure, as a vector of zeros of the same length as the boundary conditions (always considered as initialized with the right length).
         * * add the nodes to the surface, by setting to 1 (correct orientation) or -1 (reverse orientation) the element with index (nodetag-1), for all nodetag that GMSH considers in the physical group.
         *
         * Naming conventions:
         * * Dirichlet boundary: "d";
         * * Neumann boundary: "n";
         * * Names are case-insensitive, instead any capital letter in the name means that the group has the wrong orientation.
         */

        size_t start, end = 0;          // Start and end point of the current word.
        std::vector<size_t> nodesGroup; // Tags of the nodes in the group.
        std::vector<double> ignore;     // Needed as an input for gmsh but discarded
        gmsh::model::mesh::getNodesForPhysicalGroup(dimension, tag, nodesGroup, ignore); ignore.clear();

        Logger::log(Logger::Color::yellow, Logger::Level::debug, true, "Parsing physical name %s", name.c_str());
        for(;;) {
            start = end;
            while((start < name.size()) && isspace(name[start]))
                ++start;
            if(start == name.size())
                break;
            end = start;        // Now start is the next non-whitespace. Let's find the end of the word.
            while((end < name.size()) && !isspace(name[end]))
                ++end;
            std::string domain = name.substr(start, end-start);
            std::string lower = domain;
            for(size_t idx = 0; idx < lower.length(); ++idx) 
                lower[idx] = tolower(lower[idx]);
            int short orientation = (domain==lower) ? 1 : -1;
            Logger::log(Logger::Color::cyan, Logger::Level::ioBasic, true, "Physical group %s (dim %d, tag %d): found domain %s (%s) orientation %hd", name.c_str(), dimension, tag, lower.c_str(), domain.c_str(), orientation);
            domain = lower;

            if(surfaces.find(lower) == surfaces.end()) {       // New surface, not found yet.
                surfaces[lower].resize(surfaces[dirichletName].size(), 0);
                Logger::log(Logger::Color::green, Logger::Level::debug, true, "New domain");
            }
            #ifdef USE_OPENMP
            #pragma omp parallel for simd
            #endif
            for(size_t idx = 0; idx < nodesGroup.size(); ++idx) {
                if(surfaces[domain][nodesGroup[idx]-1] != 1)            // When in multiple groups, assume oriented. Then triangles are oriented correctly.
                    surfaces[domain][nodesGroup[idx]-1] = orientation;  // -1 in the index since GMSH tags start from 1.
                Logger::log(Logger::Color::reset, Logger::Level::debug, true, "Domain %s node %ld orientation: %hd", domain.c_str(), nodesGroup[idx], surfaces[domain][nodesGroup[idx]-1]);
            }
        }
    }
    static SurfaceList parseNodesPhysical(size_t const dofs     //!< Number of total dofs in the mesh (at least the nodes).
                                          ) {                   //! \brief Parse which nodes belong to which surface. See parsePhysical for the format. \return For each surface, for every dof, true/false depending whether it belongs to the surface.
        gmsh::vectorpair physicalGroups;        // Store the pairs (dimension, tag) for all physical groups
        gmsh::model::getPhysicalGroups(physicalGroups);

        SurfaceList surfaces;           // What dofs belongs to each surface.
        surfaces[dirichletName].resize(dofs, 0);// See parsePhysicalName. Put also triangles, to be filled later.
        surfaces[neumannName].resize(dofs, 0);  // ^

        for(std::pair<int,int> const &dimensionTag: physicalGroups) { // Parse the physical information into surfaces.
            std::string physicalName;       // Name of the physical tags, as output from gmsh.
            gmsh::model::getPhysicalName(dimensionTag.first, dimensionTag.second, physicalName);
            parsePhysicalName(physicalName, dimensionTag.second, dimensionTag.first, surfaces);
        }
        Logger::log(Logger::Color::green, Logger::Level::debug, true, "Found %ld surfaces", surfaces.size());
        for(auto surface : surfaces)
            Logger::log(Logger::Color::green, Logger::Level::debug, true, "Surface %s size %ld", surface.first.c_str(), surface.second.size());
        return surfaces;
    }

    static std::vector<size_t> reorderElements(Mesh::MembershipVector const & dirichletBoundary,//!< Membership to the dirichlet boundary.
                                               Mesh::MembershipVector const & neumannBoundary,  //!< Membership to the neumann boundary.
                                               size_t const nodesCount,                         //!< Total number of nodes.
                                               Mesh & meshData                                  //!<[out] \p triangleDofs and \p vertexDofs, the number of triangle and vertex degrees of freedom are filled. The membership data is reordered.
                                               ) {                                              //! Reorder the vertices and triangles so that the unknowns are the beginning of their vector section. \return Reverse order map: index j (of the appropriate array segment) stores the position of the element with tag j+1.
        std::vector<size_t> reordering(dirichletBoundary.size(),0);

        Logger::log(Logger::Color::green, Logger::Level::debug, "Computing reorder map\n");
        size_t nextDofPosition; // Index for the next DoF
        // Vertices: if in Dirichlet boundary, they are known, push them to the end; otherwise to the beginning.
        nextDofPosition = 0;
        meshData.vertexDofs_ = nodesCount;
        for(size_t vertex = 0; vertex < nodesCount; ++vertex)
            if(dirichletBoundary[vertex] != 0)
                reordering[vertex] = --meshData.vertexDofs_;
            else
                reordering[vertex] = nextDofPosition++;
        assert(nextDofPosition == meshData.vertexDofs_);

        // Triangles: if in Neumann boundary, they are known, push them to the end; otherwise to the beginning. 
        nextDofPosition = 0;                                            // Restart counting from 0
        meshData.triangleDofs_ = neumannBoundary.size() - nodesCount;   // Total number of triangles (totalDofs - nodes).
        for(size_t triangle = nodesCount; triangle < neumannBoundary.size(); ++triangle) 
            if(neumannBoundary[triangle] != 0)
                reordering[triangle] = --meshData.triangleDofs_;
            else
                reordering[triangle] = nextDofPosition++;
        assert(nextDofPosition == meshData.triangleDofs_);

        for(size_t idx = 0; idx < reordering.size(); ++idx) {
            Logger::log(Logger::Color::green, Logger::Level::debug, "Reorder %3ld %3ld.\n", idx, reordering[idx]);            
        }

        Logger::log(Logger::Color::green, Logger::Level::debug, "Reordering surfaces.\n");
        Mesh::MembershipVector tmp(reordering.size(), 0);
        for(Mesh::MembershipVector &surface: meshData.membership_) {     // Also reorder membership data.
            #ifdef USE_OPENMP
            #pragma omp parallel for simd
            #endif
            for(size_t idx = 0; idx < nodesCount; ++idx) {
                size_t position = reordering[idx];
                if(position >= meshData.vertexDofs_)
                    position += meshData.triangleDofs_;
                tmp[position] = surface[idx];
            }
            #ifdef USE_OPENMP
            #pragma omp parallel for simd
            #endif
            for(size_t idx = nodesCount; idx < tmp.size(); ++idx) {
                size_t position = reordering[idx];
                if(position < meshData.triangleDofs_)
                    position += meshData.vertexDofs_;
                else
                    position += nodesCount;
                tmp[position] = surface[idx];
            }
            tmp.swap(surface);
        }
        return reordering;
    }

    static surface3d* buildSurface(std::vector<size_t> const & reordering,      //!< Reverse reordering map. See reorderElements.
                                   std::vector<double> const & nodes,           //!< Coordinates for nodes, flattened.
                                   std::vector<size_t> const & edgesInTriangle, //!< [out] Edges of all triangles, as pairs of node tags. Triangle t has 6 indices starting at 6t.
                                   std::vector<size_t> const & edgeTags,        //!< [out] Edges of all triangles, as edge tags. Triangle t has 3 indices starting at 3t.
                                   std::vector<size_t> const & triangles,       //!< Triangles, as list of nodes (flattened).
                                   std::vector<size_t> const & nodeTags         //!< Tags of the nodes, in the same order as \p nodes.
                                   ) {                                          //! Consturct a surface in H2Lib format. \return The new surface.
        std::size_t edgesCount = 0;
        #ifdef USE_OPENMP
        #pragma omp parallel for reduction (max:edgesCount)
        #endif
        for(size_t edge = 0; edge < edgeTags.size(); ++edge) {  // Count the unique edges
            if(edgeTags[edge] > edgesCount)
                edgesCount = edgeTags[edge];
        }
        surface3d *fullSurface = new_surface3d(nodeTags.size(), edgesCount, (triangles.size()/3));
        Logger::log(Logger::Color::magenta, Logger::Level::ioBasic, "Created surface with %ld triangles, %ld edges, %ld nodes...", fullSurface->triangles, fullSurface->edges, fullSurface->vertices);

        Logger::log(Logger::Color::magenta, Logger::Level::ioBasic, "Filling node coordinates...");
        #ifdef USE_OPENMP
        #pragma omp parallel for simd
        #endif
        for(size_t vertex = 0; vertex < fullSurface->vertices; ++vertex)
            for(size_t coordinate = 0; coordinate < 3; ++coordinate)
                fullSurface->x[reordering[nodeTags[vertex]-1]][coordinate] = nodes[vertex*3 + coordinate];      // node with tag t goes to reordering[t-1]. The tag of the vertex at position 3p...3p+2 is given by nodetags[p].

        Logger::log(Logger::Color::magenta, Logger::Level::ioBasic, "Filling triangles and edges...");
        #ifdef USE_OPENMP
        #pragma omp parallel for simd
        #endif
        for(size_t triangle = 0; triangle < fullSurface->triangles; ++triangle) {
            size_t triangleMeshIdx = reordering[triangle + fullSurface->vertices];      // Triangle which is located in triangles at position 3t..3t+2 goes at reordering[t+vertices]
            for(size_t vertex = 0; vertex < 3; ++vertex) {
                size_t edgeTag = edgeTags[3*triangle + vertex]-1;       // Despite the name, make it 0-based.
                // Logger::log(Logger::Color::magenta, Logger::Level::ioBasic, "\nEdge %ld/%ld (%ld)", edgeTag, edgesCount, fullSurface->vertices);
                fullSurface->t[triangleMeshIdx][vertex] = reordering[triangles[triangle*3 + vertex]-1]; // Fill with reordered verted index. '-1' for tags to indices.
                size_t edgeVertex0 = reordering[edgesInTriangle[6*triangle + 2*vertex]-1], edgeVertex1 = reordering[edgesInTriangle[6*triangle + 2*vertex + 1]-1];
                if(edgeVertex0 < edgeVertex1) { // Each edge is filled more than once (exactly twice or more?). This condition ensures that the same value is written all times (avoids problems in case of write races).
                    fullSurface->e[edgeTag][0] = edgeVertex0;
                    fullSurface->e[edgeTag][1] = edgeVertex1;
                } else {
                    fullSurface->e[edgeTag][0] = edgeVertex1;
                    fullSurface->e[edgeTag][1] = edgeVertex0;
                }
                // Logger::log(Logger::Color::magenta, Logger::Level::ioBasic, "\nTriangle %ld (%ld), vertex %ld (%ld), updated edge %ld (%ld, %ld)", triangle, triangleMeshIdx, vertex, reordering[triangles[triangle*3 + vertex]-1], edgeTag, fullSurface->e[edgeTag][0], fullSurface->e[edgeTag][1]);
                
                [[maybe_unused]] bool found = false;
                for(size_t vertex2 = 0; vertex2 < 3; ++vertex2) // Determine to which vertex this edge is actually opposite.
                    if((fullSurface->e[edgeTag][0] != reordering[triangles[triangle*3 + vertex2]-1])&&
                       (fullSurface->e[edgeTag][1] != reordering[triangles[triangle*3 + vertex2]-1])) {
                        fullSurface->s[triangleMeshIdx][vertex2] = edgeTag;
                        found = true;
                        break;
                    }
                assert(found);
            }
        }
        Logger::log(Logger::Color::magenta, Logger::Level::ioBasic, "%d errors.\n", check_surface3d(fullSurface));
        return fullSurface;
    }

    static inline void getElements(std::array<std::vector<size_t>, 3> & tags,   //!< [out] Tags for all elements (nodes, triangles, lines).
                                   std::vector<double> &nodes,                  //!< [out] Nodes coordinates.
                                   std::vector<size_t> & edgesInTriangle,       //!< [out] Edges of all triangles, as pairs of node tags. Triangle t has 6 indices starting at 6t.
                                   std::vector<size_t> & edgeTags,              //!< [out] Edges of all triangles, as edge tags. Triangle t has 3 indices starting at 3t.
                                   std::vector<size_t> & triangles              //!< [out] List of tags identifying triangles.
                                   ) {                                          //! Parse and fill all element tags, and their definition.
        Logger::log(Logger::Color::magenta, Logger::Level::ioBasic, "Fetching data from GMSH APIs...");
        std::vector<double> ignoreDouble;       // Used for GMSH functions that require an output argument, which is however discarded.
        std::vector<int> ignoreInt;             // ^
        gmsh::model::mesh::renumberNodes();     // Make sure we need a vector, not a map.
        gmsh::model::mesh::renumberElements();  // ^

        int const gmshTriangleId = gmsh::model::mesh::getElementType("Triangle", 1);
        gmsh::model::mesh::getElementEdgeNodes(gmshTriangleId, edgesInTriangle);
        gmsh::model::mesh::createEdges();
        gmsh::model::mesh::getEdges(edgesInTriangle, edgeTags, ignoreInt); ignoreInt.clear();   // Tags of edges in the same order as edges identified as pairs. [3t <= .  < 3t+3] is triangle t.
        gmsh::model::mesh::getNodes(tags[nodesIdx], nodes, ignoreDouble); ignoreDouble.clear(); // Get all nodes
        // Logger::log(Logger::Color::magenta, Logger::Level::ioBasic, "Triangle id: %d; 6*triangles: %ld; 3*triangles: %ld; 3*nodes: %ld, nodes: %ld\n", gmshTriangleId, edgesInTriangle.size(), edgeTags.size(), nodes.size(), tags[nodesIdx].size());
        std::vector<int> elementTypes;
        std::vector<std::vector<size_t>> allElements;
        std::vector<std::vector<size_t>> allTags;
        gmsh::model::mesh::getElements(elementTypes, allTags, allElements);             // Get all elements, their types and tags.
        for(size_t idx = 0; idx < elementTypes.size(); ++idx)                           // Only get triangles
            if(elementTypes[idx] == gmshTriangleId) {
                tags[trianglesIdx] = allTags[idx];
                triangles = allElements[idx];
            }
    }
    static inline void splitTriangles(std::vector<size_t> const &triangles,     //!< Triangles description.
                                      size_t const nodesCount,                  //!< Number of nodes (offset for the triangles).
                                      SurfaceList &surfaces                     //!< [in,out] For each surface, for every mesh degree of freedom, true if it belongs to the surface.
                                      ) {                                       //! Fill information about which triangle is part of wich surface.
        Logger::log(Logger::Color::green, Logger::Level::debug, "Determining surfaces for triangles\n");
        for(size_t node = 0; node < nodesCount; ++node) {
            for(SurfaceList::value_type &surface: surfaces) {
                Logger::log(Logger::Color::yellow, Logger::Level::debug, false, "%s(v%3ld): %2hd\t", std::get<0>(surface).c_str(), node, std::get<1>(surface)[node]);
            }
            Logger::log(Logger::Color::yellow, Logger::Level::debug, "\n");
        }
        Logger::log(Logger::Color::yellow, Logger::Level::debug, "Triangle membership/orientation:\n");
        #ifdef USE_OPENMP
        #pragma omp parallel for simd
        #endif
        for(size_t triangle = 0; triangle < triangles.size()/3; ++triangle) {
            for(SurfaceList::value_type &surface: surfaces) {
                std::get<1>(surface)[triangle + nodesCount] = 1;     // Element to fill: does the triangle belong to the surface? (Iff all vertices do.)
                for(size_t point = 0; point < 3; ++point) {
                    short int orientation = std::get<1>(surface)[triangles[triangle*3 + point]-1];
                    if(orientation != 1) {      // So it's not correctly oriented, or not in the surface at all.
                        std::get<1>(surface)[triangle + nodesCount] = orientation;
                        if(orientation == 0)  // make sure the triangle is not in the surface if a vertex isn't.
                            break;
                    }
                }
                Logger::log(Logger::Color::yellow, Logger::Level::debug, false, "%s(%3ld): %2hd\t", std::get<0>(surface).c_str(), triangle, std::get<1>(surface)[triangle + nodesCount]);
            }
            Logger::log(Logger::Color::yellow, Logger::Level::debug, "\n");
        }
    }

    MeshReader::MeshReader() {
        gmsh::initialize();
        enabled_[getUnderlying(View::offset)] = saveOffset();
        enabled_[getUnderlying(View::rhs)] = saveRhs();
        enabled_[getUnderlying(View::solution)] = saveSolution();
        enabled_[getUnderlying(View::error)] = saveError();
        for(size_t idx = 0; idx < size<View>() ; ++idx) {
            if(enabled_[idx]) {
                viewTags_[2*idx] = gmsh::view::add(names[idx] + " Dirichlet");
                viewTags_[2*idx+1] = gmsh::view::add(names[idx] + " Neumann");
            } else {
                viewTags_[2*idx] = viewTags_[2*idx+1] = 0;
            }
        }
    }
    MeshReader::~MeshReader() {
        for(size_t idx = 0; idx < size<View>() ; ++idx) {
            Logger::log(Logger::Color::blue, Logger::Level::debug, false, "Writing view %ld. ", idx);
            if(enabled_[idx]) {
                std::string outname = fileName_.substr(fileName_.rfind('/') + 1);
                gmsh::view::write(viewTags_[2*idx],   outputDirectory() + names[idx] + "D_" + outname);
                gmsh::view::write(viewTags_[2*idx+1], outputDirectory() + names[idx] + "N_" + outname);
            }
        }
        Logger::log(Logger::Color::blue, Logger::Level::debug, true, "");
        gmsh::finalize();
    }
    void MeshReader::parseMesh(std::string const filename) {
        fileName_ = filename;
        Logger::log(Logger::Color::magenta, Logger::Level::ioBasic, true, "Loading mesh file %s", filename.c_str());
        gmsh::option::setNumber("General.Verbosity", 2);// Suppres info messages.
        gmsh::open(filename);   // Load file.
        gmsh::model::getCurrent(modelName_);    // Get the current model name.
        std::vector<size_t> edgesInTriangle, edgeTags, triangles;
        std::vector<double> nodes; std::array<std::vector<size_t>, 3> tags;
        getElements(tags, nodes, edgesInTriangle, edgeTags, triangles);  // Query GMSH about elements.
        SurfaceList surfaces = parseNodesPhysical((nodes.size() + triangles.size()) / 3);   // Generate information about nodes domain splitting.
        splitTriangles(triangles, nodes.size()/3, surfaces);       // Decide which triangles belong to which surface.
        for(SurfaceList::value_type &surface : surfaces)
            if((surface.first != dirichletName) && (surface.first != neumannName)) {
                mesh_.membership_.push_back({});
                mesh_.membership_.back().swap(surface.second);
                Logger::log(Logger::Color::green, Logger::Level::debug, true, "Added surface %s size %ld", surface.first.c_str(), mesh_.membership_.back().size());
            }

        reordering_ = reorderElements(surfaces[dirichletName], surfaces[neumannName], tags[nodesIdx].size(), mesh_);
        for(std::vector<int short> const & surfaceMembership : mesh_.membership_)
            Logger::log(Logger::Color::green, Logger::Level::debug, "Membership size (after reorder): %ld\n", surfaceMembership.size());

        Logger::log(Logger::Color::green, Logger::Level::debug, "Idx: Reorder Orient:\n");
        for(size_t idx = 0; idx < reordering_.size(); ++idx) {
            Logger::log(Logger::Color::green, Logger::Level::debug, "%3ld: %3ld\t %2hd\n", idx, reordering_[idx], mesh_.membership_.front()[idx]);
        }

        mesh_.mesh_ = buildSurface(reordering_, nodes, edgesInTriangle, edgeTags, triangles, tags[nodesIdx]);

        for(std::vector<int short> const & surfaceMembership : mesh_.membership_)
            Logger::log(Logger::Color::green, Logger::Level::debug, "Membership size (final): %ld\n", surfaceMembership.size());

        nodeTag_.resize(mesh_.vertices());
        elementTag_.resize(mesh_.triangles());
        #ifdef USE_OPENMP
        #pragma omp parallel for simd
        #endif
        for(size_t node = 0; node < mesh_.vertices(); ++node)
            nodeTag_[reordering_[tags[nodesIdx][node]-1]] = tags[nodesIdx][node];
        #ifdef USE_OPENMP
        #pragma omp parallel for simd
        #endif
        for(size_t triangle = 0; triangle < mesh_.triangles(); ++triangle)
            elementTag_[reordering_[triangle+mesh_.vertices()]] = tags[trianglesIdx][triangle];
    }
    void MeshReader::saveView(View const type, avector const *const data) {
        if(!enabled_[getUnderlying(type)])
            return;

        size_t &step = step_[getUnderlying(type)];
        int viewTagD = viewTags_[2*getUnderlying(type)];
        int viewTagN = viewTags_[2*getUnderlying(type)+1];
        Logger::log(Logger::Color::cyan, Logger::Level::ioBasic, true, "pletting type %d, step %ld", getUnderlying(type), step);

        std::vector<double> realDData(mesh_.vertices()), imagDData(mesh_.vertices()),
        realNData(mesh_.triangles()), imagNData(mesh_.triangles());

        #ifdef USE_OPENMP
        #pragma omp parallel for simd
        #endif
        for(size_t idx = 0; idx < mesh_.vertices(); ++idx) {
            size_t realIdx = idx;
            if(realIdx >= mesh_.vertexDofs_) {
                if(data->dim < mesh_.vertices() + mesh_.triangles()) {
                    realDData[idx] = 0;
                    imagDData[idx] = 0;
                    continue;
                }
                realIdx += mesh_.triangleDofs_;
            } 
            realDData[idx] = creal(getentry_avector(data, realIdx));
            imagDData[idx] = cimag(getentry_avector(data, realIdx));
        }
        #ifdef USE_OPENMP
        #pragma omp parallel for simd
        #endif
        for(size_t idx = 0; idx < mesh_.triangles(); ++idx) {
            size_t realIdx = idx;
            if(realIdx < mesh_.triangleDofs_) {
                realIdx += mesh_.vertexDofs_;
            } else {
                if(data->dim < mesh_.vertices() + mesh_.triangles()) {
                    realNData[idx] = 0;
                    imagNData[idx] = 0;
                    continue;
                }
                realIdx += mesh_.vertices();
            }
            realNData[idx] = creal(getentry_avector(data, realIdx));
            imagNData[idx] = cimag(getentry_avector(data, realIdx));
        }

        gmsh::view::addHomogeneousModelData(viewTagD, step  , modelName_, "NodeData", nodeTag_, realDData);
        gmsh::view::addHomogeneousModelData(viewTagD, step+1, modelName_, "NodeData", nodeTag_, imagDData);
        gmsh::view::addHomogeneousModelData(viewTagN, step  , modelName_, "ElementData", elementTag_, realNData);
        gmsh::view::addHomogeneousModelData(viewTagN, step+1, modelName_, "ElementData", elementTag_, imagNData);
        step += 2;
    }

    Mesh::~Mesh() {}
    Mesh::Mesh() {}
    Mesh::Mesh(Mesh && other)
        : membership_(other.membership_),
      vertexDofs_(other.vertexDofs_),
      triangleDofs_(other.triangleDofs_) {
      mesh_ = other.mesh_;
        other.mesh_ = nullptr;
    }
    Mesh& Mesh::operator=(Mesh &&other) {
        if(&other == this)
            return *this;

        membership_ = other.membership_;
        vertexDofs_ = other.vertexDofs_;
        triangleDofs_ = other.triangleDofs_;
        mesh_ = other.mesh_;
        other.mesh_ = nullptr;
        return *this;
    }
    size_t Mesh::triangles() const {
        assert(mesh_ != nullptr);
        return mesh_->triangles;
    }
    size_t Mesh::vertices() const {
        assert(mesh_ != nullptr);
        return mesh_->vertices;
    }
    SpaceVector<real const> Mesh::vertex(size_t const vertexIdx) const {
        assert(mesh_ != nullptr);
        return {mesh_->x[vertexIdx][0], mesh_->x[vertexIdx][1], mesh_->x[vertexIdx][2]};
    }
    SpaceVector<real const> Mesh::normal(size_t const triangleIdx) const {
        assert(mesh_ != nullptr);
        return {mesh_->n[triangleIdx][0], mesh_->n[triangleIdx][1], mesh_->n[triangleIdx][2]};
    }
    std::array<size_t,3> Mesh::triangle(size_t const triangleIdx) const {
        assert(mesh_ != nullptr);
        return {mesh_->t[triangleIdx][0], mesh_->t[triangleIdx][1], mesh_->t[triangleIdx][2]};
    }

    void Mesh::writeMesh() {
        if(!saveMesh())
            return;
        std::string fileName = std::string(getenv("HOME")) + "/Documents/dev/data/out/meshes/" + meshFile().substr(meshFile().rfind('/') + 1);
        Logger::log(Logger::Color::yellow, Logger::Level::debug, false, "Writing mesh to file %s: ", fileName.c_str());
        FILE *meshFile = fopen(fileName.c_str(), "w");
        if(meshFile == NULL) {
            Logger::log(Logger::Color::red, Logger::Level::warning, true, "error %d (%s) opening file.", errno, std::strerror(errno));
            return;
        }
        Logger::log(Logger::Color::green, Logger::Level::debug, true, "file open...");

        fprintf(meshFile, "%ld %ld %ld %ld\n", vertices(), triangles(), vertexDofs_, triangleDofs_);
        for(uint vertex = 0; vertex < vertices(); ++vertex)
            fprintf(meshFile, "%.2lf %.2lf %.2lf\n", mesh_->x[vertex][0], mesh_->x[vertex][1], mesh_->x[vertex][2]);
        for(uint triangle = 0; triangle < triangles(); ++triangle)
            fprintf(meshFile, "%d %d %d\n", mesh_->t[triangle][0], mesh_->t[triangle][1], mesh_->t[triangle][2]);

        fclose(meshFile);
    }
}
