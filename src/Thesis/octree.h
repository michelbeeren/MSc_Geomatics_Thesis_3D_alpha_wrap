// //
// // Created by Michel Beeren on 26/01/2026.
// //
//

#ifndef THESIS_OCTREE_H
#define THESIS_OCTREE_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <memory>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <deque>

namespace Thesis {

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Traits = CGAL::AABB_traits_3<K, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;

struct OctreeCell {
    CGAL::Bbox_3 bbox;
    int depth = 0;
    int ix = 0, iy = 0, iz = 0; // Integer coords for neighbor finding

    std::vector<Mesh::Face_index> faces;
    std::array<std::unique_ptr<OctreeCell>, 8> children;

    OctreeCell() { for (auto& c : children) c = nullptr; }
    bool is_leaf() const {
        for (const auto& c : children) if (c) return false;
        return true;
    }
};

// Key for the hash map to find neighbors by (depth, x, y, z)
struct CellKey {
    int depth, ix, iy, iz;
    bool operator==(const CellKey& o) const {
        return depth == o.depth && ix == o.ix && iy == o.iy && iz == o.iz;
    }
};

struct CellKeyHash {
    size_t operator()(const CellKey& k) const noexcept {
        size_t h = std::hash<int>{}(k.depth);
        h ^= std::hash<int>{}(k.ix) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>{}(k.iy) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>{}(k.iz) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

using LeafMap = std::unordered_map<CellKey, OctreeCell*, CellKeyHash>;
OctreeCell* find_leaf_neighbor_any_depth(const LeafMap& map, int depth, int nx, int ny, int nz);

// Core Functions
void collect_faces(OctreeCell& cell, const Tree& tree);
void subdivide(OctreeCell& cell);
void refine_cell(OctreeCell& cell, const Mesh& mesh, const Tree& tree, const std::map<Mesh::Face_index, K::Vector_3>& face_normals);
void balance_2to1(OctreeCell& root, const Tree& tree);
void collect_leaves_at_max_depth(const OctreeCell& cell, int target_depth, std::vector<CGAL::Bbox_3>& out, bool only_intersecting = true);

// Utilities
void collect_leaves_as_bboxes(const OctreeCell& cell, std::vector<CGAL::Bbox_3>& out, bool only_intersecting = false);
void write_off_wireframe(const std::vector<CGAL::Bbox_3>& cells, const std::string& filename);

} // namespace Thesis

#endif