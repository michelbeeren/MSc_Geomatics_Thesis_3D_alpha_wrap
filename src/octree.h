//
// Created by Michel Beeren on 26/01/2026.
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

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Traits = CGAL::AABB_traits_3<K, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;

struct OctreeCell {
    CGAL::Bbox_3 bbox;
    int depth = 0; // Stored during refinement!
    std::vector<Mesh::Face_index> faces;
    std::array<std::unique_ptr<OctreeCell>, 8> children;

    OctreeCell() { for (auto& c : children) c = nullptr; }
    bool is_leaf() const {
        for (const auto& c : children) if (c) return false;
        return true;
    }
};

class Octree {
public:
    std::unique_ptr<OctreeCell> root;

    Octree(const CGAL::Bbox_3& root_bbox) {
        root = std::make_unique<OctreeCell>();
        root->bbox = root_bbox;
        root->depth = 0;
    }

    // This is the function your Alpha Wrap code will call
    double nearest_refinement_depth(const Point_3& p) const {
        if (!root) return 0.0;
        return search_recursive(root.get(), p);
    }

private:
    double search_recursive(const OctreeCell* cell, const Point_3& p) const {
        if (cell->is_leaf()) {
            return static_cast<double>(cell->depth); // Return the pre-stored depth
        }

        // Logic to find which child the point belongs to
        const CGAL::Bbox_3& b = cell->bbox;
        double mx = 0.5 * (b.xmin() + b.xmax());
        double my = 0.5 * (b.ymin() + b.ymax());
        double mz = 0.5 * (b.zmin() + b.zmax());

        int index = 0;
        if (p.x() >= mx) index += 1;
        if (p.y() >= my) index += 2;
        if (p.z() >= mz) index += 4;

        if (cell->children[index]) {
            return search_recursive(cell->children[index].get(), p);
        }

        return static_cast<double>(cell->depth);
    }
};
void collect_faces( OctreeCell& cell, const Tree& tree);
double normal_variance( const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals);
bool has_sharp_feature( const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals, double angle_threshold_rad );
double distance_to_surface( const OctreeCell& cell, const Tree& tree );
bool surface_is_simple( const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals, double max_angle_rad );
bool needs_refinement( const OctreeCell& cell, const std::map<Mesh::Face_index, K::Vector_3>& face_normals);
void subdivide(OctreeCell& cell);
void refine_cell(OctreeCell& cell, const Tree& tree, const std::map<Mesh::Face_index, K::Vector_3>& face_normals);
// struct LeafCell;
// void collect_leaves( const OctreeCell& cell, std::vector<LeafCell>& leaves );
// void write_off_wireframe( const std::vector<LeafCell>& cells, const std::string& filename);
bool is_leaf(const OctreeCell& cell);
// void collect_intersecting_leaves( const OctreeCell& cell, std::vector<LeafCell>& out );
bool intersects_mesh_bbox(const CGAL::Bbox_3& bbox, const Tree& tree);
// void collect_cubical_leaves_at_depth( const OctreeCell& cell, int target_depth, std::vector<LeafCell>& out );
void find_max_depth(const OctreeCell& cell, int& max_depth);
void find_max_intersecting_leaf_depth( const OctreeCell& cell, int& max_depth );
// void collect_finest_intersecting_leaves( const OctreeCell& cell, int target_depth, std::vector<LeafCell>& out );
// void collect_intersecting_leaves_verified( const OctreeCell& cell, const Tree& tree, std::vector<LeafCell>& out );
// Change all occurrences of std::vector<LeafCell> to std::vector<CGAL::Bbox_3>
void collect_leaves(const OctreeCell& cell, std::vector<CGAL::Bbox_3>& leaves);
void write_off_wireframe(const std::vector<CGAL::Bbox_3>& cells, const std::string& filename);
void collect_intersecting_leaves(const OctreeCell& cell, std::vector<CGAL::Bbox_3>& out);
void collect_cubical_leaves_at_depth(const OctreeCell& cell, int target_depth, std::vector<CGAL::Bbox_3>& out);
void collect_finest_intersecting_leaves(const OctreeCell& cell, int target_depth, std::vector<CGAL::Bbox_3>& out);
void collect_intersecting_leaves_verified(const OctreeCell& cell, const Tree& tree, std::vector<CGAL::Bbox_3>& out);

#endif //THESIS_OCTREE_H