//
// Created by Michel Beeren on 26/01/2026.
//
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <map>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Real_timer.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <boost/optional.hpp>
#include <boost/variant/get.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
using face_descriptor = Mesh::Face_index;
using Ray_3   = K::Ray_3;
using Segment_3 = K::Segment_3;
using Primitive   = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using AABB_traits = CGAL::AABB_traits<K, Primitive>;
using Tree        = CGAL::AABB_tree<AABB_traits>;

#include "octree.h"

// --------------------------------OCTREE REFINEMENT-----------------------------------

struct OctreeCell {
    CGAL::Bbox_3 bbox;
    int depth = 0;
    std::vector<Mesh::Face_index> faces;
    std::array<std::unique_ptr<OctreeCell>, 8> children;

    OctreeCell() {
        for (auto& c : children)
            c = nullptr;
    }
};

void collect_faces( OctreeCell& cell, const Tree& tree) {
    cell.faces.clear();
    std::vector<Primitive::Id> hits;
    tree.all_intersected_primitives(cell.bbox, std::back_inserter(hits));

    for (auto id : hits)
        cell.faces.push_back(id);
}

double normal_variance( const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals) {
    K::Vector_3 mean(0,0,0);

    for (auto f : faces) {
        auto n = normals.at(f);
        n = n / std::sqrt(n.squared_length());
        mean += n;
    }

    mean /= faces.size();

    double variance = 0.0;
    for (auto f : faces) {
        auto n = normals.at(f);
        n = n / std::sqrt(n.squared_length());
        variance += (n - mean).squared_length();
    }

    return variance / faces.size();
}

bool has_sharp_feature( const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals, double angle_threshold_rad ) {
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& n1 = normals.at(faces[i]);

        for (size_t j = i + 1; j < faces.size(); ++j) {
            const auto& n2 = normals.at(faces[j]);

            double cos_angle =
                (n1 * n2) / std::sqrt(n1.squared_length() * n2.squared_length());

            cos_angle = std::clamp(cos_angle, -1.0, 1.0);

            double angle = std::acos(cos_angle);

            if (angle > angle_threshold_rad)
                return true;
        }
    }
    return false;
}

double distance_to_surface( const OctreeCell& cell, const Tree& tree ) {
    const auto& b = cell.bbox;

    K::Point_3 center(
        0.5 * (b.xmin() + b.xmax()),
        0.5 * (b.ymin() + b.ymax()),
        0.5 * (b.zmin() + b.zmax())
    );

    return std::sqrt(tree.squared_distance(center));
}

bool surface_is_simple( const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals, double max_angle_rad ) {
    if (faces.size() < 2)
        return true;

    const auto& n0 = normals.at(faces[0]);

    for (auto f : faces) {
        const auto& n = normals.at(f);
        double cosang = (n * n0) /
            std::sqrt(n.squared_length() * n0.squared_length());

        cosang = std::clamp(cosang, -1.0, 1.0);

        if (std::acos(cosang) > max_angle_rad)
            return false;
    }
    return true;
}

bool needs_refinement( const OctreeCell& cell, const std::map<Mesh::Face_index, K::Vector_3>& face_normals) {
    constexpr int max_depth = 9;
    constexpr double cone_angle =
        10.0 * M_PI / 180.0; // strict!

    if (cell.depth >= max_depth)
        return false;

    // 1. Empty cell → stop
    if (cell.faces.empty())
        return false;

    // 2. Surface is locally simple → stop
    if (surface_is_simple(cell.faces, face_normals, cone_angle))
        return false;

    // 3. Otherwise → refine
    return true;
}

void subdivide(OctreeCell& cell)
{
    // ------------------------------------------------------------------
    // Preconditions (important for sanity)
    // ------------------------------------------------------------------
    const CGAL::Bbox_3& b = cell.bbox;

    const double dx = b.xmax() - b.xmin();
    const double dy = b.ymax() - b.ymin();
    const double dz = b.zmax() - b.zmin();

    // This MUST hold for a true octree
    // (you can comment this out later)
    if (std::abs(dx - dy) > 1e-9 || std::abs(dy - dz) > 1e-9) {
        std::cerr << "[subdivide] Warning: non-cubic cell at depth "
                  << cell.depth << " ("
                  << dx << ", " << dy << ", " << dz << ")\n";
    }

    // ------------------------------------------------------------------
    // Compute midpoints
    // ------------------------------------------------------------------
    const double mx = 0.5 * (b.xmin() + b.xmax());
    const double my = 0.5 * (b.ymin() + b.ymax());
    const double mz = 0.5 * (b.zmin() + b.zmax());

    // ------------------------------------------------------------------
    // Helper to create a child cell
    // ------------------------------------------------------------------
    auto make_child =
        [&](double xmin, double ymin, double zmin,
            double xmax, double ymax, double zmax)
    {
        auto child = std::make_unique<OctreeCell>();
        child->bbox  = CGAL::Bbox_3(xmin, ymin, zmin,
                                    xmax, ymax, zmax);
        child->depth = cell.depth + 1;
        return child;
    };

    // ------------------------------------------------------------------
    // Create the 8 octants (standard octree layout)
    // ------------------------------------------------------------------
    cell.children[0] = make_child(b.xmin(), b.ymin(), b.zmin(), mx, my, mz);
    cell.children[1] = make_child(mx,       b.ymin(), b.zmin(), b.xmax(), my, mz);
    cell.children[2] = make_child(mx,       my,       b.zmin(), b.xmax(), b.ymax(), mz);
    cell.children[3] = make_child(b.xmin(), my,       b.zmin(), mx, b.ymax(), mz);

    cell.children[4] = make_child(b.xmin(), b.ymin(), mz, mx, my, b.zmax());
    cell.children[5] = make_child(mx,       b.ymin(), mz, b.xmax(), my, b.zmax());
    cell.children[6] = make_child(mx,       my,       mz, b.xmax(), b.ymax(), b.zmax());
    cell.children[7] = make_child(b.xmin(), my,       mz, mx, b.ymax(), b.zmax());
}

void refine_cell(OctreeCell& cell, const Tree& tree, const std::map<Mesh::Face_index, K::Vector_3>& face_normals) {
    collect_faces(cell, tree);

    if (!needs_refinement(cell, face_normals))
        return;

    subdivide(cell);

    for (auto& child : cell.children)
        refine_cell(*child, tree, face_normals);
}

struct LeafCell {
    CGAL::Bbox_3 bbox;
};

void collect_leaves( const OctreeCell& cell, std::vector<LeafCell>& leaves ) {
    bool is_leaf = true;
    for (const auto& c : cell.children)
        if (c) is_leaf = false;

    if (is_leaf) {
        leaves.push_back({cell.bbox});
        return;
    }

    for (const auto& c : cell.children)
        if (c) collect_leaves(*c, leaves);
}

void write_off_wireframe( const std::vector<LeafCell>& cells, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: cannot open output file " << filename << "\n";
        return;
    }

    // set write precision higher
    out << std::fixed << std::setprecision(17);

    constexpr int VERTS_PER_CELL = 8;
    constexpr int FACES_PER_CELL = 6; // 6 quad faces per box

    const std::size_t numVertices = cells.size() * VERTS_PER_CELL;
    const std::size_t numFaces    = cells.size() * FACES_PER_CELL;

    // ------------------------------------------------------------------
    // OFF header
    // ------------------------------------------------------------------
    out << "OFF\n";
    out << numVertices << " " << numFaces << " 0\n";

    // ------------------------------------------------------------------
    // Write vertices
    // ------------------------------------------------------------------
    for (const auto& c : cells) {
        const auto& b = c.bbox;

        const double x0 = b.xmin(), x1 = b.xmax();
        const double y0 = b.ymin(), y1 = b.ymax();
        const double z0 = b.zmin(), z1 = b.zmax();

        out << x0 << " " << y0 << " " << z0 << "\n"; // 0
        out << x1 << " " << y0 << " " << z0 << "\n"; // 1
        out << x1 << " " << y1 << " " << z0 << "\n"; // 2
        out << x0 << " " << y1 << " " << z0 << "\n"; // 3
        out << x0 << " " << y0 << " " << z1 << "\n"; // 4
        out << x1 << " " << y0 << " " << z1 << "\n"; // 5
        out << x1 << " " << y1 << " " << z1 << "\n"; // 6
        out << x0 << " " << y1 << " " << z1 << "\n"; // 7
    }

    // ------------------------------------------------------------------
    // Write quad faces (boxes)
    // ------------------------------------------------------------------
    std::size_t v = 0;
    for (std::size_t i = 0; i < cells.size(); ++i) {

        // bottom
        out << "4 " << v+0 << " " << v+1 << " " << v+2 << " " << v+3 << "\n";
        // top
        out << "4 " << v+4 << " " << v+5 << " " << v+6 << " " << v+7 << "\n";

        // sides
        out << "4 " << v+0 << " " << v+1 << " " << v+5 << " " << v+4 << "\n";
        out << "4 " << v+1 << " " << v+2 << " " << v+6 << " " << v+5 << "\n";
        out << "4 " << v+2 << " " << v+3 << " " << v+7 << " " << v+6 << "\n";
        out << "4 " << v+3 << " " << v+0 << " " << v+4 << " " << v+7 << "\n";

        v += VERTS_PER_CELL;
    }

    out.close();
}

bool is_leaf(const OctreeCell& cell) {
    for (const auto& c : cell.children)
        if (c)
            return false;
    return true;
}

void collect_intersecting_leaves( const OctreeCell& cell, std::vector<LeafCell>& out ) {
    // Case 1: leaf cell
    if (is_leaf(cell)) {
        // Only keep it if it intersects the mesh
        if (!cell.faces.empty()) {
            out.push_back({cell.bbox});
        }
        return;
    }

    // Case 2: internal node → recurse
    for (const auto& c : cell.children) {
        if (c)
            collect_intersecting_leaves(*c, out);
    }
}

bool intersects_mesh_bbox(const CGAL::Bbox_3& bbox, const Tree& tree) {
    std::vector<Primitive::Id> hits;
    tree.all_intersected_primitives(bbox, std::back_inserter(hits));
    return !hits.empty();
}

void collect_cubical_leaves_at_depth( const OctreeCell& cell, int target_depth, std::vector<LeafCell>& out ) {
    // Leaf at the right depth
    if (cell.depth == target_depth) {
        bool is_leaf = true;
        for (const auto& c : cell.children)
            if (c) { is_leaf = false; break; }

        if (is_leaf && !cell.faces.empty()) {
            out.push_back({cell.bbox});
        }
        return;
    }

    // Recurse
    for (const auto& c : cell.children)
        if (c)
            collect_cubical_leaves_at_depth(*c, target_depth, out);
}

void find_max_depth(const OctreeCell& cell, int& max_depth) {
    max_depth = std::max(max_depth, cell.depth);
    for (const auto& c : cell.children)
        if (c)
            find_max_depth(*c, max_depth);
}

void find_max_intersecting_leaf_depth( const OctreeCell& cell, int& max_depth ) {
    bool is_leaf = true;
    for (const auto& c : cell.children)
        if (c) { is_leaf = false; break; }

    if (is_leaf) {
        if (!cell.faces.empty()) {
            max_depth = std::max(max_depth, cell.depth);
        }
        return;
    }

    for (const auto& c : cell.children)
        if (c)
            find_max_intersecting_leaf_depth(*c, max_depth);
}

void collect_finest_intersecting_leaves( const OctreeCell& cell, int target_depth, std::vector<LeafCell>& out ) {
    bool is_leaf = true;
    for (const auto& c : cell.children)
        if (c) { is_leaf = false; break; }

    if (is_leaf) {
        if (cell.depth == target_depth && !cell.faces.empty()) {
            out.push_back({cell.bbox});
        }
        return;
    }

    for (const auto& c : cell.children)
        if (c)
            collect_finest_intersecting_leaves(*c, target_depth, out);
}

void collect_intersecting_leaves_verified( const OctreeCell& cell, const Tree& tree, std::vector<LeafCell>& out ){
    bool is_leaf = true;
    for (const auto& c : cell.children)
        if (c) { is_leaf = false; break; }

    if (is_leaf) {
        if (intersects_mesh_bbox(cell.bbox, tree)) {
            out.push_back({cell.bbox});
        }
        return;
    }

    for (const auto& c : cell.children)
        if (c)
            collect_intersecting_leaves_verified(*c, tree, out);
}