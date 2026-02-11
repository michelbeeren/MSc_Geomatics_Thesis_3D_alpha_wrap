#include "octree.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace Thesis {

void collect_faces(OctreeCell& cell, const Tree& tree) {
    cell.faces.clear();
    std::vector<Primitive::Id> hits;
    tree.all_intersected_primitives(cell.bbox, std::back_inserter(hits));
    for (auto id : hits) cell.faces.push_back(id);
}

double normal_variance(const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals) {
    K::Vector_3 mean(0,0,0);
    for (auto f : faces) {
        auto n = normals.at(f);
        mean += n / std::sqrt(n.squared_length());
    }
    mean /= static_cast<double>(faces.size());

    double variance = 0.0;
    for (auto f : faces) {
        auto n = normals.at(f);
        n = n / std::sqrt(n.squared_length());
        variance += (n - mean).squared_length();
    }
    return variance / static_cast<double>(faces.size());
}

bool has_sharp_feature(const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals, double angle_threshold_rad) {
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& n1 = normals.at(faces[i]);
        for (size_t j = i + 1; j < faces.size(); ++j) {
            const auto& n2 = normals.at(faces[j]);
            double cos_angle = (n1 * n2) / std::sqrt(n1.squared_length() * n2.squared_length());
            if (std::acos(std::clamp(cos_angle, -1.0, 1.0)) > angle_threshold_rad) return true;
        }
    }
    return false;
}

double distance_to_surface(const OctreeCell& cell, const Tree& tree) {
    const auto& b = cell.bbox;
    K::Point_3 center(0.5 * (b.xmin() + b.xmax()), 0.5 * (b.ymin() + b.ymax()), 0.5 * (b.zmin() + b.zmax()));
    return std::sqrt(tree.squared_distance(center));
}

bool surface_is_simple(const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals, double max_angle_rad) {
    if (faces.size() < 2) return true;
    const auto& n0 = normals.at(faces[0]);
    for (auto f : faces) {
        const auto& n = normals.at(f);
        double cosang = (n * n0) / std::sqrt(n.squared_length() * n0.squared_length());
        if (std::acos(std::clamp(cosang, -1.0, 1.0)) > max_angle_rad) return false;
    }
    return true;
}

bool needs_refinement(const OctreeCell& cell, const std::map<Mesh::Face_index, K::Vector_3>& face_normals) {
    if (cell.depth >= 6 || cell.faces.empty()) return false;
    return !surface_is_simple(cell.faces, face_normals, 10.0 * M_PI / 180.0);
}

void subdivide(OctreeCell& cell) {
    const CGAL::Bbox_3& b = cell.bbox;
    double mx = 0.5 * (b.xmin() + b.xmax());
    double my = 0.5 * (b.ymin() + b.ymax());
    double mz = 0.5 * (b.zmin() + b.zmax());

    auto make_child = [&](double xmin, double ymin, double zmin, double xmax, double ymax, double zmax) {
        auto child = std::make_unique<OctreeCell>();
        child->bbox = CGAL::Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
        child->depth = cell.depth + 1;
        return child;
    };

    cell.children[0] = make_child(b.xmin(), b.ymin(), b.zmin(), mx, my, mz);
    cell.children[1] = make_child(mx, b.ymin(), b.zmin(), b.xmax(), my, mz);
    cell.children[2] = make_child(mx, my, b.zmin(), b.xmax(), b.ymax(), mz);
    cell.children[3] = make_child(b.xmin(), my, b.zmin(), mx, b.ymax(), mz);
    cell.children[4] = make_child(b.xmin(), b.ymin(), mz, mx, my, b.zmax());
    cell.children[5] = make_child(mx, b.ymin(), mz, b.xmax(), my, b.zmax());
    cell.children[6] = make_child(mx, my, mz, b.xmax(), b.ymax(), b.zmax());
    cell.children[7] = make_child(b.xmin(), my, mz, mx, b.ymax(), b.zmax());
}

void refine_cell(OctreeCell& cell, const Tree& tree, const std::map<Mesh::Face_index, K::Vector_3>& face_normals) {
    collect_faces(cell, tree);
    if (!needs_refinement(cell, face_normals)) return;
    subdivide(cell);
    for (auto& child : cell.children) refine_cell(*child, tree, face_normals);
}

bool is_leaf(const OctreeCell& cell) {
    for (const auto& c : cell.children) if (c) return false;
    return true;
}

void collect_leaves(const OctreeCell& cell, std::vector<CGAL::Bbox_3>& leaves) {
    if (is_leaf(cell)) { leaves.push_back(cell.bbox); return; }
    for (const auto& c : cell.children) if (c) collect_leaves(*c, leaves);
}

void write_off_wireframe(const std::vector<CGAL::Bbox_3>& cells, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) return;
    out << "OFF\n" << cells.size() * 8 << " " << cells.size() * 6 << " 0\n";
    for (const auto& b : cells) {
        out << b.xmin() << " " << b.ymin() << " " << b.zmin() << "\n" << b.xmax() << " " << b.ymin() << " " << b.zmin() << "\n";
        out << b.xmax() << " " << b.ymax() << " " << b.zmin() << "\n" << b.xmin() << " " << b.ymax() << " " << b.zmin() << "\n";
        out << b.xmin() << " " << b.ymin() << " " << b.zmax() << "\n" << b.xmax() << " " << b.ymin() << " " << b.zmax() << "\n";
        out << b.xmax() << " " << b.ymax() << " " << b.zmax() << "\n" << b.xmin() << " " << b.ymax() << " " << b.zmax() << "\n";
    }
    std::size_t v = 0;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        out << "4 " << v+0 << " " << v+1 << " " << v+2 << " " << v+3 << "\n";
        out << "4 " << v+4 << " " << v+5 << " " << v+6 << " " << v+7 << "\n";
        out << "4 " << v+0 << " " << v+1 << " " << v+5 << " " << v+4 << "\n";
        out << "4 " << v+1 << " " << v+2 << " " << v+6 << " " << v+5 << "\n";
        out << "4 " << v+2 << " " << v+3 << " " << v+7 << " " << v+6 << "\n";
        out << "4 " << v+3 << " " << v+0 << " " << v+4 << " " << v+7 << "\n";
        v += 8;
    }
}

void collect_intersecting_leaves(const OctreeCell& cell, std::vector<CGAL::Bbox_3>& out) {
    if (is_leaf(cell)) { if (!cell.faces.empty()) out.push_back(cell.bbox); return; }
    for (const auto& c : cell.children) if (c) collect_intersecting_leaves(*c, out);
}

bool intersects_mesh_bbox(const CGAL::Bbox_3& bbox, const Tree& tree) {
    // Explicitly check if the optional has a value
    return !!tree.any_intersection(bbox);
    // OR: return tree.any_intersection(bbox).has_value();
}

void find_max_depth(const OctreeCell& cell, int& max_depth) {
    max_depth = std::max(max_depth, cell.depth);
    for (const auto& c : cell.children) if (c) find_max_depth(*c, max_depth);
}

void find_max_intersecting_leaf_depth(const OctreeCell& cell, int& max_depth) {
    if (is_leaf(cell)) { if (!cell.faces.empty()) max_depth = std::max(max_depth, cell.depth); return; }
    for (const auto& c : cell.children) if (c) find_max_intersecting_leaf_depth(*c, max_depth);
}

void collect_cubical_leaves_at_depth(const OctreeCell& cell, int target_depth, std::vector<CGAL::Bbox_3>& out) {
    if (cell.depth == target_depth) { if (is_leaf(cell) && !cell.faces.empty()) out.push_back(cell.bbox); return; }
    for (const auto& c : cell.children) if (c) collect_cubical_leaves_at_depth(*c, target_depth, out);
}

void collect_finest_intersecting_leaves(const OctreeCell& cell, int target_depth, std::vector<CGAL::Bbox_3>& out) {
    if (is_leaf(cell)) { if (cell.depth == target_depth && !cell.faces.empty()) out.push_back(cell.bbox); return; }
    for (const auto& c : cell.children) if (c) collect_finest_intersecting_leaves(*c, target_depth, out);
}

void collect_intersecting_leaves_verified(const OctreeCell& cell, const Tree& tree, std::vector<CGAL::Bbox_3>& out) {
    if (is_leaf(cell)) { if (intersects_mesh_bbox(cell.bbox, tree)) out.push_back(cell.bbox); return; }
    for (const auto& c : cell.children) if (c) collect_intersecting_leaves_verified(*c, tree, out);
}

} // namespace Thesis