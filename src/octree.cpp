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

    // Update your is_concave_edge to include a sharpness threshold
bool is_sharp_concave_edge(const Mesh& mesh,
                               Mesh::Halfedge_index h,
                               const std::map<Mesh::Face_index, K::Vector_3>& normals,
                               double min_angle_rad) {
    if (mesh.is_border(h) || mesh.is_border(mesh.opposite(h))) return false;

    Mesh::Face_index f1 = mesh.face(h);
    Mesh::Face_index f2 = mesh.face(mesh.opposite(h));

    const auto& n1 = normals.at(f1);
    const auto& n2 = normals.at(f2);

    // 1. Sharpness Check: Angle between normals
    // n1 and n2 are unit vectors (or should be).
    // Dot product = cos(theta). If cos(theta) < cos(threshold), the angle is bigger.
    double cos_theta = (n1 * n2) / std::sqrt(n1.squared_length() * n2.squared_length());
    if (cos_theta > std::cos(min_angle_rad)) {
        return false; // Not sharp enough
    }

    // 2. Concavity Check: Does it bend inward?
    Point_3 p_opposite = mesh.point(mesh.target(mesh.next(mesh.opposite(h))));
    Point_3 p_on_edge = mesh.point(mesh.target(h));
    K::Vector_3 v = p_opposite - p_on_edge;

    return (v * n1) > 0.000001;
}

bool is_concave_edge(const Mesh& mesh,
                         Mesh::Halfedge_index h,
                         const std::map<Mesh::Face_index, K::Vector_3>& normals) {
    // A border edge cannot be concave (it has no neighbor)
    if (mesh.is_border(h) || mesh.is_border(mesh.opposite(h))) return false;

    Mesh::Face_index f1 = mesh.face(h);
    Mesh::Face_index f2 = mesh.face(mesh.opposite(h));

    const auto& n1 = normals.at(f1);

    // Get a point on the neighboring face f2 that is NOT on the shared edge
    // mesh.next(mesh.opposite(h)) gives us the halfedge leading to the 'opposite' vertex
    Point_3 p_opposite = mesh.point(mesh.target(mesh.next(mesh.opposite(h))));
    Point_3 p_on_edge = mesh.point(mesh.target(h));

    // Vector from the edge to the opposite vertex
    K::Vector_3 v = p_opposite - p_on_edge;

    // If the dot product is positive, the neighbor face is pointing "inward"
    // relative to the current face's normal.
    return (v * n1) > 0.000001; // Small epsilon for float precision
}

bool surface_is_simple_concave(const std::vector<Mesh::Face_index>& faces,
                                const std::map<Mesh::Face_index, K::Vector_3>& normals,
                                double max_angle_rad,
                                const std::map<Mesh::Face_index, std::set<Mesh::Face_index>>& adjacency_map) {
    if (faces.size() < 2) return true;  // No adjacent faces to compare

    // Iterate through each face and its neighbors
    for (const auto& f : faces) {
        const auto& n0 = normals.at(f);

        // Get the neighbors of face f
        auto neighbors = adjacency_map.at(f);

        // Check each neighboring face
        for (const auto& neighbor : neighbors) {
            const auto& n1 = normals.at(neighbor);

            // Compute the cosine of the angle between the normals of the two faces
            double cosang = (n0 * n1) / std::sqrt(n0.squared_length() * n1.squared_length());
            double angle = std::acos(std::clamp(cosang, -1.0, 1.0));

            // Debug: Print angle and cosine to help diagnose the issue
            std::cout << "Angle between faces: " << angle * 180.0 / M_PI << " degrees\n";

            // If the angle is too small (close to zero), skip this neighbor and continue with the next one
            double angle_degrees_ = angle * 180.0 / M_PI;
            if (angle_degrees_ < 10) {
                std::cout << "Skipping refinement due to small angle: " << angle_degrees_ << " degrees\n";
                continue;  // Skip refinement for nearly parallel faces (angle < 2 degrees)
            }

            // If the angle is larger than the threshold
            if (angle > max_angle_rad) {
                // If normals are facing towards each other (dot product < 0), it's concave
                if ((n0 * n1) < 0) {
                    std::cout << "Concave sharp feature detected: " << angle_degrees_ << " degrees\n";
                    return false;  // Concave sharp feature detected, split the cell
                }
            }
        }
    }

    return true;  // No concave feature found
}

bool needs_refinement(const OctreeCell& cell,
                      const Mesh& mesh,
                      const std::map<Mesh::Face_index, K::Vector_3>& face_normals,
                      double threshold_rad) {

    if (cell.depth >= 9 || cell.faces.empty()) return false;

    for (auto f : cell.faces) {
        for (auto h : mesh.halfedges_around_face(mesh.halfedge(f))) {

            // Check if the edge is sharp AND concave
            if (is_sharp_concave_edge(mesh, h, face_normals, threshold_rad)) {

                // Spatial check: is the concave edge inside this specific cell?
                Point_3 p1 = mesh.point(mesh.source(h));
                Point_3 p2 = mesh.point(mesh.target(h));

                if (CGAL::do_intersect(cell.bbox, K::Segment_3(p1, p2))) {
                    return true;
                }
            }
        }
    }
    return false;
}

    // UPDATE THIS: You must pass 'mesh' down through refine_cell now
void refine_cell(OctreeCell& cell,
                     const Mesh& mesh, // Added Mesh parameter
                     const Tree& tree,
                     const std::map<Mesh::Face_index, K::Vector_3>& face_normals) {

    collect_faces(cell, tree);

    // Using 10 degrees as an example threshold
    double threshold = 10.0 * M_PI / 180.0;

    if (!needs_refinement(cell, mesh, face_normals, threshold)) return;

    subdivide(cell);
    for (auto& child : cell.children) {
        if (child) refine_cell(*child, mesh, tree, face_normals);
    }
}


// bool needs_refinement(const OctreeCell& cell,
//                           const Mesh& mesh,
//                           const std::map<Mesh::Face_index, K::Vector_3>& face_normals) {
//     if (cell.depth >= 9 || cell.faces.empty()) return false;
//
//     for (auto f : cell.faces) {
//         // Iterate over the 3 halfedges of this face
//         for (auto h : mesh.halfedges_around_face(mesh.halfedge(f))) {
//
//             // 1. Check if this edge is concave
//             if (is_concave_edge(mesh, h, face_normals)) {
//
//                 // 2. Spatial Check: Is this concave edge actually inside the current cell?
//                 Point_3 p1 = mesh.point(mesh.source(h));
//                 Point_3 p2 = mesh.point(mesh.target(h));
//
//                 // We refine only if the concave edge segment intersects the cell bbox
//                 if (CGAL::do_intersect(cell.bbox, K::Segment_3(p1, p2))) {
//                     return true; // Concave feature found inside this cell!
//                 }
//             }
//         }
//     }
//     return false;
// }

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

// void refine_cell(OctreeCell& cell, const Tree& tree, const std::map<Mesh::Face_index, K::Vector_3>& face_normals, const std::map<Mesh::Face_index, std::set<Mesh::Face_index>>& adjacency_map) {
//     collect_faces(cell, tree);
//     if (!needs_refinement(cell, face_normals, adjacency_map)) return;
//     subdivide(cell);
//     for (auto& child : cell.children) refine_cell(*child, tree, face_normals, adjacency_map);
// }

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
    out << std::fixed << std::setprecision(17);
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