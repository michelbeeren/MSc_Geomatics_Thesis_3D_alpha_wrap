#include "octree.h"
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace Thesis {

OctreeCell* find_leaf_neighbor_any_depth(const LeafMap& map, int depth, int nx, int ny, int nz) {
        int d = depth;
        int x = nx, y = ny, z = nz;

        while (d >= 0) {
            auto it = map.find({d, x, y, z});
            if (it != map.end()) return it->second;

            // Move to the parent's coordinates at a coarser level
            x >>= 1;
            y >>= 1;
            z >>= 1;
            --d;
        }
        return nullptr;
    }

void collect_faces(OctreeCell& cell, const Tree& tree) {
    cell.faces.clear();
    std::vector<Primitive::Id> hits;
    tree.all_intersected_primitives(cell.bbox, std::back_inserter(hits));
    for (auto id : hits) cell.faces.push_back(id);
}

void subdivide(OctreeCell& cell) {
    const CGAL::Bbox_3& b = cell.bbox;
    double mid[3] = {
        0.5 * (b.xmin() + b.xmax()),
        0.5 * (b.ymin() + b.ymax()),
        0.5 * (b.zmin() + b.zmax())
    };

    for (int i = 0; i < 8; ++i) {
        // Use bits to determine if this child is in the 'lower' or 'upper' half
        int dx = (i & 1);      // Bit 0: X
        int dy = (i & 2) >> 1; // Bit 1: Y
        int dz = (i & 4) >> 2; // Bit 2: Z

        double xmin = dx ? mid[0] : b.xmin();
        double xmax = dx ? b.xmax() : mid[0];
        double ymin = dy ? mid[1] : b.ymin();
        double ymax = dy ? b.ymax() : mid[1];
        double zmin = dz ? mid[2] : b.zmin();
        double zmax = dz ? b.zmax() : mid[2];

        auto ch = std::make_unique<OctreeCell>();
        ch->bbox = CGAL::Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
        ch->depth = cell.depth + 1;

        // Coordinates follow the bits exactly
        ch->ix = (cell.ix << 1) | dx;
        ch->iy = (cell.iy << 1) | dy;
        ch->iz = (cell.iz << 1) | dz;

        cell.children[i] = std::move(ch);
    }
}

bool is_sharp_concave_edge(const Mesh& mesh, Mesh::Halfedge_index h, const std::map<Mesh::Face_index, K::Vector_3>& normals, double min_angle_rad) {
    if (mesh.is_border(h) || mesh.is_border(mesh.opposite(h))) return false;

    const auto& n1 = normals.at(mesh.face(h));
    const auto& n2 = normals.at(mesh.face(mesh.opposite(h)));

    double dot = (n1 * n2) / std::sqrt(n1.squared_length() * n2.squared_length());
    if (dot > std::cos(min_angle_rad)) return false; // Too flat

    Point_3 p_opp = mesh.point(mesh.target(mesh.next(mesh.opposite(h))));
    Point_3 p_edge = mesh.point(mesh.target(h));
    return ((p_opp - p_edge) * n1) > 0.000001; // Concavity check
}

void refine_cell(OctreeCell& cell, const Mesh& mesh, const Tree& tree, const std::map<Mesh::Face_index, K::Vector_3>& face_normals) {
    collect_faces(cell, tree);
    if (cell.depth >= 9 || cell.faces.empty()) return;

    bool needs_split = false;
    double threshold = 10.0 * M_PI / 180.0;

    for (auto f : cell.faces) {
        for (auto h : mesh.halfedges_around_face(mesh.halfedge(f))) {
            if (is_sharp_concave_edge(mesh, h, face_normals, threshold)) {
                if (CGAL::do_intersect(cell.bbox, K::Segment_3(mesh.point(mesh.source(h)), mesh.point(mesh.target(h))))) {
                    needs_split = true;
                    break;
                }
            }
        }
        if (needs_split) break;
    }

    if (needs_split) {
        subdivide(cell);
        for (auto& child : cell.children) refine_cell(*child, mesh, tree, face_normals);
    }
}

// Helper for balancing: Find a leaf at specific coords or coarser
OctreeCell* find_leaf(const LeafMap& map, int d, int x, int y, int z) {
    while (d >= 0) {
        auto it = map.find({d, x, y, z});
        if (it != map.end()) return it->second;
        x >>= 1; y >>= 1; z >>= 1; d--;
    }
    return nullptr;
}

    void balance_2to1(OctreeCell& root, const Tree& tree) {
    LeafMap leaves;
    std::deque<OctreeCell*> q;

    // Helper to add leaves to map and queue
    auto register_leaf = [&](auto self, OctreeCell* c) -> void {
        if (c->is_leaf()) {
            leaves[{c->depth, c->ix, c->iy, c->iz}] = c;
            q.push_back(c);
        } else {
            for (auto& ch : c->children) self(self, ch.get());
        }
    };
    register_leaf(register_leaf, &root);

    const int offsets[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};

    while (!q.empty()) {
        OctreeCell* c = q.front();
        q.pop_front();
        if (!c->is_leaf()) continue;

        for (auto& off : offsets) {
            int nx = c->ix + off[0], ny = c->iy + off[1], nz = c->iz + off[2];
            int max_coord = (1 << c->depth) - 1;
            if (nx < 0 || ny < 0 || nz < 0 || nx > max_coord || ny > max_coord || nz > max_coord) continue;

            OctreeCell* n = find_leaf_neighbor_any_depth(leaves, c->depth, nx, ny, nz);
            if (!n) continue;

            // 2:1 Rule: If neighbor is too coarse, split it
            if (c->depth > n->depth + 1) {
                leaves.erase({n->depth, n->ix, n->iy, n->iz}); // Remove old leaf
                subdivide(*n);
                for (auto& child : n->children) {
                    collect_faces(*child, tree);
                    leaves[{child->depth, child->ix, child->iy, child->iz}] = child.get();
                    q.push_back(child.get()); // Check new children
                }
                q.push_back(c); // Re-check current cell against new finer neighbors
            }
        }
    }
}

void collect_leaves_as_bboxes(const OctreeCell& cell, std::vector<CGAL::Bbox_3>& out, bool only_intersecting) {
    if (cell.is_leaf()) {
        if (!only_intersecting || !cell.faces.empty()) out.push_back(cell.bbox);
        return;
    }
    for (auto& ch : cell.children) collect_leaves_as_bboxes(*ch, out, only_intersecting);
}

void write_off_wireframe(const std::vector<CGAL::Bbox_3>& cells, const std::string& filename) {
    std::ofstream out(filename);
    // Set precision to 17
    out << std::fixed << std::setprecision(17);
    out << "OFF\n" << cells.size() * 8 << " " << cells.size() * 6 << " 0\n";
    for (const auto& b : cells) {
        out << b.xmin() << " " << b.ymin() << " " << b.zmin() << "\n" << b.xmax() << " " << b.ymin() << " " << b.zmin() << "\n"
            << b.xmax() << " " << b.ymax() << " " << b.zmin() << "\n" << b.xmin() << " " << b.ymax() << " " << b.zmin() << "\n"
            << b.xmin() << " " << b.ymin() << " " << b.zmax() << "\n" << b.xmax() << " " << b.ymin() << " " << b.zmax() << "\n"
            << b.xmax() << " " << b.ymax() << " " << b.zmax() << "\n" << b.xmin() << " " << b.ymax() << " " << b.zmax() << "\n";
    }
    for (size_t i = 0; i < cells.size(); ++i) {
        size_t v = i * 8;
        out << "4 " << v << " " << v+1 << " " << v+2 << " " << v+3 << "\n"
            << "4 " << v+4 << " " << v+5 << " " << v+6 << " " << v+7 << "\n"
            << "4 " << v << " " << v+1 << " " << v+5 << " " << v+4 << "\n"
            << "4 " << v+1 << " " << v+2 << " " << v+6 << " " << v+5 << "\n"
            << "4 " << v+2 << " " << v+3 << " " << v+7 << " " << v+6 << "\n"
            << "4 " << v+3 << " " << v << " " << v+4 << " " << v+7 << "\n";
    }
}

void collect_leaves_at_max_depth(const OctreeCell& cell, int target_depth, std::vector<CGAL::Bbox_3>& out, bool only_intersecting) {
    // Base Case: We reached the target depth
    if (cell.depth == target_depth) {
        if (cell.is_leaf()) {
            if (!only_intersecting || !cell.faces.empty()) {
                out.push_back(cell.bbox);
            }
        }
        return;
    }

    // Recursive Case: Keep going deeper
    if (!cell.is_leaf()) {
        for (const auto& child : cell.children) {
            if (child) {
                collect_leaves_at_max_depth(*child, target_depth, out, only_intersecting);
            }
        }
    }
}

} // namespace Thesis

// #include "octree.h"
// #include <iostream>
// #include <fstream>
// #include <iomanip>
// #include <cmath>
// #include <algorithm>
//
// namespace Thesis {
//
// void collect_faces(OctreeCell& cell, const Tree& tree) {
//     cell.faces.clear();
//     std::vector<Primitive::Id> hits;
//     tree.all_intersected_primitives(cell.bbox, std::back_inserter(hits));
//     for (auto id : hits) cell.faces.push_back(id);
// }
//
// double normal_variance(const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals) {
//     K::Vector_3 mean(0,0,0);
//     for (auto f : faces) {
//         auto n = normals.at(f);
//         mean += n / std::sqrt(n.squared_length());
//     }
//     mean /= static_cast<double>(faces.size());
//
//     double variance = 0.0;
//     for (auto f : faces) {
//         auto n = normals.at(f);
//         n = n / std::sqrt(n.squared_length());
//         variance += (n - mean).squared_length();
//     }
//     return variance / static_cast<double>(faces.size());
// }
//
// bool has_sharp_feature(const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals, double angle_threshold_rad) {
//     for (size_t i = 0; i < faces.size(); ++i) {
//         const auto& n1 = normals.at(faces[i]);
//         for (size_t j = i + 1; j < faces.size(); ++j) {
//             const auto& n2 = normals.at(faces[j]);
//             double cos_angle = (n1 * n2) / std::sqrt(n1.squared_length() * n2.squared_length());
//             if (std::acos(std::clamp(cos_angle, -1.0, 1.0)) > angle_threshold_rad) return true;
//         }
//     }
//     return false;
// }
//
// double distance_to_surface(const OctreeCell& cell, const Tree& tree) {
//     const auto& b = cell.bbox;
//     K::Point_3 center(0.5 * (b.xmin() + b.xmax()), 0.5 * (b.ymin() + b.ymax()), 0.5 * (b.zmin() + b.zmax()));
//     return std::sqrt(tree.squared_distance(center));
// }
//
// bool surface_is_simple(const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals, double max_angle_rad) {
//     if (faces.size() < 2) return true;
//     const auto& n0 = normals.at(faces[0]);
//     for (auto f : faces) {
//         const auto& n = normals.at(f);
//         double cosang = (n * n0) / std::sqrt(n.squared_length() * n0.squared_length());
//         if (std::acos(std::clamp(cosang, -1.0, 1.0)) > max_angle_rad) return false;
//     }
//     return true;
// }
//
//     // Update your is_concave_edge to include a sharpness threshold
// bool is_sharp_concave_edge(const Mesh& mesh,
//                                Mesh::Halfedge_index h,
//                                const std::map<Mesh::Face_index, K::Vector_3>& normals,
//                                double min_angle_rad) {
//     if (mesh.is_border(h) || mesh.is_border(mesh.opposite(h))) return false;
//
//     Mesh::Face_index f1 = mesh.face(h);
//     Mesh::Face_index f2 = mesh.face(mesh.opposite(h));
//
//     const auto& n1 = normals.at(f1);
//     const auto& n2 = normals.at(f2);
//
//     // 1. Sharpness Check: Angle between normals
//     // n1 and n2 are unit vectors (or should be).
//     // Dot product = cos(theta). If cos(theta) < cos(threshold), the angle is bigger.
//     double cos_theta = (n1 * n2) / std::sqrt(n1.squared_length() * n2.squared_length());
//     if (cos_theta > std::cos(min_angle_rad)) {
//         return false; // Not sharp enough
//     }
//
//     // 2. Concavity Check: Does it bend inward?
//     Point_3 p_opposite = mesh.point(mesh.target(mesh.next(mesh.opposite(h))));
//     Point_3 p_on_edge = mesh.point(mesh.target(h));
//     K::Vector_3 v = p_opposite - p_on_edge;
//
//     return (v * n1) > 0.000001;
// }
//
// bool is_concave_edge(const Mesh& mesh,
//                          Mesh::Halfedge_index h,
//                          const std::map<Mesh::Face_index, K::Vector_3>& normals) {
//     // A border edge cannot be concave (it has no neighbor)
//     if (mesh.is_border(h) || mesh.is_border(mesh.opposite(h))) return false;
//
//     Mesh::Face_index f1 = mesh.face(h);
//     Mesh::Face_index f2 = mesh.face(mesh.opposite(h));
//
//     const auto& n1 = normals.at(f1);
//
//     // Get a point on the neighboring face f2 that is NOT on the shared edge
//     // mesh.next(mesh.opposite(h)) gives us the halfedge leading to the 'opposite' vertex
//     Point_3 p_opposite = mesh.point(mesh.target(mesh.next(mesh.opposite(h))));
//     Point_3 p_on_edge = mesh.point(mesh.target(h));
//
//     // Vector from the edge to the opposite vertex
//     K::Vector_3 v = p_opposite - p_on_edge;
//
//     // If the dot product is positive, the neighbor face is pointing "inward"
//     // relative to the current face's normal.
//     return (v * n1) > 0.000001; // Small epsilon for float precision
// }
//
// bool surface_is_simple_concave(const std::vector<Mesh::Face_index>& faces,
//                                 const std::map<Mesh::Face_index, K::Vector_3>& normals,
//                                 double max_angle_rad,
//                                 const std::map<Mesh::Face_index, std::set<Mesh::Face_index>>& adjacency_map) {
//     if (faces.size() < 2) return true;  // No adjacent faces to compare
//
//     // Iterate through each face and its neighbors
//     for (const auto& f : faces) {
//         const auto& n0 = normals.at(f);
//
//         // Get the neighbors of face f
//         auto neighbors = adjacency_map.at(f);
//
//         // Check each neighboring face
//         for (const auto& neighbor : neighbors) {
//             const auto& n1 = normals.at(neighbor);
//
//             // Compute the cosine of the angle between the normals of the two faces
//             double cosang = (n0 * n1) / std::sqrt(n0.squared_length() * n1.squared_length());
//             double angle = std::acos(std::clamp(cosang, -1.0, 1.0));
//
//             // Debug: Print angle and cosine to help diagnose the issue
//             std::cout << "Angle between faces: " << angle * 180.0 / M_PI << " degrees\n";
//
//             // If the angle is too small (close to zero), skip this neighbor and continue with the next one
//             double angle_degrees_ = angle * 180.0 / M_PI;
//             if (angle_degrees_ < 10) {
//                 std::cout << "Skipping refinement due to small angle: " << angle_degrees_ << " degrees\n";
//                 continue;  // Skip refinement for nearly parallel faces (angle < 2 degrees)
//             }
//
//             // If the angle is larger than the threshold
//             if (angle > max_angle_rad) {
//                 // If normals are facing towards each other (dot product < 0), it's concave
//                 if ((n0 * n1) < 0) {
//                     std::cout << "Concave sharp feature detected: " << angle_degrees_ << " degrees\n";
//                     return false;  // Concave sharp feature detected, split the cell
//                 }
//             }
//         }
//     }
//
//     return true;  // No concave feature found
// }
//
// bool needs_refinement(const OctreeCell& cell,
//                       const Mesh& mesh,
//                       const std::map<Mesh::Face_index, K::Vector_3>& face_normals,
//                       double threshold_rad) {
//
//     if (cell.depth >= 9 || cell.faces.empty()) return false;
//
//     for (auto f : cell.faces) {
//         for (auto h : mesh.halfedges_around_face(mesh.halfedge(f))) {
//
//             // Check if the edge is sharp AND concave
//             if (is_sharp_concave_edge(mesh, h, face_normals, threshold_rad)) {
//
//                 // Spatial check: is the concave edge inside this specific cell?
//                 Point_3 p1 = mesh.point(mesh.source(h));
//                 Point_3 p2 = mesh.point(mesh.target(h));
//
//                 if (CGAL::do_intersect(cell.bbox, K::Segment_3(p1, p2))) {
//                     return true;
//                 }
//             }
//         }
//     }
//     return false;
// }
//
//     // UPDATE THIS: You must pass 'mesh' down through refine_cell now
// void refine_cell(OctreeCell& cell,
//                      const Mesh& mesh, // Added Mesh parameter
//                      const Tree& tree,
//                      const std::map<Mesh::Face_index, K::Vector_3>& face_normals) {
//
//     collect_faces(cell, tree);
//
//     // Using 10 degrees as an example threshold
//     double threshold = 10.0 * M_PI / 180.0;
//
//     if (!needs_refinement(cell, mesh, face_normals, threshold)) return;
//
//     subdivide(cell);
//     for (auto& child : cell.children) {
//         if (child) refine_cell(*child, mesh, tree, face_normals);
//     }
// }
//
//
// // bool needs_refinement(const OctreeCell& cell,
// //                           const Mesh& mesh,
// //                           const std::map<Mesh::Face_index, K::Vector_3>& face_normals) {
// //     if (cell.depth >= 9 || cell.faces.empty()) return false;
// //
// //     for (auto f : cell.faces) {
// //         // Iterate over the 3 halfedges of this face
// //         for (auto h : mesh.halfedges_around_face(mesh.halfedge(f))) {
// //
// //             // 1. Check if this edge is concave
// //             if (is_concave_edge(mesh, h, face_normals)) {
// //
// //                 // 2. Spatial Check: Is this concave edge actually inside the current cell?
// //                 Point_3 p1 = mesh.point(mesh.source(h));
// //                 Point_3 p2 = mesh.point(mesh.target(h));
// //
// //                 // We refine only if the concave edge segment intersects the cell bbox
// //                 if (CGAL::do_intersect(cell.bbox, K::Segment_3(p1, p2))) {
// //                     return true; // Concave feature found inside this cell!
// //                 }
// //             }
// //         }
// //     }
// //     return false;
// // }
//
// void subdivide(OctreeCell& cell) {
//     const CGAL::Bbox_3& b = cell.bbox;
//     double mx = 0.5 * (b.xmin() + b.xmax());
//     double my = 0.5 * (b.ymin() + b.ymax());
//     double mz = 0.5 * (b.zmin() + b.zmax());
//
//     auto make_child = [&](int child,
//                           double xmin, double ymin, double zmin,
//                           double xmax, double ymax, double zmax)
//     {
//         auto ch = std::make_unique<OctreeCell>();
//         ch->bbox = CGAL::Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
//         ch->depth = cell.depth + 1;
//
//         int dx = (child & 1) ? 1 : 0;
//         int dy = (child & 2) ? 1 : 0;
//         int dz = (child & 4) ? 1 : 0;
//
//         ch->ix = cell.ix * 2 + dx;
//         ch->iy = cell.iy * 2 + dy;
//         ch->iz = cell.iz * 2 + dz;
//
//         return ch;
//     };
//
//     cell.children[0] = make_child(0, b.xmin(), b.ymin(), b.zmin(), mx, my, mz);
//     cell.children[1] = make_child(1, mx,      b.ymin(), b.zmin(), b.xmax(), my, mz);
//     cell.children[2] = make_child(2, mx,      my,      b.zmin(), b.xmax(), b.ymax(), mz);
//     cell.children[3] = make_child(3, b.xmin(), my,      b.zmin(), mx, b.ymax(), mz);
//     cell.children[4] = make_child(4, b.xmin(), b.ymin(), mz,      mx, my, b.zmax());
//     cell.children[5] = make_child(5, mx,      b.ymin(), mz,      b.xmax(), my, b.zmax());
//     cell.children[6] = make_child(6, mx,      my,      mz,      b.xmax(), b.ymax(), b.zmax());
//     cell.children[7] = make_child(7, b.xmin(), my,      mz,      mx, b.ymax(), b.zmax());
// }
//
// void collect_leaves(OctreeCell& c, LeafMap& map) {
//     if (c.is_leaf()) {
//         map[{c.depth, c.ix, c.iy, c.iz}] = &c;
//         return;
//     }
//     for (auto& ch : c.children) if (ch) collect_leaves(*ch, map);
// }
//
// OctreeCell* find_leaf_neighbor_any_depth(const LeafMap& map, int depth, int nx, int ny, int nz)
// {
//     // Try same depth, then walk up to coarser levels
//     int d = depth;
//     int x = nx, y = ny, z = nz;
//
//     while (d >= 0) {
//         auto it = map.find({d, x, y, z});
//         if (it != map.end()) return it->second;
//
//         // go one level coarser
//         x >>= 1; y >>= 1; z >>= 1;
//         --d;
//     }
//     return nullptr;
// }
//
// void balance_2to1(OctreeCell& root, const Mesh& mesh, const Tree& tree, const std::map<Mesh::Face_index, K::Vector_3>& face_normals)
// {
//     LeafMap leaves;
//     collect_leaves(root, leaves);
//
//     std::deque<OctreeCell*> q;
//     // q.reserve(leaves.size());
//     for (auto& kv : leaves) q.push_back(kv.second);
//
//     auto rebuild_leaf_map = [&](){
//         leaves.clear();
//         collect_leaves(root, leaves);
//     };
//
//     // 6-face neighbor offsets
//     const int off[6][3] = {
//         {+1,0,0},{-1,0,0},{0,+1,0},{0,-1,0},{0,0,+1},{0,0,-1}
//     };
//
//     while (!q.empty()) {
//         OctreeCell* c = q.front();
//         q.pop_front();
//         if (!c || !c->is_leaf()) continue;
//
//         for (auto& o : off) {
//             int nx = c->ix + o[0];
//             int ny = c->iy + o[1];
//             int nz = c->iz + o[2];
//
//             // outside domain at this depth?
//             int maxI = (1 << c->depth) - 1;
//             if (nx < 0 || ny < 0 || nz < 0 || nx > maxI || ny > maxI || nz > maxI)
//                 continue;
//
//             OctreeCell* n = find_leaf_neighbor_any_depth(leaves, c->depth, nx, ny, nz);
//             if (!n || !n->is_leaf()) continue;
//
//             int diff = c->depth - n->depth;
//             if (diff > 1) {
//                 // neighbor is too coarse, refine neighbor
//                 if (n->depth < 9) {
//                     // IMPORTANT: do NOT apply needs_refinement here.
//                     // This is "forced refinement" to satisfy 2:1 balance.
//                     subdivide(*n);
//
//                     // (optional but recommended) re-collect faces for new children:
//                     for (auto& ch : n->children) {
//                         if (ch) collect_faces(*ch, tree);
//                     }
//
//                     rebuild_leaf_map();
//
//                     // recheck around both regions
//                     q.push_back(c);
//                     for (auto& ch : n->children) if (ch) q.push_back(ch.get());
//                 }
//             } else if (diff < -1) {
//                 // c is too coarse compared to neighbor; refine c
//                 if (c->depth < 9) {
//                     subdivide(*c);
//                     for (auto& ch : c->children) {
//                         if (ch) collect_faces(*ch, tree);
//                     }
//                     rebuild_leaf_map();
//
//                     for (auto& ch : c->children) if (ch) q.push_back(ch.get());
//                 }
//             }
//         }
//     }
// }
//
//
// void subdivide2(OctreeCell& cell) {
//     const CGAL::Bbox_3& b = cell.bbox;
//     double mx = 0.5 * (b.xmin() + b.xmax());
//     double my = 0.5 * (b.ymin() + b.ymax());
//     double mz = 0.5 * (b.zmin() + b.zmax());
//
//     auto make_child = [&](double xmin, double ymin, double zmin, double xmax, double ymax, double zmax) {
//         auto child = std::make_unique<OctreeCell>();
//         child->bbox = CGAL::Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
//         child->depth = cell.depth + 1;
//         return child;
//     };
//
//     cell.children[0] = make_child(b.xmin(), b.ymin(), b.zmin(), mx, my, mz);
//     cell.children[1] = make_child(mx, b.ymin(), b.zmin(), b.xmax(), my, mz);
//     cell.children[2] = make_child(mx, my, b.zmin(), b.xmax(), b.ymax(), mz);
//     cell.children[3] = make_child(b.xmin(), my, b.zmin(), mx, b.ymax(), mz);
//     cell.children[4] = make_child(b.xmin(), b.ymin(), mz, mx, my, b.zmax());
//     cell.children[5] = make_child(mx, b.ymin(), mz, b.xmax(), my, b.zmax());
//     cell.children[6] = make_child(mx, my, mz, b.xmax(), b.ymax(), b.zmax());
//     cell.children[7] = make_child(b.xmin(), my, mz, mx, b.ymax(), b.zmax());
// }
//
// // void refine_cell(OctreeCell& cell, const Tree& tree, const std::map<Mesh::Face_index, K::Vector_3>& face_normals, const std::map<Mesh::Face_index, std::set<Mesh::Face_index>>& adjacency_map) {
// //     collect_faces(cell, tree);
// //     if (!needs_refinement(cell, face_normals, adjacency_map)) return;
// //     subdivide(cell);
// //     for (auto& child : cell.children) refine_cell(*child, tree, face_normals, adjacency_map);
// // }
//
// bool is_leaf(const OctreeCell& cell) {
//     for (const auto& c : cell.children) if (c) return false;
//     return true;
// }
//
// void collect_leaves(const OctreeCell& cell, std::vector<CGAL::Bbox_3>& leaves) {
//     if (is_leaf(cell)) { leaves.push_back(cell.bbox); return; }
//     for (const auto& c : cell.children) if (c) collect_leaves(*c, leaves);
// }
//
// void write_off_wireframe(const std::vector<CGAL::Bbox_3>& cells, const std::string& filename) {
//     std::ofstream out(filename);
//     if (!out) return;
//     out << std::fixed << std::setprecision(17);
//     out << "OFF\n" << cells.size() * 8 << " " << cells.size() * 6 << " 0\n";
//     for (const auto& b : cells) {
//         out << b.xmin() << " " << b.ymin() << " " << b.zmin() << "\n" << b.xmax() << " " << b.ymin() << " " << b.zmin() << "\n";
//         out << b.xmax() << " " << b.ymax() << " " << b.zmin() << "\n" << b.xmin() << " " << b.ymax() << " " << b.zmin() << "\n";
//         out << b.xmin() << " " << b.ymin() << " " << b.zmax() << "\n" << b.xmax() << " " << b.ymin() << " " << b.zmax() << "\n";
//         out << b.xmax() << " " << b.ymax() << " " << b.zmax() << "\n" << b.xmin() << " " << b.ymax() << " " << b.zmax() << "\n";
//     }
//     std::size_t v = 0;
//     for (std::size_t i = 0; i < cells.size(); ++i) {
//         out << "4 " << v+0 << " " << v+1 << " " << v+2 << " " << v+3 << "\n";
//         out << "4 " << v+4 << " " << v+5 << " " << v+6 << " " << v+7 << "\n";
//         out << "4 " << v+0 << " " << v+1 << " " << v+5 << " " << v+4 << "\n";
//         out << "4 " << v+1 << " " << v+2 << " " << v+6 << " " << v+5 << "\n";
//         out << "4 " << v+2 << " " << v+3 << " " << v+7 << " " << v+6 << "\n";
//         out << "4 " << v+3 << " " << v+0 << " " << v+4 << " " << v+7 << "\n";
//         v += 8;
//     }
// }
//
// void collect_intersecting_leaves(const OctreeCell& cell, std::vector<CGAL::Bbox_3>& out) {
//     if (is_leaf(cell)) { if (!cell.faces.empty()) out.push_back(cell.bbox); return; }
//     for (const auto& c : cell.children) if (c) collect_intersecting_leaves(*c, out);
// }
//
// bool intersects_mesh_bbox(const CGAL::Bbox_3& bbox, const Tree& tree) {
//     // Explicitly check if the optional has a value
//     return !!tree.any_intersection(bbox);
//     // OR: return tree.any_intersection(bbox).has_value();
// }
//
// void find_max_depth(const OctreeCell& cell, int& max_depth) {
//     max_depth = std::max(max_depth, cell.depth);
//     for (const auto& c : cell.children) if (c) find_max_depth(*c, max_depth);
// }
//
// void find_max_intersecting_leaf_depth(const OctreeCell& cell, int& max_depth) {
//     if (is_leaf(cell)) { if (!cell.faces.empty()) max_depth = std::max(max_depth, cell.depth); return; }
//     for (const auto& c : cell.children) if (c) find_max_intersecting_leaf_depth(*c, max_depth);
// }
//
// void collect_cubical_leaves_at_depth(const OctreeCell& cell, int target_depth, std::vector<CGAL::Bbox_3>& out) {
//     if (cell.depth == target_depth) { if (is_leaf(cell) && !cell.faces.empty()) out.push_back(cell.bbox); return; }
//     for (const auto& c : cell.children) if (c) collect_cubical_leaves_at_depth(*c, target_depth, out);
// }
//
// void collect_finest_intersecting_leaves(const OctreeCell& cell, int target_depth, std::vector<CGAL::Bbox_3>& out) {
//     if (is_leaf(cell)) { if (cell.depth == target_depth && !cell.faces.empty()) out.push_back(cell.bbox); return; }
//     for (const auto& c : cell.children) if (c) collect_finest_intersecting_leaves(*c, target_depth, out);
// }
//
// void collect_intersecting_leaves_verified(const OctreeCell& cell, const Tree& tree, std::vector<CGAL::Bbox_3>& out) {
//     if (is_leaf(cell)) { if (intersects_mesh_bbox(cell.bbox, tree)) out.push_back(cell.bbox); return; }
//     for (const auto& c : cell.children) if (c) collect_intersecting_leaves_verified(*c, tree, out);
// }
//
// } // namespace Thesis