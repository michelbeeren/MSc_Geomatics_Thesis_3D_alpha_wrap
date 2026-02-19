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
