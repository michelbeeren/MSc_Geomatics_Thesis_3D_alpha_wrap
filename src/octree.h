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

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;

using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Traits    = CGAL::AABB_traits_3<K, Primitive>;
using Tree      = CGAL::AABB_tree<Traits>;

struct OctreeCell;
void collect_faces( OctreeCell& cell, const Tree& tree);
double normal_variance( const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals);
bool has_sharp_feature( const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals, double angle_threshold_rad );
double distance_to_surface( const OctreeCell& cell, const Tree& tree );
bool surface_is_simple( const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals, double max_angle_rad );
bool needs_refinement( const OctreeCell& cell, const std::map<Mesh::Face_index, K::Vector_3>& face_normals);
void subdivide(OctreeCell& cell);
void refine_cell(OctreeCell& cell, const Tree& tree, const std::map<Mesh::Face_index, K::Vector_3>& face_normals);
struct LeafCell;
void collect_leaves( const OctreeCell& cell, std::vector<LeafCell>& leaves );
void write_off_wireframe( const std::vector<LeafCell>& cells, const std::string& filename);
bool is_leaf(const OctreeCell& cell);
void collect_intersecting_leaves( const OctreeCell& cell, std::vector<LeafCell>& out );
bool intersects_mesh_bbox(const CGAL::Bbox_3& bbox, const Tree& tree);
void collect_cubical_leaves_at_depth( const OctreeCell& cell, int target_depth, std::vector<LeafCell>& out );
void find_max_depth(const OctreeCell& cell, int& max_depth);
void find_max_intersecting_leaf_depth( const OctreeCell& cell, int& max_depth );
void collect_finest_intersecting_leaves( const OctreeCell& cell, int target_depth, std::vector<LeafCell>& out );
void collect_intersecting_leaves_verified( const OctreeCell& cell, const Tree& tree, std::vector<LeafCell>& out );

#endif //THESIS_OCTREE_H