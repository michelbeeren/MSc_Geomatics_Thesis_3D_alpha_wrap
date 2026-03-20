//
// Created by Michel Beeren on 20/02/2026.
//

#ifndef THESIS_HAUSDORFF_H
#define THESIS_HAUSDORFF_H

#include <string>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Real_timer.h>
#include <CGAL/alpha_wrap_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
using Vector_3 = K::Vector_3;

namespace PMP = CGAL::Polygon_mesh_processing;

void hausdorff_distance(Mesh& wrapped_mesh, const Mesh& original_mesh, const double filter);
void write_hausdorff_distance(const Mesh& wrapped_mesh, const std::string& filename);
std::vector<double> point_to_mesh_distances(const std::vector<Point_3>& sample_points, const Mesh& target_mesh);
void write_distances_to_csv(const std::string& filename, const std::vector<double>& distances);

#endif //THESIS_HAUSDORFF_H