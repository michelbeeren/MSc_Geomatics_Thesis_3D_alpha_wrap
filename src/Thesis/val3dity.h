//
// Created by Michel Beeren on 26/01/2026.
//


#ifndef THESIS_VAL3DITY_H
#define THESIS_VAL3DITY_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;

void write_cityjson_from_mesh(const Mesh& mesh, const std::string& filename);
bool report_is_valid(const std::string& report_path);
bool run_val3dity_and_check(const std::string& input_path, const std::string& report_path, const std::string& primitive = "");
bool valid_mesh_boolean(const Mesh& mesh_);
bool valid_file_boolean(const std::string& input_path);

#endif //THESIS_VAL3DITY_H