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
void write_off_from_mesh(const Mesh& mesh, const std::string& filename);
void write_obj_from_mesh(const Mesh& mesh, const std::string& filename);
bool check_report_for_validity(const std::string& report_path);
bool valid_mesh_boolean(const Mesh& mesh_);

#endif //THESIS_VAL3DITY_H