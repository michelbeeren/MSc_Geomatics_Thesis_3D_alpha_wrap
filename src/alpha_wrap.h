//
// Created by Michel Beeren on 26/01/2026.
//

#ifndef THESIS_ALPHA_WRAP_H
#define THESIS_ALPHA_WRAP_H

std::string generate_output_name(const std::string input_name_, const double rel_alpha_, const double rel_offset_);
Mesh _3D_alpha_wrap(const std::string output_, const double relative_alpha_, const double relative_offset_, Mesh& mesh_);
double rel_offset_to_offset(Mesh& mesh, const double relative_offset);
Point_3 random_point_inside_mesh(const Mesh& mesh);
Mesh _3D_alpha_inside_wrap(const std::string output_, const double relative_alpha_, const double relative_offset_, Mesh& mesh_);

#endif //THESIS_ALPHA_WRAP_H