//
// Created by Michel Beeren on 02/02/2026.
//

#ifndef THESIS_MAT_H
#define THESIS_MAT_H

std::vector<Point_3> _surface_sampling(const Mesh& mesh, const double relative_alpha_);
void write_points_as_off(const std::string& filename, const std::vector<Point_3>& points);

#endif //THESIS_MAT_H