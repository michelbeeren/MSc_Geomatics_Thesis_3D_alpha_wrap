//
// Created by Michel Beeren on 02/02/2026.
//

#ifndef THESIS_MAT_H
#define THESIS_MAT_H

std::vector<Point_3> _surface_sampling(const Mesh& mesh, const double relative_alpha_);
void write_points_as_off(const std::string& filename, const std::vector<Point_3>& points);
void write_coords_npy(const std::string& coords_path, const std::vector<Point_3>& points);
std::vector<Vector_3> normals_for_points_from_closest_face(const Mesh& mesh, const Tree& tree, const std::map<Mesh::Face_index, Vector_3>& face_normals, const std::vector<Point_3>& points);
void write_normals_npy(const std::string& normals_path, const std::vector<Vector_3>& normals);
void run_masbcpp_compute_ma(const std::string& in_dir, const std::string& out_dir);
void generate_lfs_ply(const std::string& in_dir, const std::string& out_dir, const std::string& output_path);
void generate_lfs_ply_concave(const std::string& in_dir, const std::string& out_dir, const std::string& output_path);

#endif //THESIS_MAT_H