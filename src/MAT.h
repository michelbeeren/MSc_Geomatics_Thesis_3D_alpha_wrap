//
// Created by Michel Beeren on 02/02/2026.
//

#ifndef THESIS_MAT_H
#define THESIS_MAT_H

std::vector<Point_3> _surface_sampling(const Mesh& mesh, const double relative_alpha_);
void write_points_as_off(const std::string& filename, const std::vector<Point_3>& points);
void write_points_as_npy(const std::string& filename, const std::vector<Point_3>& points);
void write_coords_npy(const std::string& coords_path, const std::vector<Point_3>& points);
std::vector<Vector_3> normals_for_points_from_closest_face(const Mesh& mesh, const Tree& tree, const std::map<Mesh::Face_index, Vector_3>& face_normals, const std::vector<Point_3>& points);
void write_normals_npy(const std::string& normals_path, const std::vector<Vector_3>& normals);
void run_masbcpp_compute_ma(const std::string& in_dir, const std::string& out_dir);
void write_npy_points_as_off(const std::string& npy_file, const std::string& off_file);
// void write_mat_colored_coff(const std::string& coords_npy, const std::string& radii_npy, const std::string& coff_file);
void write_mat_colored_coff(const std::string& ma_coords_out_npy,const std::string& ma_qidx_out_npy,const std::string& coords_indexed_by_qidx_npy,const std::string& coff_file);
void write_surface_colored_by_ball_proxy(const std::string& ma_coords_out_npy, const std::string& ma_qidx_out_npy, const std::string& coords_indexed_by_qidx_npy, const std::string& coff_file);
void write_normals_as_obj_lines(const std::string& obj_file,
                                const std::vector<Point_3>& pts,
                                const std::vector<Vector_3>& normals,
                                double scale);
double bbox_diag(const Mesh& mesh);
void write_surface_colored_by_power_distance(const std::string& ma_coords_out_npy,
                                             const std::string& ma_qidx_out_npy,
                                             const std::string& surface_coords_npy,
                                             const std::string& coff_file);
void generate_lfs_ply(const std::string& in_dir,
                      const std::string& out_dir,
                      const std::string& output_path);

#endif //THESIS_MAT_H