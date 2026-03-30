//
// Created by Michel Beeren on 26/01/2026.
//

#ifndef THESIS_ALPHA_WRAP_H
#define THESIS_ALPHA_WRAP_H
#include "octree.h"
using Point_container = std::vector<Point_3>;
struct MeshData {
    Mesh mesh;  // CGAL Surface_mesh type
    std::map<Mesh::Face_index, K::Vector_3> face_normals;  // Normals of the faces indexed by Face_index
    std::unique_ptr<Tree> tree;  // Optional octree structure
    std::map<Mesh::Face_index, std::set<Mesh::Face_index>> adjacency_map;  // Adjacency map of faces indexed by Face_index
};
enum class Input_kind
{
    Triangle_mesh,
    Point_cloud,
    Unknown
};

// struct LeafCell {CGAL::Bbox_3 bbox;};

Mesh _3D_alpha_wrap(const std::string filename, const double relative_alpha_, const double relative_offset_, MeshData& data_, bool max_d_in_offets, const double max_d, bool MAT, bool Octree, bool write_out_, bool validate, bool hausdorff);
Mesh _3D_alpha_wrap_tr_mesh(const std::string filename, const double relative_alpha_, const double relative_offset_, MeshData& data_, bool max_d_in_offets, const double max_d, bool write_out_, bool validate);
Mesh _3D_alpha_wrap_pc(const std::string filename, const double relative_alpha_, const double relative_offset_, bool max_d_in_offets, const double max_d, bool write_out_, bool validate);
double rel_offset_to_offset(Mesh& mesh, const double relative_offset);
Mesh _3D_alpha_inside_wrap(const std::string filename, const double relative_alpha_, const double relative_offset_, Mesh& mesh_, bool write_out_, bool validate);
MeshData mesh_input(const std::string& filename, bool compute_normals = true, bool build_tree = true);
Point_container point_input(const std::string& filename);
Input_kind detect_input_kind(const std::string& filename);

#endif //THESIS_ALPHA_WRAP_H