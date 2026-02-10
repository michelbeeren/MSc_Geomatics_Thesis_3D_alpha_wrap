//
// Created by Michel Beeren on 26/01/2026.
//

#ifndef THESIS_ALPHA_WRAP_H
#define THESIS_ALPHA_WRAP_H
#include "octree.h"
struct MeshData {
    Mesh mesh;
    std::map<Mesh::Face_index, K::Vector_3> face_normals;
    std::unique_ptr<Tree> tree;  // optional
};

// struct LeafCell {CGAL::Bbox_3 bbox;};

Mesh _3D_alpha_wrap(const std::string filename, const double relative_alpha_, const double relative_offset_, MeshData& data_, bool MAT, bool Octree, bool write_out_, bool validate);
double rel_offset_to_offset(Mesh& mesh, const double relative_offset);
Mesh _3D_alpha_inside_wrap(const std::string filename, const double relative_alpha_, const double relative_offset_, Mesh& mesh_, bool write_out_, bool validate);
MeshData mesh_input(const std::string& filename, bool compute_normals = true, bool build_tree = true);

#endif //THESIS_ALPHA_WRAP_H