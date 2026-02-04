//
// Created by Michel Beeren on 26/01/2026.
//

#include <string>
#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <vector>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Real_timer.h>
#include <CGAL/alpha_wrap_3.h>
#include <map>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/variant/get.hpp>
#include <CGAL/IO/Color.h>

namespace PMP = CGAL::Polygon_mesh_processing;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;
using Segment_3 = K::Segment_3;
using Primitive   = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using AABB_traits = CGAL::AABB_traits<K, Primitive>;
using Tree        = CGAL::AABB_tree<AABB_traits>;

#include "alpha_wrap.h"
#include "val3dity.h"

// --------------------------------------MESH INPUT-------------------------------------
MeshData mesh_input(const std::string& filename, bool compute_normals, bool build_tree)
{
    // meshes input file, optionally also computes normals and builds a tree
    MeshData out;

    if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, out.mesh) ||
        CGAL::is_empty(out.mesh) ||
        !CGAL::is_triangle_mesh(out.mesh))
    {
        throw std::runtime_error("Failed to read mesh: " + filename);
    }

    if (compute_normals) {
        CGAL::Polygon_mesh_processing::compute_face_normals(
            out.mesh,
            boost::make_assoc_property_map(out.face_normals)
        );
    }

    if (build_tree) {
        out.tree = std::make_unique<Tree>(
            faces(out.mesh).begin(),
            faces(out.mesh).end(),
            out.mesh
        );
        out.tree->accelerate_distance_queries();

        // defensive check (should never trigger, but matches your wish)
        if (!out.tree) {
            throw std::runtime_error("Tree not built for: " + filename);
        }
    }

    std::cout << "Input file successfully meshed" << std::endl;
    return out;
}

// --------------------------------------ALPHA WRAP-------------------------------------
// generate output name
std::string generate_output_name(const std::string input_name_, const double rel_alpha_, const double rel_offset_)
{
    const int rel_alpha_for_txt = rel_alpha_;
    const int rel_offset_for_txt = rel_offset_;
    const std::string rel_alpha_txt = std::to_string(rel_alpha_for_txt);
    const std::string rel_offset_txt = std::to_string(rel_offset_for_txt);

    // Replace "Input" → "Output"
    std::string out_path = input_name_;
    std::string from = "Input";
    std::string to   = "Output";

    size_t pos = out_path.find(from);
    if (pos != std::string::npos) {
        out_path.replace(pos, from.length(), to);
    }

    // Remove last 4 characters (".off")
    if (out_path.size() > 4 && out_path.substr(out_path.size() - 4) == ".off") {
        out_path = out_path.substr(0, out_path.size() - 4);
    }

    // Add cuurent alpha and offset to the output name
    std::string out_name = out_path + "_a=" + rel_alpha_txt + "_offset=" + rel_offset_txt + ".off";
    return out_name;
}

// 3D alpha wrap
Mesh _3D_alpha_wrap(const std::string filename, const double relative_alpha_, const double relative_offset_, Mesh& mesh_, bool write_out_, bool validate) {
    // compute alpha and offset from a_rel and d_rel and bbox
    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh_);
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                         CGAL::square(bbox.ymax() - bbox.ymin()) +
                                         CGAL::square(bbox.zmax() - bbox.zmin()));
    // std::cout << "diagonal bbox length: " << diag_length << std::endl;
    const double alpha = diag_length / relative_alpha_;
    const double offset = diag_length / relative_offset_;
    std::cout << "--------------------3D ALPHA WRAPPING THE INPUT:----------------" << std::endl;
    std::cout << "alpha = " << alpha << " (a_rel = " << relative_alpha_ << ") and offset = " << offset << " (offset_rel = " << relative_offset_ << ")" << std::endl;
    // Construct the wrap
    CGAL::Real_timer t;
    t.start();

    Mesh wrap;
    CGAL::alpha_wrap_3(mesh_, alpha, offset, wrap);

    t.stop();
    std::cout << "🎁 succesfully alpha wrapped! 🎁: " << num_vertices(wrap) << " 🔘vertices🔘, " << num_faces(wrap) << " 📐faces📐, it took " << t.time() << " s.⏰" << std::endl;
    // std::cout << "Took " << t.time() << " s." << std::endl;

    // ------------------------------ write output mesh -------------------------------
    if (write_out_) {
        // generate output name from input name
        std::string output_ = generate_output_name(filename, relative_alpha_, relative_offset_);
        std::filesystem::path p(output_);
        std::filesystem::create_directory(p.parent_path());
    std::cout << "📝Writing📝 to: " << output_ << std::endl;
    CGAL::IO::write_polygon_mesh(output_, wrap, CGAL::parameters::stream_precision(25));}

    // ----------------------------- validate the output ------------------------------
    if (validate) {
        valid_mesh_boolean(wrap);
    }
    std::cout << "---------------------------------------------------------------" << std::endl;

    return wrap;

}

// compute offset
double rel_offset_to_offset(Mesh& mesh, const double relative_offset)
{
    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                         CGAL::square(bbox.ymax() - bbox.ymin()) +
                                         CGAL::square(bbox.zmax() - bbox.zmin()));
    // std::cout << "diagonal bbox length: " << diag_length << std::endl;
    const double offset = diag_length / relative_offset;
    return offset;
}

// -------------------------------------WRAP FROM INSIDE------------------------------------
// TODO not sure if this is really reliable point placer if input mesh is not valid
Point_3 random_point_inside_mesh(const Mesh& mesh)
{
    // containment oracle
    CGAL::Side_of_triangle_mesh<Mesh, K> inside(mesh);

    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
    CGAL::Random rng;

    while (true)
    {
        double x = rng.get_double(bbox.xmin(), bbox.xmax());
        double y = rng.get_double(bbox.ymin(), bbox.ymax());
        double z = rng.get_double(bbox.zmin(), bbox.zmax());
        Point_3 p(x, y, z);

        if (inside(p) == CGAL::ON_BOUNDED_SIDE)
            return p;
    }
}

Mesh _3D_alpha_inside_wrap(const std::string filename, const double relative_alpha_, const double relative_offset_, Mesh& mesh_, bool write_out_, bool validate) {
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh_);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));
  const double alpha = diag_length / relative_alpha_;
  const double offset = diag_length / relative_offset_;

    std::cout << "-------------3D ALPHA WRAPPING THE INPUT FROM INSIDE:----------" << std::endl;
    std::cout << "alpha = " << alpha << " (a_rel = " << relative_alpha_ << ") and offset = " << offset << " (offset_rel = " << relative_offset_ << ")" << std::endl;

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

    std::vector<Point_3> seeds = {
        random_point_inside_mesh(mesh_)
    };

  Mesh wrap;
  CGAL::alpha_wrap_3(mesh_, alpha, offset, wrap, CGAL::parameters::seed_points(std::ref(seeds)));

  t.stop();
    std::cout << "🎁 succesfully alpha wrapped! 🎁: " << num_vertices(wrap) << " 🔘vertices🔘, " << num_faces(wrap) << " 📐faces📐, it took " << t.time() << " s.⏰" << std::endl;

    if (write_out_) {
        // ---------------------------- generate output name -------------------------------
        std::string output_ = generate_output_name(filename, relative_alpha_, relative_offset_);
        std::filesystem::path p(output_);
        std::filesystem::create_directory(p.parent_path());

        // Remove last 4 characters (".off")
        if (output_.size() > 4 && output_.substr(output_.size() - 4) == ".off") {
            output_ = output_.substr(0, output_.size() - 4);
        }

        // Add current alpha and offset to the output name
        std::string out_name = output_ + "_inside_wrap.off";

        // ------------------------------ write out output file -------------------------------
        std::cout << "📝Writing📝 to: " << output_ << std::endl;
        CGAL::IO::write_polygon_mesh(out_name, wrap, CGAL::parameters::stream_precision(17));
    }

    // ------------------------------ validate wrapped mesh ------------------------------------
    if (validate) {
        valid_mesh_boolean(wrap);
    }

    std::cout << "---------------------------------------------------------------" << std::endl;
    return wrap;
}
