//
// Created by Michel Beeren on 26/01/2026.
//

#include <string>
#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <vector>
#include <iomanip>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Real_timer.h>
#include <CGAL/alpha_wrap_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;

#include "alpha_wrap.h"

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
Mesh _3D_alpha_wrap(const std::string output_, const double relative_alpha_, const double relative_offset_, Mesh& mesh_) {
    // std::cout << "omg omg it is a triangle mesh!;" << std::endl;
    // std::cout << "relative_alpha: " << relative_alpha_ << " and relative offset: " << relative_offset_ << std::endl;
    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh_);
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                         CGAL::square(bbox.ymax() - bbox.ymin()) +
                                         CGAL::square(bbox.zmax() - bbox.zmin()));
    // std::cout << "diagonal bbox length: " << diag_length << std::endl;
    const double alpha = diag_length / relative_alpha_;
    const double offset = diag_length / relative_offset_;
    std::cout << "Alpha wrapping with a_rel = " << relative_alpha_ << " and offset_rel = " << relative_offset_ << ", which results in alpha = " << alpha << " and offset = " << offset << std::endl;
    // Construct the wrap
    CGAL::Real_timer t;
    t.start();

    Mesh wrap;
    CGAL::alpha_wrap_3(mesh_, alpha, offset, wrap);

    t.stop();
    std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces, it took " << t.time() << " s." << std::endl;
    // std::cout << "Took " << t.time() << " s." << std::endl;

    // ------------------------------ write output mesh -------------------------------
    const std::string output_name = output_;
    std::cout << "Writing to " << output_name << std::endl;
    CGAL::IO::write_polygon_mesh(output_name, wrap, CGAL::parameters::stream_precision(17));

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

// 3D alpha wrap from inside TODO point inside generation can be better
Mesh _3D_alpha_inside_wrap(const std::string output_, const double relative_alpha_, const double relative_offset_, Mesh& mesh_) {
  // std::cout << "omg omg it is a triangle mesh!;" << std::endl;
  // std::cout << "relative_alpha: " << relative_alpha_ << " and relative offset: " << relative_offset_ << std::endl;
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh_);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));
  // std::cout << "diagonal bbox length: " << diag_length << std::endl;
  const double alpha = diag_length / relative_alpha_;
  const double offset = diag_length / relative_offset_;
    std::cout << "Alpha wrapping from inside with a_rel = " << relative_alpha_ << " and offset_rel = " << relative_offset_ << ", which results in alpha = " << alpha << " and offset = " << offset << std::endl;

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

    // // point now placed in middle of bounding box
    // std::vector<Point_3> seeds = {
    //     Point_3((bbox.xmin() + bbox.xmax()) / 2,(bbox.ymin() + bbox.ymax()) / 2,(bbox.zmin() + bbox.zmax()) / 2)
    // };

    // with random point generated inside input mesh -> not reliable method still
    std::vector<Point_3> seeds = {
        random_point_inside_mesh(mesh_)
    };

  Mesh wrap;
  CGAL::alpha_wrap_3(mesh_, alpha, offset, wrap, CGAL::parameters::seed_points(std::ref(seeds)));

  t.stop();
  std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces, it took " << t.time() << " s." << std::endl;
  // std::cout << "Took " << t.time() << " s." << std::endl;

    // ------------------------------ write out output file -------------------------------
    std::string out_path = output_;
    // Remove last 4 characters (".off")
    if (out_path.size() > 4 && out_path.substr(out_path.size() - 4) == ".off") {
        out_path = out_path.substr(0, out_path.size() - 4);
    }

    // Add cuurent alpha and offset to the output name
    std::string out_name = out_path + "_inside_wrap.off";

    std::cout << "Writing to " << out_name << std::endl;
    CGAL::IO::write_polygon_mesh(out_name, wrap, CGAL::parameters::stream_precision(17));

    // std::string offset_output = "../data/Test_3D_Alphawrap/Output/3DBAG_Buildings/offset_testtest.off";
    // Mesh offset_mesh = offset_mesh_by_5cm(wrap, offset);

    // std::cout << "Writing to " << offset_output << std::endl;
    // CGAL::IO::write_polygon_mesh(offset_output, offset_mesh, CGAL::parameters::stream_precision(17));

    return wrap;
}
