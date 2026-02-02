//
// Created by Michel Beeren on 02/02/2026.
//

#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/measure.h>   // area()
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <iostream>
namespace PMP = CGAL::Polygon_mesh_processing;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
#include "MAT.h"

// Blue noise sample input mesh
std::vector<Point_3> _surface_sampling(const Mesh& mesh, const double relative_alpha_) {
    std::vector<Point_3> surface_points;

    // compute alpha and offset from a_rel and d_rel and bbox
    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                         CGAL::square(bbox.ymax() - bbox.ymin()) +
                                         CGAL::square(bbox.zmax() - bbox.zmin()));
    // std::cout << "diagonal bbox length: " << diag_length << std::endl;
    const double alpha = diag_length / relative_alpha_;
    const double density = 10.0 / (std::sqrt(3.0) * alpha * alpha);

    PMP::sample_triangle_mesh(
        mesh,
        std::back_inserter(surface_points),
        CGAL::parameters::use_random_uniform_sampling(true)
                        .number_of_points_per_area_unit(density)
                        .do_sample_vertices(false)
                        .do_sample_edges(false)
                        .do_sample_faces(true)
    );

    std::cout << "Sampled " << surface_points.size();

    return surface_points;
}

// For visualization of the blue noise samling
void write_points_as_off(const std::string& filename, const std::vector<Point_3>& points)
{
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Cannot open file " << filename << "\n";
        return;
    }

    out << std::setprecision(17) << std::fixed; // <<<<<< IMPORTANT

    out << "OFF\n";
    out << points.size() << " 0 0\n";

    for (const auto& p : points) {
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
}




