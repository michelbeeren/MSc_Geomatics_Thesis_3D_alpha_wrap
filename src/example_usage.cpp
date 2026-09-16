#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Real_timer.h>

namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;

int main(int argc, char** argv)
{
    // ==================== INPUT FILE PATH ====================
    const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("../example_data/input/3DBAG_simpel_house.obj");

    // ==================== OUTPUT FILE PATH ===================
    const std::string output_name = "../example_data/output/output.obj";

    // ===================== SET PARAMETERS ====================
    const double alpha = 3;
    const double offset = 0.03;
    const double tau = 2*offset;

    // Create mesh from input
    Mesh mesh;
    if (!PMP::IO::read_polygon_mesh(filename, mesh) || CGAL::is_empty(mesh)) {
        std::cerr << "Failed to read input mesh: " << filename << std::endl;
        return EXIT_FAILURE;
    }

    CGAL::Real_timer t;
    t.start();

    // create alpha wrap output mesh
    Mesh wrap;
    CGAL::alpha_wrap_3(mesh, alpha, offset, wrap, CGAL::parameters::max_distance_to_input(tau));
    t.stop();
    std::cout << "🎁 Successfully alpha wrapped! 🎁: " << num_vertices(wrap) << " 🔘vertices🔘, " << num_faces(wrap) << " 📐faces📐, it took " << t.time() << " s.⏰" << std::endl;

    // write output mesh
    std::filesystem::path p(output_name);
    std::filesystem::create_directories(p.parent_path());
    CGAL::IO::write_polygon_mesh(output_name, wrap, CGAL::parameters::stream_precision(20));
    std::cout<<"output mesh written to" << output_name <<std::endl;
}
