//
// Created by Michel Beeren on 04/05/2026.
//

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Real_timer.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <iostream>
#include <string>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
using face_descriptor = Mesh::Face_index;
using Ray_3       = K::Ray_3;
using Segment_3   = K::Segment_3;
using Primitive   = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using AABB_traits = CGAL::AABB_traits<K, Primitive>;
using Tree        = CGAL::AABB_tree<AABB_traits>;

struct Arguments
{
    std::string input_path;
    std::string output_path;
    double alpha = 0.0;
    double offset = 0.0;
    double tau = 0.0;
};

// ToDo what input types are possible???
void print_usage(const char* modified_3d_alpha_wrapping)
{
    std::cerr
      << "Usage:\n"
      << "  " << modified_3d_alpha_wrapping
      << " <input> <output_mesh> --alpha <value> --offset <value> --tau <value>\n\n"
      << "Example:\n"
      << "  " << modified_3d_alpha_wrapping
      << " input.off output.off --alpha 0.03 --offset 0.005 --tau 0.5\n";
}

bool parse_arguments(int argc, char** argv, Arguments& args)
{
    if(argc < 9)
    {
        return false;
    }

    args.input_path = argv[1];
    args.output_path = argv[2];

    for(int i = 3; i < argc; ++i)
    {
        const std::string key = argv[i];

        if(key == "--alpha" && i + 1 < argc)
        {
            args.alpha = std::atof(argv[++i]);
        }
        else if(key == "--offset" && i + 1 < argc)
        {
            args.offset = std::atof(argv[++i]);
        }
        else if(key == "--tau" && i + 1 < argc)
        {
            args.tau = std::atof(argv[++i]);
        }
        else
        {
            std::cerr << "Unknown or incomplete argument: " << key << "\n";
            return false;
        }
    }

    if(args.alpha <= 0.0)
    {
        std::cerr << "Error: alpha must be positive.\n";
        return false;
    }

    if(args.offset <= 0.0)
    {
        std::cerr << "Error: offset must be positive.\n";
        return false;
    }

    if(args.tau <= 1.0)
    {
        std::cerr << "Error: tau must be higher than 1.\n";
        return false;
    }

    return true;
}


int main(int argc, char** argv)
{
    Arguments args;

    if(!parse_arguments(argc, argv, args))
    {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    Mesh input_mesh;

    std::cout << "Reading input mesh: " << args.input_path << "\n";

    if(!CGAL::IO::read_polygon_mesh(args.input_path, input_mesh))
    {
        std::cerr << "Error: could not read input mesh.\n";
        return EXIT_FAILURE;
    }

    if(input_mesh.is_empty())
    {
        std::cerr << "Error: input mesh is empty.\n";
        return EXIT_FAILURE;
    }

    Mesh wrap_mesh;

    std::cout << "Running modified alpha wrap...\n";
    std::cout << "  alpha  = " << args.alpha << "\n";
    std::cout << "  offset = " << args.offset << "\n";
    std::cout << "  tau    = " << args.tau << "\n";

    // Replace this call with the exact signature of your modified function.
    //
    // Option A: if you added a tau overload:
    CGAL::alpha_wrap_3(input_mesh, args.alpha, args.offset, args.tau, wrap_mesh);

    // Option B: if your modified algorithm still uses the original CGAL signature,
    // then use this instead:
    //
    // CGAL::alpha_wrap_3(input_mesh, args.alpha, args.offset, wrap_mesh);

    std::cout << "Writing output mesh: " << args.output_path << "\n";

    if(!CGAL::IO::write_polygon_mesh(args.output_path, wrap_mesh,
                                     CGAL::parameters::stream_precision(17)))
    {
        std::cerr << "Error: could not write output mesh.\n";
        return EXIT_FAILURE;
    }

    std::cout << "Done.\n";

    return EXIT_SUCCESS;
}