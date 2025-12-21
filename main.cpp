#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_wrap_3.h>
#include <CGAL/IO/facets_in_complex_3_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>

int main() {
    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point = K::Point_3;
    using Mesh = CGAL::Surface_mesh<Point>;

    // Load a surface mesh from a file
    std::ifstream input("../Test_3D_Alphawrap/Input/cube2.off");
    Mesh input_mesh;
    if (!input || !(input >> input_mesh) || input_mesh.is_empty()) {
        std::cerr << "Invalid or missing input mesh file." << std::endl;
        return EXIT_FAILURE;
    }

    // Parameters
    const double relative_alpha = 0.4; // Controls wrapping tightness
    const double offset = 0.1;           // Distance around geometry

    Mesh wrap;
    CGAL::alpha_wrap_3(input_mesh, relative_alpha, offset, wrap);

    std::ofstream output("../Test_3D_Alphawrap/Output/cube_wrapped2.off");
    output << wrap;
    std::cout << "Wrapped mesh saved to wrapped.off" << std::endl;

    return EXIT_SUCCESS;
}
