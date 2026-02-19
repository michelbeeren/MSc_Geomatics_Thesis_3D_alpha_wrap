// ========================== INCLUDE STUFF ===============================
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
#include <filesystem>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/IO/Color.h>
#include <fstream>
#include <vector>

// ========================= NAMESPACES/USING =============================
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

// ============================= STRUCTS ===================================
// struct OctreeCell {
//     CGAL::Bbox_3 bbox;
//     int depth = 0;
//     std::vector<Mesh::Face_index> faces;
//     std::array<std::unique_ptr<OctreeCell>, 8> children;
//
//     OctreeCell() {
//         for (auto& c : children)
//             c = nullptr;
//     }
// };
// struct LeafCell {
//     CGAL::Bbox_3 bbox;
// };

// ====================== OWN CREATED HEADERS ==============================
#include "val3dity.h"
#include "alpha_wrap.h"
#include "edge_refinement.h"
#include "octree.h"
#include "MAT.h"

// ==================================================================================
// ==================================== MAIN ========================================
// ==================================================================================
int main(int argc, char** argv)
{
  // Read the input
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("../data/Input/3DBAG_Buildings/bouwkunde.off");
    std::cout << "------------------------------------------------------------" << std::endl;
  std::cout << "Reading input: " << filename << std::endl;

  const double relative_alpha = 2000; //2000. //20. //1000.
  const double relative_offset = 10000.; // 7000. //600. //12000.

    // ----------------------MESH INPUT FILE (optional: compute normals and tree)------------------
    auto data = mesh_input(filename, true, true); // set both to false if you do not want to compute normals + tree
    Mesh& mesh = data.mesh;
    auto face_normals = data.face_normals;
    std::cout << "face_normals.size = " << face_normals.size() << std::endl;
    Tree& tree = *data.tree;

    // ----------------------------IS INPUT MESH VALID?------------------------------------
    valid_mesh_boolean(mesh);


    // ------------------------------ALPHA WRAP INPUT---------------------------------------
    Mesh alpha_wrap = _3D_alpha_wrap(filename,relative_alpha,relative_offset, data, false, true, true, false); // set both to false if you do not want to write out the file and test if valid
    // -------------------------- alpha wrap from inside -----------------------------------
    // Mesh alpha_inside_wrap = _3D_alpha_inside_wrap( filename,relative_alpha,relative_offset, mesh, true, false);

      // // -------------------------------OCTREE REFINEMENT---------------------------
    // // Create root octree cell
    // OctreeCell root;
    // // compute bbox
    // CGAL::Bbox_3 tight = CGAL::Polygon_mesh_processing::bbox(mesh);
    //
    // // Center
    // double cx = 0.5 * (tight.xmin() + tight.xmax());
    // double cy = 0.5 * (tight.ymin() + tight.ymax());
    // double cz = 0.5 * (tight.zmin() + tight.zmax());
    //
    // // Largest half-extent
    // double half = std::max({
    //     0.5 * (tight.xmax() - tight.xmin()),
    //     0.5 * (tight.ymax() - tight.ymin()),
    //     0.5 * (tight.zmax() - tight.zmin())
    // });
    //
    // // Optional looseness
    // half *= 1.1;
    //
    // // TRUE cube
    // root.bbox = CGAL::Bbox_3(
    //     cx - half, cy - half, cz - half,
    //     cx + half, cy + half, cz + half
    // );
    // root.depth = 0;
    //
    // // Run Octree refinement
    // refine_cell(root, tree, face_normals);
    //
    // // Collect all leaf cells
    // std::vector<LeafCell> leaves_all; // use these if you want all cubes
    // collect_leaves(root, leaves_all);
    //
    // // Collect leaf cells that intersect the input
    // std::vector<LeafCell> leaves;
    // collect_intersecting_leaves(root, leaves);
    //
    // // Collect finest level of detail leaves
    // int finest_depth = 0;
    // find_max_intersecting_leaf_depth(root, finest_depth);
    //
    // std::vector<LeafCell> cubes;
    // collect_finest_intersecting_leaves(root, finest_depth, cubes);
    //
    // // collect intersecting cubes 2.0
    // std::vector<LeafCell> cubes2;
    // collect_intersecting_leaves_verified(root, tree, cubes2);
    //
    // std::cout << "Leaf cells: " << cubes2.size() << std::endl; // can also use 'leaves_all' or 'leaves' or 'cubes'
    //
    // // Write .off wireframe file
    // std::string octree_refinement = "../data/Output/3DBAG_Buildings/Octree_refinement.off";
    //
    // write_off_wireframe(cubes2, octree_refinement); // can also use 'leaves_all' or 'leaves' or 'cubes'
    //
    // std::cout << "Octree wireframe written to:\n"
    //           << octree_refinement << std::endl;

    //     // ----------------- SHARPENING EDGES --------------------
    // std::string weight_output_wrap =
    //     "../data/Output/demo/refined.ply";
    // // std::string weight_output_inner_wrap =
    // //     "../data/Test_3D_Alphawrap/Output/3DBAG_Buildings/d_innerwrap_input_bk.ply";
    //
    // // // offset wrapped from inside mesh
    // double _2_offset = rel_offset_to_offset(mesh, relative_offset);
    // double offset = (rel_offset_to_offset(mesh, relative_offset))/2;
    // // Mesh offset_inner_wrap = offset_mesh(alpha_inside_wrap, _2_offset);
    // // std::cout << "Succesfully offsetted wrapped from inside mesh" << std::endl;
    //
    // Mesh offset_input = offset_mesh(mesh, offset);
    //
    // // tag + colorise outer wrap
    // CGAL::Real_timer t;
    // t.start();
    // Mesh distance_alpha = add_midpoint_distance_tag(alpha_wrap, offset_input, 220);
    // // std::cout << "Succesfully adding midpoint distance to alpha wrap mesh" << std::endl;
    //
    // // std::cout << "Writing to " << weight_output_wrap << std::endl;
    // {
    //   std::ofstream out(weight_output_wrap, std::ios::binary);
    //   if (!out) {
    //       std::cerr << "Cannot open " << weight_output_wrap << " for writing\n";
    //   } else {
    //       CGAL::IO::write_PLY(
    //           out,
    //           distance_alpha,
    //           CGAL::parameters::stream_precision(17)
    //       );
    //   }
    // }
    // std::string refining =
    // "../data/Output/demo/sharpe_concave_corners.ply";
    // Mesh refined_edges = refine_round_edges(distance_alpha, offset_input);
    // t.stop();
    // std::cout << "Took " << t.time() << " s." << std::endl;
    // std::cout << "Succesfully refined edges" << std::endl;
    //
    // std::cout << "Writing to " << refining << std::endl;
    // {
    //   std::ofstream out(refining, std::ios::binary);
    //   if (!out) {
    //       std::cerr << "Cannot open " << refining << " for writing\n";
    //   } else {
    //       CGAL::IO::write_PLY(
    //           out,
    //           refined_edges,
    //           CGAL::parameters::stream_precision(17)
    //       );
    //   }
    // }
    // // ----------------------------IS MESH VALID?------------------------------------
    // // std::cout << "Testing if input mesh is valid:" << std::endl;
    // valid_mesh_boolean(refined_edges);
    // std::cout << "------------------------------------------------------------" << std::endl;



    // // tag + colorise inner wrap
    // Mesh distance_alpha_inner = add_midpoint_distance_tag(offset_inner_wrap, offset_input, 220);
    // std::cout << "Succesfully adding midpoint distance to alpha wrap from inside mesh" << std::endl;
    //
    // std::cout << "Writing to " << weight_output_inner_wrap << std::endl;
    // {
    //   std::ofstream out(weight_output_inner_wrap, std::ios::binary);
    //   if (!out) {
    //       std::cerr << "Cannot open " << weight_output_inner_wrap << " for writing\n";
    //   } else {
    //       CGAL::IO::write_PLY(
    //           out,
    //           distance_alpha_inner,
    //           CGAL::parameters::stream_precision(17)
    //       );
    //   }
    // }

  return EXIT_SUCCESS;
}