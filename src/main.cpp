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

// ====================== OWN CREATED HEADERS ==============================
#include "val3dity.h"
#include "alpha_wrap.h"
#include "edge_refinement.h"
#include "octree.h"
#include "MAT.h"
#include "hausdorff.h"

// ==================================================================================
// ==================================== MAIN ========================================
// ==================================================================================
int main(int argc, char** argv)
{
  // Read the input
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("../data/Input/3DBAG_Buildings/aula.off");
    std::cout << "------------------------------------------------------------" << std::endl;
  std::cout << "Reading input: " << filename << std::endl;

  const double relative_alpha = 20; //2000. //20. //1000.
  const double relative_offset = 6000.; // 7000. //600. //12000.

  bool beeren_method = true;
  const double max_d_to_input_in_offsets_ = 1.5;
  bool write_output_ = true;
  bool validate_ = true;

  const Input_kind kind = detect_input_kind(filename);

  if (kind == Input_kind::Point_cloud) {
    std::cout << "⚽️🏀⚾️🥎🎾🏐🎱 input type = Point Cloud ⚽️🏀⚾️🥎🎾🏐🎱" << std::endl;
    if (beeren_method) {
      std::cout << "🐻🐻🐻 BEEREN METHOD 🐻🐻🐻" << std::endl;
      std::cout << "Ⓜ️🅰️❎ 📐distance📐 to input in offsets = " << max_d_to_input_in_offsets_ << std::endl;
    }
    else if (!beeren_method) {
      std::cout << "🏃🏼‍♀️‍➡️🏃🏽‍♀️‍➡️🏃🏾‍♀️‍➡️🏃🏿‍♀️‍➡️Running normal algorithm" << std::endl;
    }
    Mesh alpha_wrap_pc = _3D_alpha_wrap_pc(filename,relative_alpha,relative_offset,beeren_method,max_d_to_input_in_offsets_,write_output_,validate_);
  }
  else if (kind == Input_kind::Triangle_mesh) {
    std::cout << "🔺⚠️🔼▲ input type = triangle mesh ▲🔼⚠️🔺" << std::endl;
    if (beeren_method) {
      std::cout << "🐻🐻🐻 BEEREN METHOD 🐻🐻🐻" << std::endl;
      std::cout << "Ⓜ️🅰️❎ 📐distance📐 to input in offsets = " << max_d_to_input_in_offsets_ << std::endl;
    }
    else if (!beeren_method) {
      std::cout << "🏃🏼‍♀️‍➡️🏃🏽‍♀️‍➡️🏃🏾‍♀️‍➡️🏃🏿‍♀️‍➡️Running normal algorithm" << std::endl;
    }
    auto data = mesh_input(filename, true, true); // set both to false if you do not want to compute normals + tree
    Mesh& mesh = data.mesh;
    auto face_normals = data.face_normals;
    std::cout << "face_normals.size = " << face_normals.size() << std::endl;
    Tree& tree = *data.tree;
    valid_mesh_boolean(mesh); // is input mesh valid
    Mesh alpha_wrap = _3D_alpha_wrap_tr_mesh(filename,relative_alpha,relative_offset,data, beeren_method, max_d_to_input_in_offsets_, write_output_, validate_);
  }
  else
  {
    throw std::runtime_error("Input is neither a valid triangle mesh nor a valid point cloud: " + filename);
  }

    // ------------------------------ALPHA WRAP INPUT---------------------------------------
    // Mesh alpha_wrap = _3D_alpha_wrap(filename,relative_alpha,relative_offset, data, true, false, false, true, true, false); // set both to false if you do not want to write out the file and test if valid


    // ---------------------------------STATISTICS---------------------------------------
    // std::vector<Point_3> samples = _surface_sampling(mesh, 200.0);
    // std::vector<double> distances = point_to_mesh_distances(samples, alpha_wrap);
    // std::cout << "distances.size() = " << distances.size() << std::endl;
    // std::cout << "samples.size() = " << samples.size() << std::endl;
    // write_distances_to_csv("../data/Output/stats/distances_bk_normal.csv", distances);

    // -------------------------- alpha wrap from inside -----------------------------------
    // Mesh alpha_inside_wrap = _3D_alpha_inside_wrap( filename,relative_alpha,relative_offset, mesh, true, false);

  return EXIT_SUCCESS;
}