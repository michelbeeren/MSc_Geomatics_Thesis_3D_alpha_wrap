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
#include <memory>

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
#include "resuls.h"

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

  // exploder_to_off("../data/Input/3DBAG_Buildings/aula.off","../data/Input/exploded/aula_1.off",1,0.0);

  bool beeren_method = true;
  const double max_d_to_input_in_offsets_ = 2;
  bool write_output_ = true;
  bool validate_ = false;
  const bool run_statistics_sweeps = false; // set true to generate sweep CSV files

  Mesh input_mesh_for_statistics;
  bool has_input_mesh_for_statistics = false;
  std::unique_ptr<MeshData> triangle_data;

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
    triangle_data = std::make_unique<MeshData>(mesh_input(filename, true, true)); // set both to false if you do not want to compute normals + tree
    MeshData& data = *triangle_data;
    Mesh& mesh = data.mesh;
    input_mesh_for_statistics = mesh;
    has_input_mesh_for_statistics = true;
    auto face_normals = data.face_normals;
    std::cout << "face_normals.size = " << face_normals.size() << std::endl;
    Tree& tree = *data.tree;
    valid_mesh_boolean(mesh); // is input mesh valid
    Mesh alpha_wrap = _3D_alpha_wrap_tr_mesh(filename, relative_alpha, relative_offset, data,
                                             beeren_method, max_d_to_input_in_offsets_,
                                             write_output_, validate_);
  }
  else
  {
    throw std::runtime_error("Input is neither a valid triangle mesh nor a valid point cloud: " + filename);
  }

  // if (triangle_data) {
  //   Mesh alpha_wrap_octree = _3D_alpha_wrap(filename, relative_alpha, relative_offset, *triangle_data,
  //                                           beeren_method, max_d_to_input_in_offsets_,
  //                                           false, true, write_output_, validate_, false);
  // }

  if (run_statistics_sweeps) {
    if (!has_input_mesh_for_statistics) {
      throw std::runtime_error("Statistics sweeps require a triangle-mesh input.");
    }

    // statistics_over_relative_alpha_to_csv(
    //     {20.0, 30, 40, 50, 60, 70, 80, 90, 100},
    //     relative_offset, max_d_to_input_in_offsets_, input_mesh_for_statistics,
    //     false, false, true, false,
    //     "../data/Output/statistics/sweep_alpha_o1000_n10_20.csv",10);

    // statistics_over_relative_offset_to_csv(
    //     relative_alpha,
    //     {500.0, 1000.0, 2000.0},
    //     max_d_to_input_in_offsets_, input_mesh_for_statistics,
    //     true, true, true, false,
    //     "../data/Output/statistics/sweep_offset.csv");

    // statistics_over_tau_to_csv(
    //     relative_alpha, relative_offset,
    //     {1.1, 1.12, 1.15, 1.2, 1.25, 1.3, 1.5, 1.7, 2.0, 2.5, 3.0, 4, 5, 6, 7.5, 8.5, 10, 12, 15, 20, 25, 30},
    //     input_mesh_for_statistics,
    //     false, true, true, false,
    //     "../data/Output/statistics/sweep_tau_o1000_n10_40.csv", 5);
  }

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
