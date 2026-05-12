#ifndef THESIS_CLI_ALPHA_WRAPPING_H
#define THESIS_CLI_ALPHA_WRAPPING_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/boost/graph/helpers.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cli_helpers
{
namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
using Point_container = std::vector<Point_3>;

enum class Input_kind
{
  Triangle_mesh,
  Point_cloud,
  Unknown
};

struct Wrap_request
{
  std::string input_path;
  double relative_alpha = 0.0;
  double relative_offset = 0.0;
  double tau = 0.0;
  bool use_noprmal_alpha_wrap = false;
};

struct Wrap_result
{
  Input_kind input_kind = Input_kind::Unknown;
  Mesh wrap_mesh;
  double absolute_alpha = 0.0;
  double absolute_offset = 0.0;
  std::size_t input_size = 0;
};

inline bool is_positive_finite(const double x)
{
  return std::isfinite(x) && x > 0.0;
}

inline double bbox_diagonal(const CGAL::Bbox_3& bbox)
{
  const double dx = bbox.xmax() - bbox.xmin();
  const double dy = bbox.ymax() - bbox.ymin();
  const double dz = bbox.zmax() - bbox.zmin();
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

inline Mesh read_mesh_or_throw(const std::string& path, bool* was_triangulated = nullptr)
{
  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(path, mesh) || CGAL::is_empty(mesh))
    throw std::runtime_error("Failed to read a polygon mesh from: " + path);

  bool triangulated = false;
  if(!CGAL::is_triangle_mesh(mesh))
  {
    triangulated = true;
    PMP::triangulate_faces(mesh);
    if(!CGAL::is_triangle_mesh(mesh))
      throw std::runtime_error("Failed to triangulate mesh from: " + path);
  }

  if(was_triangulated)
    *was_triangulated = triangulated;

  return mesh;
}

inline Point_container read_points_or_throw(const std::string& path)
{
  Point_container points;
  if(!CGAL::IO::read_points(path, std::back_inserter(points)) || points.empty())
    throw std::runtime_error("Failed to read a point cloud from: " + path);

  return points;
}

inline Input_kind detect_input_kind(const std::string& filename)
{
  {
    Mesh mesh;
    if(PMP::IO::read_polygon_mesh(filename, mesh) && !CGAL::is_empty(mesh))
      return Input_kind::Triangle_mesh;
  }

  {
    Point_container points;
    if(CGAL::IO::read_points(filename, std::back_inserter(points)) && !points.empty())
      return Input_kind::Point_cloud;
  }

  return Input_kind::Unknown;
}

inline std::pair<double, double> compute_absolute_alpha_offset(const Mesh& mesh,
                                                               const double relative_alpha,
                                                               const double relative_offset)
{
  const CGAL::Bbox_3 bbox = PMP::bbox(mesh);
  const double diag = bbox_diagonal(bbox);
  return std::make_pair(diag / relative_alpha, diag / relative_offset);
}

inline std::pair<double, double> compute_absolute_alpha_offset(const Point_container& points,
                                                               const double relative_alpha,
                                                               const double relative_offset)
{
  const CGAL::Bbox_3 bbox = CGAL::bbox_3(points.begin(), points.end());
  const double diag = bbox_diagonal(bbox);
  return std::make_pair(diag / relative_alpha, diag / relative_offset);
}

inline Mesh wrap_triangle_mesh(const Mesh& input_mesh,
                               const double alpha,
                               const double offset,
                               const double tau,
                               const bool use_noprmal_alpha_wrap = false)
{
  Mesh wrap;
  if(use_noprmal_alpha_wrap)
  {
    CGAL::alpha_wrap_3(input_mesh, alpha, offset, wrap);
  }
  else
  {
    CGAL::alpha_wrap_3(input_mesh, alpha, offset, wrap,
                       CGAL::parameters::max_distance_to_input_in_offsets(tau));
  }
  return wrap;
}

inline Mesh wrap_point_cloud(const Point_container& points,
                             const double alpha,
                             const double offset,
                             const double tau,
                             const bool use_noprmal_alpha_wrap = false)
{
  Mesh wrap;
  if(use_noprmal_alpha_wrap)
  {
    CGAL::alpha_wrap_3(points, alpha, offset, wrap);
  }
  else
  {
    CGAL::alpha_wrap_3(points, alpha, offset, wrap,
                       CGAL::parameters::max_distance_to_input_in_offsets(tau));
  }
  return wrap;
}

inline Wrap_result run_alpha_wrap(const Wrap_request& request)
{
  if(request.input_path.empty())
    throw std::runtime_error("Input path is empty.");
  if(!is_positive_finite(request.relative_alpha))
    throw std::runtime_error("alpha must be a positive finite value.");
  if(!is_positive_finite(request.relative_offset))
    throw std::runtime_error("offset must be a positive finite value.");
  if(!(std::isfinite(request.tau) && request.tau > 1.0))
    throw std::runtime_error("tau must be finite and strictly larger than 1.0.");

  const Input_kind kind = detect_input_kind(request.input_path);
  if(kind == Input_kind::Unknown)
    throw std::runtime_error("Input is neither a valid triangle mesh nor a valid point cloud: " + request.input_path);

  if(kind == Input_kind::Triangle_mesh)
  {
    bool input_was_triangulated = false;
    Mesh input_mesh = read_mesh_or_throw(request.input_path, &input_was_triangulated);
    if(input_was_triangulated)
      std::cout << "Input is not a triangle mesh. It will first be triangulated.\n";
    const auto [alpha, offset] = compute_absolute_alpha_offset(input_mesh, request.relative_alpha, request.relative_offset);
    std::cout << "Input type: triangle mesh (" << CGAL::num_faces(input_mesh) << " faces)\n";
    std::cout << "Running beeren_method with alpha=" << alpha << ", offset=" << offset << ", tau=" << request.tau << "\n";
    return Wrap_result{kind,
                       wrap_triangle_mesh(input_mesh, alpha, offset, request.tau, request.use_noprmal_alpha_wrap),
                       alpha, offset, CGAL::num_faces(input_mesh)};
  }

  Point_container points = read_points_or_throw(request.input_path);
  const auto [alpha, offset] = compute_absolute_alpha_offset(points, request.relative_alpha, request.relative_offset);
  std::cout << "Input type: point cloud (" << points.size() << " points)\n";
  std::cout << "Running beeren_method with alpha=" << alpha << ", offset=" << offset << ", tau=" << request.tau << "\n";
  return Wrap_result{kind,
                     wrap_point_cloud(points, alpha, offset, request.tau, request.use_noprmal_alpha_wrap),
                     alpha, offset, points.size()};
}
  
inline bool write_output_mesh(const std::string& output_path, const Mesh& mesh)
{
  std::filesystem::path p(output_path);
  if(p.has_parent_path())
    std::filesystem::create_directories(p.parent_path());

  return CGAL::IO::write_polygon_mesh(output_path, mesh, CGAL::parameters::stream_precision(17));
}

} // namespace cli_helpers

#endif // THESIS_CLI_ALPHA_WRAPPING_H
