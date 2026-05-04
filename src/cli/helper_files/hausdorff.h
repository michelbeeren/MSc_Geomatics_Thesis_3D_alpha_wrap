#ifndef THESIS_CLI_HAUSDORFF_H
#define THESIS_CLI_HAUSDORFF_H

#include "alpha_wrapping.h"

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/IO/Color.h>

#include <fstream>
#include <algorithm>
#include <limits>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace cli_helpers
{

inline std::vector<double> point_to_mesh_distances(const std::vector<Point_3>& sample_points,
                                                   const Mesh& target_mesh)
{
  using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
  using AABB_traits = CGAL::AABB_traits<K, Primitive>;
  using Tree = CGAL::AABB_tree<AABB_traits>;

  Tree tree(faces(target_mesh).begin(), faces(target_mesh).end(), target_mesh);
  tree.accelerate_distance_queries();

  std::vector<double> distances;
  distances.reserve(sample_points.size());
  for(const Point_3& p : sample_points)
    distances.push_back(std::sqrt(tree.squared_distance(p)));

  return distances;
}

inline void write_distances_to_csv(const std::string& filename, const std::vector<double>& distances)
{
  std::ofstream out(filename);
  if(!out)
    throw std::runtime_error("Could not open file: " + filename);

  out << "index,distance\n";
  for(std::size_t i = 0; i < distances.size(); ++i)
    out << i << "," << distances[i] << "\n";
}

inline void hausdorff_distance(Mesh& wrapped_mesh, const Mesh& original_mesh)
{
  using fd = Mesh::Face_index;
  using vd = Mesh::Vertex_index;
  using Color = CGAL::IO::Color;
  using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
  using AABB_traits = CGAL::AABB_traits<K, Primitive>;
  using Tree = CGAL::AABB_tree<AABB_traits>;

  Tree tree(faces(original_mesh).begin(), faces(original_mesh).end(), original_mesh);
  tree.accelerate_distance_queries();

  Mesh::Property_map<fd, Color> fcolor;
  bool created_color = false;
  boost::tie(fcolor, created_color) = wrapped_mesh.add_property_map<fd, Color>("f:color", Color(255, 255, 255));
  if(!created_color)
    return;

  double dmin = std::numeric_limits<double>::max();
  double dmax = 0.0;
  std::vector<double> face_distances(num_faces(wrapped_mesh), 0.0);

  for(fd f : wrapped_mesh.faces())
  {
    const auto h = wrapped_mesh.halfedge(f);
    const vd v0 = wrapped_mesh.source(h);
    const vd v1 = wrapped_mesh.target(h);
    const vd v2 = wrapped_mesh.target(wrapped_mesh.next(h));

    const Point_3& p0 = wrapped_mesh.point(v0);
    const Point_3& p1 = wrapped_mesh.point(v1);
    const Point_3& p2 = wrapped_mesh.point(v2);

    const Point_3 c((p0.x() + p1.x() + p2.x()) / 3.0,
                    (p0.y() + p1.y() + p2.y()) / 3.0,
                    (p0.z() + p1.z() + p2.z()) / 3.0);

    const double dist = std::sqrt(tree.squared_distance(c));
    face_distances[f] = dist;
    dmin = std::min(dmin, dist);
    dmax = std::max(dmax, dist);
  }

  const double range = (dmax > dmin) ? (dmax - dmin) : 1.0;
  for(fd f : wrapped_mesh.faces())
  {
    const double t = (face_distances[f] - dmin) / range;
    const unsigned char r = static_cast<unsigned char>(255.0 * t);
    const unsigned char b = static_cast<unsigned char>(255.0 * (1.0 - t));
    fcolor[f] = Color(r, 0, b);
  }
}

inline void write_hausdorff_distance(const Mesh& wrapped_mesh, const std::string& filename)
{
  std::string out_name = filename;
  if(out_name.size() > 4 && out_name.substr(out_name.size() - 4) == ".off")
    out_name = out_name.substr(0, out_name.size() - 4);
  out_name += "_Hausdorff.ply";

  using FaceColorMap = Mesh::Property_map<Mesh::Face_index, CGAL::IO::Color>;
  const std::optional<FaceColorMap> face_colors_opt =
      wrapped_mesh.property_map<Mesh::Face_index, CGAL::IO::Color>("f:color");
  if(!face_colors_opt)
    return;

  std::ofstream out(out_name);
  if(!out)
    return;

  CGAL::IO::write_PLY(out, wrapped_mesh, CGAL::parameters::stream_precision(17));
}

} // namespace cli_helpers

#endif // THESIS_CLI_HAUSDORFF_H
