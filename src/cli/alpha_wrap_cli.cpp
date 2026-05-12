#include "helper_files/alpha_wrapping.h"
#include "helper_files/hausdorff.h"
#include "helper_files/val3dity.h"

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;

struct Arguments
{
  std::string input_path;
  std::string output_path;
  double alpha = 0.0;
  double offset = 0.0;
  double tau = 0.0;
  bool validate = false;
  bool wrap_invalid_only = false;
  bool use_noprmal_alpha_wrap = false;
};

void print_usage(const char* executable)
{
  std::cerr << "Usage:\n"
            << "  " << executable << " <input_path> <output_path> <alpha> <offset> <tau> "
            << "[--validate] [--wrap_invalid_only] [--use_noprmal_alpha_wrap]\n\n"
            << "Notes:\n"
            << "  - alpha and offset are absolute values\n"
            << "  - tau is max_distance_to_input_in_offsets and must be > 1.0\n"
            << "  - add --validate to run mesh validation\n"
            << "  - add --wrap_invalid_only to validate each connected component, wrap invalid ones,\n"
            << "    and wrap intersecting disconnected component groups together\n"
            << "  - add --use_noprmal_alpha_wrap to use CGAL alpha_wrap_3 without tau parameter\n"
            << "  - --wrap_invalid_only is only supported for triangle-mesh inputs (not point clouds)\n";
}

bool parse_arguments(const int argc, char** argv, Arguments& args)
{
  if(argc < 6)
    return false;

  args.input_path = argv[1];
  args.output_path = argv[2];

  try
  {
    args.alpha = std::stod(argv[3]);
    args.offset = std::stod(argv[4]);
    args.tau = std::stod(argv[5]);
  }
  catch(const std::exception&)
  {
    return false;
  }

  for(int i = 6; i < argc; ++i)
  {
    const std::string flag(argv[i]);
    if(flag == "--validate")
      args.validate = true;
    else if(flag == "--wrap_invalid_only")
      args.wrap_invalid_only = true;
    else if(flag == "--use_noprmal_alpha_wrap" || flag == "--use_normal_alpha_wrap")
      args.use_noprmal_alpha_wrap = true;
    else
      return false;
  }

  return true;
}

struct Selective_wrap_result
{
  cli_helpers::Mesh wrap_mesh;
  std::size_t component_count = 0;
  std::size_t invalid_component_count = 0;
  std::size_t wrapped_group_count = 0;
  std::size_t kept_component_count = 0;
};

class Disjoint_set
{
public:
  explicit Disjoint_set(const std::size_t n)
      : parent_(n), rank_(n, 0)
  {
    std::iota(parent_.begin(), parent_.end(), 0);
  }

  std::size_t find(const std::size_t x)
  {
    if(parent_[x] != x)
      parent_[x] = find(parent_[x]);
    return parent_[x];
  }

  void unite(const std::size_t a, const std::size_t b)
  {
    std::size_t ra = find(a);
    std::size_t rb = find(b);
    if(ra == rb)
      return;
    if(rank_[ra] < rank_[rb])
      std::swap(ra, rb);
    parent_[rb] = ra;
    if(rank_[ra] == rank_[rb])
      ++rank_[ra];
  }

private:
  std::vector<std::size_t> parent_;
  std::vector<unsigned int> rank_;
};

double bbox_distance_lower_bound(const CGAL::Bbox_3& a, const CGAL::Bbox_3& b)
{
  const double dx = (a.xmax() < b.xmin()) ? (b.xmin() - a.xmax()) :
                    (b.xmax() < a.xmin()) ? (a.xmin() - b.xmax()) : 0.0;
  const double dy = (a.ymax() < b.ymin()) ? (b.ymin() - a.ymax()) :
                    (b.ymax() < a.ymin()) ? (a.ymin() - b.ymax()) : 0.0;
  const double dz = (a.zmax() < b.zmin()) ? (b.zmin() - a.zmax()) :
                    (b.zmax() < a.zmin()) ? (a.zmin() - b.zmax()) : 0.0;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

class Component_distance_cache
{
public:
  using Primitive = CGAL::AABB_face_graph_triangle_primitive<cli_helpers::Mesh>;
  using Traits = CGAL::AABB_traits<cli_helpers::K, Primitive>;
  using Tree = CGAL::AABB_tree<Traits>;

  explicit Component_distance_cache(const std::vector<cli_helpers::Mesh>& components)
      : components_(components),
        distance_cache_(components.size(), std::vector<double>(components.size(), -1.0))
  {
    bboxes_.reserve(components_.size());
    trees_.reserve(components_.size());
    for(const auto& mesh : components_)
    {
      bboxes_.push_back(PMP::bbox(mesh));
      auto tree = std::make_unique<Tree>(faces(mesh).begin(), faces(mesh).end(), mesh);
      tree->accelerate_distance_queries();
      trees_.push_back(std::move(tree));
    }
  }

  double lower_bound(const std::size_t i, const std::size_t j) const
  {
    return bbox_distance_lower_bound(bboxes_[i], bboxes_[j]);
  }

  double distance(const std::size_t i, const std::size_t j)
  {
    if(i == j)
      return 0.0;
    const std::size_t a = (std::min)(i, j);
    const std::size_t b = (std::max)(i, j);
    double& cached = distance_cache_[a][b];
    if(cached >= 0.0)
      return cached;
    cached = compute_distance(a, b);
    return cached;
  }

private:
  double min_vertex_to_mesh_distance(const std::size_t source, const std::size_t target) const
  {
    const cli_helpers::Mesh& src = components_[source];
    const Tree& trg_tree = *trees_[target];
    double min_sq = std::numeric_limits<double>::infinity();
    for(const auto v : src.vertices())
    {
      const double sq = trg_tree.squared_distance(src.point(v));
      if(sq < min_sq)
        min_sq = sq;
    }
    if(!std::isfinite(min_sq))
      return std::numeric_limits<double>::infinity();
    return std::sqrt(min_sq);
  }

  double compute_distance(const std::size_t i, const std::size_t j)
  {
    if(PMP::do_intersect(components_[i], components_[j]))
      return 0.0;
    const double d_ij = min_vertex_to_mesh_distance(i, j);
    const double d_ji = min_vertex_to_mesh_distance(j, i);
    return (std::min)(d_ij, d_ji);
  }

  const std::vector<cli_helpers::Mesh>& components_;
  std::vector<CGAL::Bbox_3> bboxes_;
  std::vector<std::unique_ptr<Tree>> trees_;
  std::vector<std::vector<double>> distance_cache_;
};

struct Base_group
{
  std::vector<std::size_t> members;
  bool is_special = false; // intersecting group and/or invalid component
};

Selective_wrap_result run_wrap_invalid_components_only(const Arguments& args)
{
  const cli_helpers::Input_kind kind = cli_helpers::detect_input_kind(args.input_path);
  if(kind == cli_helpers::Input_kind::Point_cloud)
    throw std::runtime_error("--wrap_invalid_only can only be used with triangle-mesh inputs, not point clouds.");
  if(kind != cli_helpers::Input_kind::Triangle_mesh)
    throw std::runtime_error("Input is neither a valid triangle mesh nor a valid point cloud: " + args.input_path);

  bool input_was_triangulated = false;
  const cli_helpers::Mesh input_mesh = cli_helpers::read_mesh_or_throw(args.input_path, &input_was_triangulated);
  if(input_was_triangulated)
    std::cout << "Input is not a triangle mesh. It will first be triangulated.\n";
  std::vector<cli_helpers::Mesh> components;
  PMP::split_connected_components(input_mesh, components);

  const double alpha = args.alpha;
  const double offset = args.offset;

  std::cout << "Finding invalid groups in input\n";

  // 1) Per-component validity
  std::vector<bool> component_validity(components.size(), false);
  std::size_t invalid_component_count = 0;
  for(std::size_t i = 0; i < components.size(); ++i)
  {
    const cli_helpers::Mesh& component = components[i];
    const bool is_valid = cli_helpers::valid_mesh_boolean(component, false);
    component_validity[i] = is_valid;
    if(!is_valid)
      ++invalid_component_count;
  }

  // 2) Base grouping by geometric intersections
  Disjoint_set intersection_groups(components.size());
  std::size_t intersecting_pairs = 0;
  std::vector<CGAL::Bbox_3> component_bboxes;
  component_bboxes.reserve(components.size());
  for(const auto& component : components)
    component_bboxes.push_back(PMP::bbox(component));

  for(std::size_t i = 0; i < components.size(); ++i)
  {
    for(std::size_t j = i + 1; j < components.size(); ++j)
    {
      if(!CGAL::do_overlap(component_bboxes[i], component_bboxes[j]))
        continue;
      if(PMP::do_intersect(components[i], components[j]))
      {
        intersection_groups.unite(i, j);
        ++intersecting_pairs;
      }
    }
  }

  std::unordered_map<std::size_t, std::vector<std::size_t>> grouped_by_root;
  grouped_by_root.reserve(components.size());
  for(std::size_t i = 0; i < components.size(); ++i)
    grouped_by_root[intersection_groups.find(i)].push_back(i);

  std::vector<Base_group> base_groups;
  base_groups.reserve(grouped_by_root.size());
  std::unordered_map<std::size_t, std::size_t> root_to_group_idx;
  root_to_group_idx.reserve(grouped_by_root.size());

  for(const auto& [root, members] : grouped_by_root)
  {
    const std::size_t idx = base_groups.size();
    root_to_group_idx[root] = idx;
    Base_group group;
    group.members = members;
    bool has_invalid = false;
    for(const std::size_t c : members)
    {
      if(!component_validity[c])
      {
        has_invalid = true;
        break;
      }
    }
    group.is_special = (members.size() > 1) || has_invalid;
    base_groups.push_back(std::move(group));
  }

  std::vector<std::size_t> component_to_base_group(components.size(), 0);
  for(std::size_t c = 0; c < components.size(); ++c)
    component_to_base_group[c] = root_to_group_idx[intersection_groups.find(c)];

  (void)intersecting_pairs;

  // 3) Proximity-based merging starting from each special group
  Component_distance_cache distance_cache(components);
  Disjoint_set proximity_groups(base_groups.size());
  std::size_t proximity_merges = 0;

  for(std::size_t g = 0; g < base_groups.size(); ++g)
  {
    if(!base_groups[g].is_special)
      continue;

    double best_distance = std::numeric_limits<double>::infinity();
    std::size_t best_component = components.size();

    for(const std::size_t source_component : base_groups[g].members)
    {
      for(std::size_t other_component = 0; other_component < components.size(); ++other_component)
      {
        if(component_to_base_group[other_component] == g)
          continue;
        const double lower = distance_cache.lower_bound(source_component, other_component);
        if(lower >= best_distance)
          continue;
        const double d = distance_cache.distance(source_component, other_component);
        if(d < best_distance)
        {
          best_distance = d;
          best_component = other_component;
        }
      }
    }

    if(best_component == components.size())
      continue;

    const std::size_t target_group = component_to_base_group[best_component];
    const bool target_is_special = base_groups[target_group].is_special;
    const double threshold = target_is_special ? 2.5 * offset : 1.5 * offset;

    if(best_distance < threshold)
    {
      const std::size_t before_a = proximity_groups.find(g);
      const std::size_t before_b = proximity_groups.find(target_group);
      if(before_a != before_b)
      {
        proximity_groups.unite(g, target_group);
        ++proximity_merges;
      }
    }
  }

  (void)proximity_merges;

  // 4) Build final groups after proximity merging
  std::unordered_map<std::size_t, std::vector<std::size_t>> final_group_to_base_groups;
  final_group_to_base_groups.reserve(base_groups.size());
  for(std::size_t g = 0; g < base_groups.size(); ++g)
    final_group_to_base_groups[proximity_groups.find(g)].push_back(g);

  std::vector<std::vector<std::size_t>> groups_to_wrap;
  std::vector<std::vector<std::size_t>> groups_to_keep;
  groups_to_wrap.reserve(final_group_to_base_groups.size());
  groups_to_keep.reserve(final_group_to_base_groups.size());
  std::size_t kept_component_count = 0;

  for(const auto& [root, grouped_base_ids] : final_group_to_base_groups)
  {
    (void)root;
    bool should_wrap_group = false;
    for(const std::size_t bg : grouped_base_ids)
    {
      if(base_groups[bg].is_special)
      {
        should_wrap_group = true;
        break;
      }
    }
    if(should_wrap_group)
      groups_to_wrap.push_back(grouped_base_ids);
    else
    {
      groups_to_keep.push_back(grouped_base_ids);
      for(const std::size_t bg : grouped_base_ids)
        kept_component_count += base_groups[bg].members.size();
    }
  }

  std::cout << "Groups selected for wrapping (invalid/intersecting/too_close): "
            << groups_to_wrap.size()
            << " | Connected components kept unchanged: "
            << kept_component_count << "\n";
  std::cout << "Starting wrap of selected groups with alpha=" << alpha
            << ", offset=" << offset << ", tau=" << args.tau << "\n";

  cli_helpers::Mesh output_mesh;
  for(const auto& grouped_base_ids : groups_to_keep)
  {
    for(const std::size_t bg : grouped_base_ids)
    {
      for(const std::size_t c : base_groups[bg].members)
        CGAL::copy_face_graph(components[c], output_mesh);
    }
  }

  for(const auto& grouped_base_ids : groups_to_wrap)
  {
    cli_helpers::Mesh group_input;
    for(const std::size_t bg : grouped_base_ids)
    {
      for(const std::size_t c : base_groups[bg].members)
        CGAL::copy_face_graph(components[c], group_input);
    }

    const cli_helpers::Mesh wrapped_component =
        cli_helpers::wrap_triangle_mesh(group_input, alpha, offset, args.tau, args.use_noprmal_alpha_wrap);
    CGAL::copy_face_graph(wrapped_component, output_mesh);
  }

  return {output_mesh, components.size(), invalid_component_count, groups_to_wrap.size(), kept_component_count};
}

void rewrap_intersections_until_clean(cli_helpers::Mesh& mesh, const Arguments& args)
{
  constexpr std::size_t max_intersection_iterations = 50;

  for(std::size_t iteration = 1; iteration <= max_intersection_iterations; ++iteration)
  {
    std::vector<cli_helpers::Mesh> components;
    PMP::split_connected_components(mesh, components);
    if(components.size() < 2)
      return;

    std::vector<bool> intersects(components.size(), false);
    std::vector<CGAL::Bbox_3> bboxes;
    bboxes.reserve(components.size());
    for(const auto& component : components)
      bboxes.push_back(PMP::bbox(component));

    for(std::size_t i = 0; i < components.size(); ++i)
    {
      for(std::size_t j = i + 1; j < components.size(); ++j)
      {
        if(!CGAL::do_overlap(bboxes[i], bboxes[j]))
          continue;
        if(PMP::do_intersect(components[i], components[j]))
        {
          intersects[i] = true;
          intersects[j] = true;
        }
      }
    }

    std::size_t intersecting_component_count = 0;
    for(const bool v : intersects)
      if(v) ++intersecting_component_count;

    if(intersecting_component_count == 0)
      return;

    std::cout << "Post-wrap intersection check: found " << intersecting_component_count
              << " intersecting connected components; rewrapping them together\n";

    cli_helpers::Mesh keep_mesh;
    cli_helpers::Mesh intersect_mesh;
    for(std::size_t i = 0; i < components.size(); ++i)
    {
      if(intersects[i])
        CGAL::copy_face_graph(components[i], intersect_mesh);
      else
        CGAL::copy_face_graph(components[i], keep_mesh);
    }

    const cli_helpers::Mesh rewrapped_intersections =
        cli_helpers::wrap_triangle_mesh(intersect_mesh, args.alpha, args.offset, args.tau, args.use_noprmal_alpha_wrap);

    cli_helpers::Mesh repaired_mesh;
    CGAL::copy_face_graph(keep_mesh, repaired_mesh);
    CGAL::copy_face_graph(rewrapped_intersections, repaired_mesh);
    mesh = std::move(repaired_mesh);
  }

  throw std::runtime_error("Post-wrap intersection cleanup did not converge within 50 iterations.");
}

void validate_and_repair_until_valid(cli_helpers::Mesh& mesh, const Arguments& args)
{
  constexpr std::size_t max_validation_iterations = 50;

  for(std::size_t iteration = 1; iteration <= max_validation_iterations; ++iteration)
  {
    std::vector<cli_helpers::Mesh> validation_components;
    PMP::split_connected_components(mesh, validation_components);
    const std::size_t total_components = validation_components.size();

    std::cout << "Validation: identified " << total_components
              << " connected components --> validating each\n";

    cli_helpers::Mesh valid_components_mesh;
    cli_helpers::Mesh invalid_components_mesh;
    std::size_t valid_components = 0;
    std::size_t invalid_components = 0;

    for(std::size_t i = 0; i < total_components; ++i)
    {
      const bool is_valid = cli_helpers::valid_mesh_boolean(validation_components[i], false);
      if(is_valid)
      {
        ++valid_components;
        CGAL::copy_face_graph(validation_components[i], valid_components_mesh);
      }
      else
      {
        ++invalid_components;
        CGAL::copy_face_graph(validation_components[i], invalid_components_mesh);
      }
    }

    std::cout << "Validation result: [" << valid_components << "/" << total_components << "] VALID.\n";


    std::cout << "Rewrapping " << invalid_components << " invalid connected components together\n";
    const cli_helpers::Mesh rewrapped_invalid =
        cli_helpers::wrap_triangle_mesh(invalid_components_mesh, args.alpha, args.offset, args.tau, args.use_noprmal_alpha_wrap);

    cli_helpers::Mesh repaired_mesh;
    CGAL::copy_face_graph(valid_components_mesh, repaired_mesh);
    CGAL::copy_face_graph(rewrapped_invalid, repaired_mesh);
    mesh = std::move(repaired_mesh);
  }

  throw std::runtime_error("Validation-repair loop did not converge within 50 iterations.");
}

int main(int argc, char** argv)
{
  Arguments args;
  if(!parse_arguments(argc, argv, args))
  {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  try
  {
    if(!cli_helpers::is_positive_finite(args.alpha))
      throw std::runtime_error("alpha must be a positive finite value.");
    if(!cli_helpers::is_positive_finite(args.offset))
      throw std::runtime_error("offset must be a positive finite value.");
    if(!(std::isfinite(args.tau) && args.tau > 1.0))
      throw std::runtime_error("tau must be finite and strictly larger than 1.0.");

    const auto wrap_start = std::chrono::steady_clock::now();

    cli_helpers::Mesh output_mesh;
    if(args.wrap_invalid_only)
    {
      Selective_wrap_result result = run_wrap_invalid_components_only(args);
      output_mesh = std::move(result.wrap_mesh);
      std::cout << "Checking wrapped output for remaining intersections\n";
      rewrap_intersections_until_clean(output_mesh, args);
    }
    else
    {
      const cli_helpers::Input_kind kind = cli_helpers::detect_input_kind(args.input_path);
      if(kind == cli_helpers::Input_kind::Unknown)
        throw std::runtime_error("Input is neither a valid triangle mesh nor a valid point cloud: " + args.input_path);

      if(kind == cli_helpers::Input_kind::Triangle_mesh)
      {
        bool input_was_triangulated = false;
        cli_helpers::Mesh input_mesh = cli_helpers::read_mesh_or_throw(args.input_path, &input_was_triangulated);
        if(input_was_triangulated)
          std::cout << "Input is not a triangle mesh. It will first be triangulated.\n";
        std::cout << "Input type: triangle mesh (" << CGAL::num_faces(input_mesh) << " faces)\n";
        std::cout << "Starting alpha wrapping with parameters alpha=" << args.alpha
                  << ", offset=" << args.offset << ", tau=" << args.tau << "\n";
        output_mesh = cli_helpers::wrap_triangle_mesh(input_mesh, args.alpha, args.offset, args.tau, args.use_noprmal_alpha_wrap);
      }
      else
      {
        cli_helpers::Point_container points = cli_helpers::read_points_or_throw(args.input_path);
        std::cout << "Input type: point cloud (" << points.size() << " points)\n";
        std::cout << "Starting alpha wrapping with parameters alpha=" << args.alpha
                  << ", offset=" << args.offset << ", tau=" << args.tau << "\n";
        output_mesh = cli_helpers::wrap_point_cloud(points, args.alpha, args.offset, args.tau, args.use_noprmal_alpha_wrap);
      }
    }

    if(args.validate)
    {
      std::cout << "Wrapping finished. Output will now be validated.\n";
      validate_and_repair_until_valid(output_mesh, args);
    }

    const auto wrap_end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> wrap_elapsed = wrap_end - wrap_start;

    std::cout << "Wrap result: " << num_vertices(output_mesh) << " vertices, "
              << num_faces(output_mesh) << " faces\n";
    std::cout << "Total wrapping time: " << std::fixed << std::setprecision(3)
              << wrap_elapsed.count() << " s\n";
    std::cout << "Writing output mesh: " << args.output_path << "\n";
    if(!cli_helpers::write_output_mesh(args.output_path, output_mesh))
    {
      std::cerr << "Error: could not write output mesh.\n";
      return EXIT_FAILURE;
    }

    std::cout << "Done.\n";
    return EXIT_SUCCESS;
  }
  catch(const std::exception& e)
  {
    std::cerr << "Error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
}
