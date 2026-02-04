//
// Created by Michel Beeren on 04/02/2026.
//

#ifndef CGAL_ALPHA_WRAP_3_FEATURE_SIZE_MAT_H
#define CGAL_ALPHA_WRAP_3_FEATURE_SIZE_MAT_H

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <memory>

#include <boost/property_map/property_map.hpp>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

namespace CGAL {
namespace Alpha_wrap_3 {
namespace internal {

// Header-only: loads ASCII PLY (x y z feature_size) and provides nearest feature size.
template <class Kernel>
class Feature_size_MAT
{
public:
  using Point_3 = typename Kernel::Point_3;

  Feature_size_MAT() = default;

  explicit Feature_size_MAT(std::string ply_path)
    : m_ply_path(std::move(ply_path))
  {}

  void set_ply_path(std::string ply_path)
  {
    m_ply_path = std::move(ply_path);
    clear();
  }

  // Return feature_size of nearest PLY point to query.
  double nearest_feature_size(const Point_3& q) const
  {
    ensure_loaded();

    Neighbor_search search(*m_tree, q, 1);
    auto it = search.begin();
    if (it == search.end())
      return 1.0; // fallback

    const std::size_t id = it->first.id;
    return (id < m_feature_sizes.size()) ? m_feature_sizes[id] : 1.0;
  }

private:
  struct FS_Point
  {
    Point_3 p;
    std::size_t id;
  };

  // property map: FS_Point -> Point_3
  struct FS_Point_map
  {
    typedef FS_Point key_type;
    typedef Point_3 value_type;
    typedef const value_type& reference;
    typedef boost::readable_property_map_tag category;
  };

  friend inline typename FS_Point_map::reference get(FS_Point_map, const FS_Point& v)
  {
    return v.p;
  }

  using Base_traits = CGAL::Search_traits_3<Kernel>;
  using Traits      = CGAL::Search_traits_adapter<FS_Point, FS_Point_map, Base_traits>;
  using Tree        = CGAL::Kd_tree<Traits>;
  using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits>;

private:
  void clear() const
  {
    m_loaded = false;
    m_pts.clear();
    m_feature_sizes.clear();
    m_tree.reset();
  }

  void ensure_loaded() const
  {
    if (m_loaded) return;
    if (m_ply_path.empty())
      throw std::runtime_error("Feature_size_MAT: PLY path is empty.");

    std::ifstream in(m_ply_path.c_str());
    if (!in)
      throw std::runtime_error("Feature_size_MAT: cannot open PLY: " + m_ply_path);

    // Skip header until end_header
    std::string line;
    while (std::getline(in, line))
    {
      if (line == "end_header") break;
    }

    // Read vertices: x y z feature_size
    double x, y, z, fs;
    std::size_t id = 0;

    while (in >> x >> y >> z >> fs)
    {
      FS_Point v;
      v.p  = Point_3(x, y, z);
      v.id = id;

      m_pts.push_back(v);
      m_feature_sizes.push_back(fs);
      ++id;
    }

    if (m_pts.empty())
      throw std::runtime_error("Feature_size_MAT: no vertices parsed from PLY: " + m_ply_path);

    m_tree = std::make_unique<Tree>(m_pts.begin(), m_pts.end());
    m_tree->build();

    m_loaded = true;
  }

private:
  std::string m_ply_path;

  mutable bool m_loaded = false;
  mutable std::vector<FS_Point> m_pts;
  mutable std::vector<double>   m_feature_sizes;
  mutable std::unique_ptr<Tree> m_tree;
};

} // namespace internal
} // namespace Alpha_wrap_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_FEATURE_SIZE_MAT_H
