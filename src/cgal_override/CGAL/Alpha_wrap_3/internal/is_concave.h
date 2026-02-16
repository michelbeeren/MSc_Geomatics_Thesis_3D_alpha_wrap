#ifndef CGAL_ALPHA_WRAP_3_IS_CONCAVE_H
#define CGAL_ALPHA_WRAP_3_IS_CONCAVE_H

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h> // Required for searching
#include <vector>

namespace CGAL {
    namespace Alpha_wrap_3 {
        namespace internal {

            template <class Kernel>
            class Octree_refinement_depth {
                using Point_3 = typename Kernel::Point_3;
                using Bbox_3  = CGAL::Bbox_3;
                using Traits  = CGAL::Search_traits_3<Kernel>;
                using Tree    = CGAL::Kd_tree<Traits>;
                using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits>;

                Tree m_tree;
                bool m_has_data;
                double m_finest_edge_length;

            public:
                Octree_refinement_depth(const std::vector<Bbox_3>& finest_cubes)
                    : m_has_data(false), m_finest_edge_length(0.0) {
                    if (!finest_cubes.empty()) {
                        // Calculate edge length once from the first cube
                        m_finest_edge_length = finest_cubes[0].xmax() - finest_cubes[0].xmin();

                        for(const auto& b : finest_cubes) {
                            m_tree.insert(Point_3((b.xmin()+b.xmax())/2.0,
                                                  (b.ymin()+b.ymax())/2.0,
                                                  (b.zmin()+b.zmax())/2.0));
                        }
                        m_tree.build();
                        m_has_data = true;
                    }
                }

                double finest_edge_length() const { return m_finest_edge_length; }

                double squared_dist_to_finest(const Point_3& p) const {
                    if (!m_has_data) return (std::numeric_limits<double>::max)();

                    // Standard Kd-tree search for the single nearest neighbor (k=1)
                    Neighbor_search search(m_tree, p, 1);

                    // The search returns a pair: [Point, Squared_Distance]
                    if (search.begin() != search.end()) {
                        return CGAL::to_double(search.begin()->second);
                    }

                    return (std::numeric_limits<double>::max)();
                }
            };

        } // namespace internal
    } // namespace Alpha_wrap_3
} // namespace CGAL
#endif