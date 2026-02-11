#ifndef CGAL_ALPHA_WRAP_3_OCTREE_REFINEMENT_DEPTH_H
#define CGAL_ALPHA_WRAP_3_OCTREE_REFINEMENT_DEPTH_H

#include "../../../../octree.h"
#include <cmath>

namespace CGAL {
    namespace Alpha_wrap_3 {
        namespace internal {

            template <class Kernel>
            class Octree_refinement_depth {
                std::vector<::Thesis::OctreeCell> m_cells;

            public:
                using Point_3 = typename Kernel::Point_3;

                Octree_refinement_depth(std::vector<::Thesis::OctreeCell> cells)
                    : m_cells(std::move(cells)) {}

                double nearest_refinement_depth(const Point_3& p) const {
                    if (m_cells.empty()) return 0.0;

                    int deepest_found = 0;
                    bool found = false;

                    for (const auto& cell : m_cells) {
                        // Check if point is inside this bbox
                        if (p.x() >= cell.bbox.xmin() && p.x() <= cell.bbox.xmax() &&
                            p.y() >= cell.bbox.ymin() && p.y() <= cell.bbox.ymax() &&
                            p.z() >= cell.bbox.zmin() && p.z() <= cell.bbox.zmax()) {

                            // If we find multiple, take the deepest one (highest resolution)
                            if (cell.depth > deepest_found || !found) {
                                deepest_found = cell.depth;
                                found = true;
                            }
                            }
                    }

                    return static_cast<double>(deepest_found);
                }
            };

        } // namespace internal
    } // namespace Alpha_wrap_3
} // namespace CGAL
#endif