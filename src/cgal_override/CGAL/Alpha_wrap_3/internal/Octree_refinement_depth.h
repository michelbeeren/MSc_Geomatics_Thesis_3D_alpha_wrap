#ifndef CGAL_ALPHA_WRAP_3_OCTREE_REFINEMENT_DEPTH_H
#define CGAL_ALPHA_WRAP_3_OCTREE_REFINEMENT_DEPTH_H

#include "../../../../octree.h"// Ensure this path is correct for your project structure
//

namespace CGAL {
namespace Alpha_wrap_3 {
namespace internal {

template <class Kernel>
class Octree_refinement_depth {
    std::unique_ptr<OctreeCell> m_root;

public:
    using Point_3 = typename Kernel::Point_3;
    // Constructor takes ownership of the built octree root
    Octree_refinement_depth(std::unique_ptr<OctreeCell> root)
        : m_root(std::move(root)) {}

    /**
     * Finds the leaf cell containing the point p and returns its depth.
     * If the point is outside the octree, it returns 0 (the root depth).
     */
    double nearest_refinement_depth(const Point_3& p) const {
        if (!m_root) return 0.0;

        // If the point is completely outside the octree, return base depth
        if (p.x() < m_root->bbox.xmin() || p.x() > m_root->bbox.xmax() ||
            p.y() < m_root->bbox.ymin() || p.y() > m_root->bbox.ymax() ||
            p.z() < m_root->bbox.zmin() || p.z() > m_root->bbox.zmax()) {
            return 0.0;
        }

        return find_depth_recursive(m_root.get(), p);
    }

private:
    double find_depth_recursive(const OctreeCell* cell, const Point_3& p) const {
        // Base case: if it's a leaf, return current depth
        bool is_leaf = true;
        for (const auto& child : cell->children) {
            if (child) {
                is_leaf = false;
                break;
            }
        }

        if (is_leaf) {
            return static_cast<double>(cell->depth);
        }

        // Recursive step: Find which child contains the point
        // Using the same midpoint logic as your subdivide() function
        const CGAL::Bbox_3& b = cell->bbox;
        double mx = 0.5 * (b.xmin() + b.xmax());
        double my = 0.5 * (b.ymin() + b.ymax());
        double mz = 0.5 * (b.zmin() + b.zmax());

        // Standard Octree indexing (as defined in your subdivide function):
        // 0: min_x, min_y, min_z | 1: max_x, min_y, min_z | 2: max_x, max_y, min_z ...
        int index = 0;
        if (p.x() >= mx) index += 1;
        if (p.y() >= my) index += 2; // Assuming index 2 & 3 have y > my based on your code
        if (p.z() >= mz) index += 4;

        // Note: Check your subdivide() layout to ensure the index matches perfectly.
        // Based on your specific subdivide implementation:
        // children[0-3] are z < mz, children[4-7] are z >= mz.

        // Safety check for null children (though a balanced tree shouldn't have them)
        if (cell->children[index]) {
            return find_depth_recursive(cell->children[index].get(), p);
        }

        return static_cast<double>(cell->depth);
    }
};

} // namespace internal
} // namespace Alpha_wrap_3
} // namespace CGAL
#endif // CGAL_ALPHA_WRAP_3_OCTREE_REFINEMENT_DEPTH_H
