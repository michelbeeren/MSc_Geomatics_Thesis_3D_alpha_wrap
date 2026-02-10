//
// Created by Michel Beeren on 10/02/2026.
//

#ifndef CGAL_ALPHA_WRAP_3_OCTREE_REFINEMENT_DEPTH_H
#define CGAL_ALPHA_WRAP_3_OCTREE_REFINEMENT_DEPTH_H

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <memory>

#include <CGAL/Bbox_3.h>  // For Bbox_3, assuming you're using CGAL bounding boxes
#include <CGAL/Point_3.h>  // For Point_3, assuming you're using CGAL points

namespace CGAL {
namespace Alpha_wrap_3 {
namespace internal {

// Class to compute the refinement depth based on bounding boxes.
template <class Kernel>
class Octree_refinement_depth
{
public:
  using Point_3 = typename Kernel::Point_3;
  using Bbox_3  = CGAL::Bbox_3;  // Assuming cells are represented as Bbox_3 objects

  // Default constructor
  Octree_refinement_depth() = default;

  // Constructor that accepts the cells (bounding boxes)
  explicit Octree_refinement_depth(const std::vector<Bbox_3>& cells)
    : m_cells(cells)
  {}

  // Function to set the cells (bounding boxes) if needed
  void set_cells(const std::vector<Bbox_3>& cells)
  {
    m_cells = cells;
  }

  // Function to return the refinement depth of the cell containing the point
  int nearest_refinement_depth(const Point_3& q) const
  {
    for (std::size_t i = 0; i < m_cells.size(); ++i)
    {
      // Check if point q lies inside the bounding box by comparing coordinates
      if (q.x() >= m_cells[i].xmin() && q.x() <= m_cells[i].xmax() &&
          q.y() >= m_cells[i].ymin() && q.y() <= m_cells[i].ymax() &&
          q.z() >= m_cells[i].zmin() && q.z() <= m_cells[i].zmax())  // Check for inclusion in Bbox_3
      {
        return get_refinement_depth(i);  // Return the refinement depth for the cell
      }
    }

    return 1;  // Default refinement depth if point is not inside any cell
  }

private:
  // Function to get the refinement depth for a given cell (based on index)
  // You may need to implement this based on your own criteria.
  int get_refinement_depth(std::size_t cell_index) const
  {
    // For now, we return a mock refinement depth based on the cell index.
    // Implement your refinement depth logic here.
    return static_cast<int>(cell_index + 1);  // Example: depth is just the index + 1
  }

private:
  std::vector<Bbox_3> m_cells;  // Cells (bounding boxes) representing the octree structure
};

} // namespace internal
} // namespace Alpha_wrap_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_OCTREE_REFINEMENT_DEPTH_H
