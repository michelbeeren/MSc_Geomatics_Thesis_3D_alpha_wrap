// Copyright (c) 2019-2023 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL: https://github.com/CGAL/cgal/blob/v6.1/Alpha_wrap_3/include/CGAL/Alpha_wrap_3/internal/Alpha_wrap_3.h $
// $Id: include/CGAL/Alpha_wrap_3/internal/Alpha_wrap_3.h b26b07a1242 $
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//                 Cedric Portaneri,
//                 Mael Rouxel-Labbé
//                 Andreas Fabri
//                 Michael Hemmer
//
#warning "USING LOCAL MODIFIED ALPHA WRAP HEADERS"
#pragma message("USING LOCAL MODIFIED ALPHA WRAP (CLION)")
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_ALPHA_WRAP_3_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_ALPHA_WRAP_3_H

#ifdef CGAL_AW3_DEBUG_PP
 #ifndef CGAL_AW3_DEBUG
  #define CGAL_AW3_DEBUG
  #define CGAL_AW3_DEBUG_INITIALIZATION
  #define CGAL_AW3_DEBUG_STEINER_COMPUTATION
  #define CGAL_AW3_DEBUG_QUEUE
  #define CGAL_AW3_DEBUG_FACET_STATUS
  #define CGAL_AW3_DEBUG_MANIFOLDNESS
 #endif
#endif

#include <CGAL/license/Alpha_wrap_3.h>
#include <CGAL/Alpha_wrap_3/internal/Feature_size_MAT.h>
// #include <CGAL/Alpha_wrap_3/internal/Octree_refinement_depth.h>
#include <CGAL/Alpha_wrap_3/internal/is_concave.h>
#include "../../../../octree.h"

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_triangulation_cell_base_3.h>
#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_triangulation_vertex_base_3.h>
#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_AABB_geom_traits.h>
#include <CGAL/Alpha_wrap_3/internal/gate_priority_queue.h>
#include <CGAL/Alpha_wrap_3/internal/geometry_utils.h>
#include <CGAL/Alpha_wrap_3/internal/oracles.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>
#include <CGAL/Robust_weighted_circumcenter_filtered_traits_3.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Default.h>
#include <CGAL/Named_function_parameters.h>
#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
 #include <CGAL/Modifiable_priority_queue.h>
#endif
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h> // only if non-manifoldness is not treated
#include <CGAL/property_map.h>
#include <CGAL/Real_timer.h>

#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <cmath>
#include <iostream>
#include <iterator>
#include <queue>
#include <string>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

namespace {

namespace AW3i = ::CGAL::Alpha_wraps_3::internal;

} // unnamed namespace

struct Wrapping_default_visitor
{
  Wrapping_default_visitor() { }

  template <typename AlphaWrapper>
  void on_alpha_wrapping_begin(const AlphaWrapper&) { }

  template <typename AlphaWrapper>
  void on_flood_fill_begin(const AlphaWrapper&) { }

  // Whether the flood filling process should continue
  template <typename Wrapper>
  constexpr bool go_further(const Wrapper&) { return true; }

  template <typename AlphaWrapper, typename Gate>
  void before_facet_treatment(const AlphaWrapper&, const Gate&) { }

  template <typename Wrapper, typename Point>
  void before_Steiner_point_insertion(const Wrapper&, const Point&) { }

  template <typename Wrapper, typename VertexHandle>
  void after_Steiner_point_insertion(const Wrapper&, VertexHandle) { }

  template <typename AlphaWrapper>
  void on_flood_fill_end(const AlphaWrapper&) { }

  template <typename AlphaWrapper>
  void on_alpha_wrapping_end(const AlphaWrapper&) { };
};

template <typename Oracle_,
          typename Triangulation_ = CGAL::Default>
class Alpha_wrapper_3
{
  using Oracle = Oracle_;

  // Triangulation
  using Base_GT = typename Oracle::Geom_traits;
  using Default_Gt = CGAL::Robust_circumcenter_filtered_traits_3<Base_GT>;

  using Default_Vb = Alpha_wrap_triangulation_vertex_base_3<Default_Gt>;
  using Default_Cb = Alpha_wrap_triangulation_cell_base_3<Default_Gt>;
  using Default_Cbt = Cell_base_with_timestamp<Default_Cb>; // for determinism
  using Default_Tds = CGAL::Triangulation_data_structure_3<Default_Vb, Default_Cbt>;
  using Default_Triangulation = CGAL::Delaunay_triangulation_3<Default_Gt, Default_Tds, Fast_location>;

public:
  using Triangulation = typename Default::Get<Triangulation_, Default_Triangulation>::type;

  // Use the geom traits from the triangulation, and trust the (advanced) user that provided it
  using Geom_traits = typename Triangulation::Geom_traits;

private:
  using Cell_handle = typename Triangulation::Cell_handle;
  using Facet = typename Triangulation::Facet;
  using Vertex_handle = typename Triangulation::Vertex_handle;
  using Locate_type = typename Triangulation::Locate_type;

  using Gate = internal::Gate<Triangulation>;

  // A sorted queue is a priority queue sorted by circumradius, and is experimentally significantly
  // slower. However, intermediate results are usable: at each point of the algorithm, the wrap
  // has a somewhat uniform mesh element size.
  //
  // An unsorted queue is a LIFO queue, and is experimentally much faster (~35%),
  // but intermediate results are not useful: a LIFO carves deep before than wide,
  // and can thus for example leave scaffolding faces till almost the end of the refinement.
#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
  using Alpha_PQ = Modifiable_priority_queue<Gate, Less_gate, Gate_ID_PM<Triangulation>, CGAL_BOOST_PAIRING_HEAP>;
#else
  using Alpha_PQ = std::stack<Gate>;
#endif

  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;
  using Ball_3 = typename Geom_traits::Ball_3;
  using Iso_cuboid_3 = typename Geom_traits::Iso_cuboid_3;

  using SC = Simple_cartesian<double>;
  using SC_Point_3 = SC::Point_3;
  using SC_Vector_3 = SC::Vector_3;
  using SC_Iso_cuboid_3 = SC::Iso_cuboid_3;
  using SC2GT = Cartesian_converter<SC, Geom_traits>;

  using Seeds = std::vector<Point_3>;

protected:
  Oracle m_oracle;
  SC_Iso_cuboid_3 m_bbox;
  double m_xmid = 0.0, m_ymid = 0.0, m_zmid = 0.0;

  // MAT support (optional)
  bool m_use_mat = false;
  std::string m_mat_path;
  mutable std::unique_ptr<CGAL::Alpha_wrap_3::internal::Feature_size_MAT<Geom_traits>> m_fs_mat_ptr;


  // Octree support (optional)
  bool m_use_octree = false;
  std::vector<CGAL::Bbox_3> m_cells;  // New member to store the octree path
  mutable std::unique_ptr<CGAL::Alpha_wrap_3::internal::Octree_refinement_depth<Geom_traits>> m_octree_ptr;

  FT m_alpha = FT(-1), m_sq_alpha = FT(-1);
  FT m_offset = FT(-1), m_sq_offset = FT(-1);

  Seeds m_seeds;

  Triangulation m_tr;

  Alpha_PQ m_queue;

public:
  Alpha_wrapper_3()
#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
      // '4096' is an arbitrary, not-too-small value for the largest ID in queue initialization
    : m_queue(4096)
#endif
  {
    // Due to the Steiner point computation being a dichotomy, the algorithm is inherently inexact
    // and passing exact kernels is explicitly disabled to ensure no misunderstanding.
    static_assert(std::is_floating_point<FT>::value);
  }

  Alpha_wrapper_3(const Oracle& oracle)
    :
      m_oracle(oracle),
      m_tr(Geom_traits(oracle.geom_traits()))
#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
      , m_queue(4096)
#endif
  {
    static_assert(std::is_floating_point<FT>::value);
  }

public:
  const Geom_traits& geom_traits() const { return m_tr.geom_traits(); }
  Oracle& oracle() { return m_oracle; }
  const Oracle& oracle() const { return m_oracle; }
  Triangulation& triangulation() { return m_tr; }
  const Triangulation& triangulation() const { return m_tr; }
  const Alpha_PQ& queue() const { return m_queue; }

  double default_alpha() const
  {
    const Bbox_3 bbox = m_oracle.bbox();
    const double diag_length = std::sqrt(square(bbox.xmax() - bbox.xmin()) +
                                         square(bbox.ymax() - bbox.ymin()) +
                                         square(bbox.zmax() - bbox.zmin()));

    return diag_length / 20.;
  }

private:
  const Point_3& circumcenter(const Cell_handle c) const
  {
    // We only cross an infinite facet once, so this isn't going to be recomputed many times
    if(m_tr.is_infinite(c))
    {
      const int inf_index = c->index(m_tr.infinite_vertex());
      c->set_circumcenter(
            geom_traits().construct_circumcenter_3_object()(m_tr.point(c, (inf_index+1)&3),
                                                            m_tr.point(c, (inf_index+2)&3),
                                                            m_tr.point(c, (inf_index+3)&3)));
    }

    return c->circumcenter(geom_traits());
  }

public:
  void set_mat_path(const std::string& path)
  {
    m_mat_path = path;
    m_use_mat  = !m_mat_path.empty();

    if(m_use_mat)
    {
      // always recreate (cheap enough compared to wrapping)
      m_fs_mat_ptr = std::make_unique<
        CGAL::Alpha_wrap_3::internal::Feature_size_MAT<Geom_traits>
      >(m_mat_path);
    }
    else
    {
      m_fs_mat_ptr.reset();
    }
  }

  void set_octree(const std::vector<CGAL::Bbox_3>& cells)
  {
    m_cells = cells; // Store the Bboxes (these should be your depth-9 cubes)
    m_use_octree = !m_cells.empty();

    if (m_use_octree) {
      // Pass the vector of Bboxes directly to your new Kd-tree based class
      m_octree_ptr = std::make_unique<CGAL::Alpha_wrap_3::internal::Octree_refinement_depth<Geom_traits>>(m_cells);
    } else {
      m_octree_ptr.reset();
    }
  }



  template <typename OutputMesh,
            typename InputNamedParameters = parameters::Default_named_parameters,
            typename OutputNamedParameters = parameters::Default_named_parameters>
  // -----------------actual 3D alpha wrap algorithm-----------------------------
  void operator()(const double alpha, // = default_alpha()
                  const double offset, // = alpha / 30.
                  OutputMesh& output_mesh,
                  const InputNamedParameters& in_np = parameters::default_values(),
                  const OutputNamedParameters& out_np = parameters::default_values())
  {
    namespace PMP = Polygon_mesh_processing;

    using parameters::get_parameter;
    using parameters::get_parameter_reference;
    using parameters::choose_parameter;

    // how vertex coordinats are stored in output_mesh
    using OVPM = typename CGAL::GetVertexPointMap<OutputMesh, OutputNamedParameters>::type;
    OVPM ovpm = choose_parameter(get_parameter(out_np, internal_np::vertex_point),
                                 get_property_map(vertex_point, output_mesh));

    // receives callbacks and is used for progress reporting, interruption, logging, ...
    using Visitor = typename internal_np::Lookup_named_param_def<
                      internal_np::visitor_t,
                      InputNamedParameters,
                      Wrapping_default_visitor // default
                    >::reference;

    Wrapping_default_visitor default_visitor;
    Visitor visitor = choose_parameter(get_parameter_reference(in_np, internal_np::visitor), default_visitor);

    // Points used to create initial cavities
    m_seeds = choose_parameter(get_parameter_reference(in_np, internal_np::seed_points), Seeds());

    // Whether or not some cells should be reflagged as "inside" after the refinement+carving loop
    // as ended, as to ensure that the result is not only combinatorially manifold, but also
    // geometrically manifold.
    //
    // -- Warning --
    // Manifoldness postprocessing will be performed even if the wrapping is interrupted (and
    // this option is enabled).
    const bool do_enforce_manifoldness = choose_parameter(get_parameter(in_np, internal_np::do_enforce_manifoldness), true);

    // Whether to keep pockets of "outside" cells that are not connected to the exterior (or to the
    // initial cavities, if used).
    //
    // -- Warning --
    // Pockets of "outside" cells will be purged even if the wrapping is interrupted (and
    // this option is enabled).
    const bool keep_inner_ccs = choose_parameter(get_parameter(in_np, internal_np::keep_inner_connected_components), false);

    // This parameter enables avoiding recomputing the triangulation from scratch when wrapping
    // the same input for multiple values of alpha (and typically the same offset values).
    //
    // -- Warning --
    // If this is enabled, the 3D triangulation will NOT be re-initialized at launch.
    // This means that the triangulation is NOT cleared, even if:
    // - you use an alpha value that is greater than what was used in a previous run; you will
    //   obtain the same result as the last run.
    // - you use a different offset value between runs, you might then get points that are not
    //   on the offset surface corresponding to that corresponding to the latter offset value.
    const bool refining = choose_parameter(get_parameter(in_np, internal_np::refine_triangulation), false);

#ifdef CGAL_AW3_TIMER
    CGAL::Real_timer t;
    t.start();
#endif

    // --------------3D alpha wrap algorithm actually begins----------
    visitor.on_alpha_wrapping_begin(*this);

    // initialise: compute offset surface, insert points, builds or updates 3D triangulation
    if(!initialize(alpha, offset, refining))
      return;

#ifdef CGAL_AW3_TIMER
    t.stop();
    std::cout << "Initialization took: " << t.time() << " s." << std::endl;
    t.reset();
    t.start();
#endif

#ifdef CGAL_AW3_DEBUG_DUMP_INTERMEDIATE_WRAPS
    dump_triangulation_faces("starting_wrap.off", true /*only_boundary_faces*/);
#endif

    // tetrahedracells are labeled inside/outside (based on alpha treshhold and intersection?) --> so probably need to change parts here
    alpha_flood_fill(visitor);

#ifdef CGAL_AW3_DEBUG_DUMP_INTERMEDIATE_WRAPS
    dump_triangulation_faces("flood_filled_wrap.off", true /*only_boundary_faces*/);
#endif

#ifdef CGAL_AW3_TIMER
    t.stop();
    std::cout << "Flood filling took: " << t.time() << " s." << std::endl;
    t.reset();
    t.start();
#endif

    if(do_enforce_manifoldness)
    {
      make_manifold();

#ifdef CGAL_AW3_DEBUG_DUMP_INTERMEDIATE_WRAPS
      dump_triangulation_faces("manifold_wrap.off", true /*only_boundary_faces*/);
#endif

#ifdef CGAL_AW3_TIMER
      t.stop();
      std::cout << "Manifoldness postprocessing took: " << t.time() << " s." << std::endl;
      t.reset();
      t.start();
#endif
    }

    if(!keep_inner_ccs)
    {
      // We could purge *before* manifold enforcement, but making the mesh manifold is
      // very cheap in most cases, so it is better to keep the code simpler.
      purge_inner_connected_components();

#ifdef CGAL_AW3_DEBUG_DUMP_INTERMEDIATE_WRAPS
      dump_triangulation_faces("purged_wrap.off", true /*only_boundary_faces*/);
#endif
    }

    extract_surface(output_mesh, ovpm, !do_enforce_manifoldness);

#ifdef CGAL_AW3_TIMER
    t.stop();
    std::cout << "Surface extraction took: " << t.time() << " s." << std::endl;
#endif

#ifdef CGAL_AW3_DEBUG
    std::cout << "Alpha wrap vertices:  " << vertices(output_mesh).size() << std::endl;
    std::cout << "Alpha wrap faces:     " << faces(output_mesh).size() << std::endl;

 #ifdef CGAL_AW3_DEBUG_DUMP_INTERMEDIATE_WRAPS
    IO::write_polygon_mesh("final_wrap.off", output_mesh, CGAL::parameters::stream_precision(17));
    dump_triangulation_faces("final_tr.off", false /*only_boundary_faces*/);
 #endif
#endif

    visitor.on_alpha_wrapping_end(*this);
    // std::cout << "ooh managed to modify the alpha wrap algorithm :)" << std::endl;
  }

  // -------------------------------------------------------------------------------------------
  // Compatibility overload (OLD): keeps CGAL's existing alpha_wrap_3.h wrappers compiling.
  // Forwards to the new overload with empty MAT_path => behavior A.
  // -------------------------------------------------------------------------------------------
  template <typename OutputMesh>
  void operator()(const double alpha,
                  OutputMesh& output_mesh)
  {
    return operator()(alpha, alpha / 30. /*offset*/, output_mesh);
  }

  // Convenience overloads
  template <typename OutputMesh>
  void operator()(OutputMesh& output_mesh)
  {
    return operator()(default_alpha(), output_mesh);
  }



  // This function is public only because it is used in the tests
  SC_Iso_cuboid_3 construct_bbox(const double offset)
  {
    // Input axis-aligned bounding box
    SC_Iso_cuboid_3 bbox = m_oracle.bbox();
    const SC_Point_3 bbox_centroid = midpoint((bbox.min)(), (bbox.max)());

    // Scale a bit to create the initial points not too close to the input
    double scaling = 1.2;
    CGAL::Aff_transformation_3<SC> scale(SCALING, scaling);
    bbox = SC_Iso_cuboid_3(scale.transform((bbox.min)()), scale.transform((bbox.max)()));

    // Translate bbox back to initial centroid
    const SC_Point_3 bbox_transformed_centroid = midpoint((bbox.min)(), (bbox.max)());
    const SC_Vector_3 diff_centroid = bbox_centroid - bbox_transformed_centroid;
    CGAL::Aff_transformation_3<SC> centroid_translate(TRANSLATION, diff_centroid);
    bbox = SC_Iso_cuboid_3(centroid_translate.transform((bbox.min)()),
                           centroid_translate.transform((bbox.max)()));

    // Add the offset
    SC_Vector_3 offset_ext = std::sqrt(3.) * offset * SC_Vector_3(1, 1, 1);
    CGAL::Aff_transformation_3<SC> translate_m(TRANSLATION, - offset_ext);
    CGAL::Aff_transformation_3<SC> translate_M(TRANSLATION,   offset_ext);
    bbox = SC_Iso_cuboid_3(translate_m.transform((bbox.min)()), translate_M.transform((bbox.max)()));

    return bbox;
  }

private:
  // The distinction between inside boundary and outside boundary is the presence of cells
  // being flagged for manifoldness: inside boundary considers those outside, and outside
  // boundary considers them inside.
  bool is_on_inside_boundary(Cell_handle ch, Cell_handle nh) const
  {
    return (ch->is_inside() != nh->is_inside()); // one is "inside", the other is not
  }

  bool is_on_outside_boundary(Cell_handle ch, Cell_handle nh) const
  {
    return (ch->is_outside() != nh->is_outside()); // one is "outside", the other is not
  }

private:
  // Adjust the bbox & insert its corners to construct the starting triangulation
  void insert_bbox_corners()
  {
    m_bbox = construct_bbox(CGAL::to_double(m_offset));

    m_xmid = 0.5 * (m_bbox.xmin() + m_bbox.xmax());
    m_ymid = 0.5 * (m_bbox.ymin() + m_bbox.ymax());
    m_zmid = 0.5 * (m_bbox.zmin() + m_bbox.zmax());

#ifdef CGAL_AW3_DEBUG_INITIALIZATION
    std::cout << "Insert Bbox vertices" << std::endl;
#endif

    // insert in dt the eight corner vertices of the input loose bounding box
    for(int i=0; i<8; ++i)
    {
      const Point_3 bp = SC2GT()(m_bbox.vertex(i));
      Vertex_handle bv = m_tr.insert(bp);
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
      std::cout << "\t" << bp << std::endl;
#endif
      bv->type() = AW3i::Vertex_type:: BBOX_VERTEX;
    }
  }

  // Two criteria:
  // - Cells that are intersecting the input are inside
  // - Cells whose circumcenter is in the offset volume are inside: this is because
  // we need to have outside cell circumcenters outside of the volume to ensure
  // that the refinement point is separated from the existing point set.
  Cell_label cavity_cell_label(const Cell_handle ch)
  {
    CGAL_precondition(!m_tr.is_infinite(ch));

    const Tetrahedron_with_outside_info<Geom_traits> tet(ch, geom_traits());
    if(m_oracle.do_intersect(tet))
      return Cell_label::INSIDE;

    const Point_3& ch_cc = circumcenter(ch);
    typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
    const Ball_3 ch_cc_offset_ball = ball(ch_cc, m_sq_offset);
    const bool is_cc_in_offset = m_oracle.do_intersect(ch_cc_offset_ball);

    return is_cc_in_offset ? Cell_label::INSIDE : Cell_label::OUTSIDE;
  }

  // Create a cavity using seeds rather than starting from the infinity.
  //
  // For each seed, insert the seeds and its neighboring vertices on a regular icosahedron.
  // The idea behind an icosahedron rather than e.g. a simple tetrahedron is twofold:
  // - Improve the odds of both starting the algorithm, but also of not immediately stopping,
  //   which could happen if a conflict zone included all the initial cavity when using a simple cavity.
  // - Base the cavity diameter on alpha, and allow its cells to intersect the input. If a single
  //   tetrahedron is used, one is forced to make it small-enough such that it does not intersect
  //   the input, and then it could be forced to be smaller than 'alpha' and the algorithm cannot start,
  //   for example if the seed is close to the offset. If we have many tetrahedra, this is far less
  //   likely to happen.
  //
  // For an edge length a, the radius of an icosahedron is
  //   r = 0.5 * a * sqrt(phi * sqrt(5)), phi being the golden ratio
  // which yields
  //   r = a * sin(2pi / 5) =~ 0.95105651629515353 * a
  // Faces of an icosahedron are equilateral triangles with size a and circumradius a / sqrt(3)
  // Since we want faces of the icosahedron to be traversable, we want a such that
  //   a / sqrt(3) > alpha
  // Hence r such that
  //   r / (sqrt(3) * sin(2pi/5)) > alpha
  //   r > alpha * sqrt(3) * sin(2pi / 5) ~= 1.6472782070926637 * alpha
  //
  // Furthermore, the triangles between edges of the icosahedron and the center of the icosahedron
  // are not equilateral triangles since a is slightly bigger than r. They are
  // slightly flattened isocele triangles with base 'a' and the circumradius is smaller.
  // The circumradius is
  //   r_iso = r² / (2 * h) = r² / (2 * sqrt(r² - (a / 2)²))
  // Since r = a * sin(2pi / 5)
  //  r_iso = r² / (2 * sqrt(r² - (r / 2 * sin(2pi / 5))²))
  //        = r / (2 * sqrt(1 - 1 / (2 * sin(2pi / 5))²))
  // So we want r_iso > alpha, i.e.
  //   r > (2 * sqrt(1 - 1 / (2 * sin(2pi / 5))²)) * alpha ~= 1.7013016167040798 * alpha
  //
  // Another way is to simply make faces incident to the seed always traversable, and then
  // we only have to ensure faces opposite of the seed are traversable (i.e., radius ~= 1.65 * alpha)
  bool initialize_with_cavities()
  {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
    std::cout << "> Dig cavities" << std::endl;
    std::cout << m_seeds.size() << " seed(s)" << std::endl;
#endif

    CGAL_precondition(!m_seeds.empty());

    // Get a double value approximating the scaling factors
//    std::cout << sqrt(3) * sin(2pi / 5) << std::endl;
//    std::cout << (2. * std::sqrt(1. - 1. / square(2 * std::sin(2 * CGAL_PI / 5)))) << std::endl;

    Iso_cuboid_3 bbox = SC2GT()(m_bbox);

    std::vector<Vertex_handle> seed_vs;
    for(const Point_3& seed_p : m_seeds)
    {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
      std::cout << "Initialize from seed " << seed_p << std::endl;
#endif

      if(bbox.has_on_unbounded_side(seed_p))
      {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
        std::cerr << "Warning: seed " << seed_p << " is outside the bounding box" << std::endl;
#endif
        continue;
      }

      // get the closest point on the input
      const Point_3 closest_pt = m_oracle.closest_point(seed_p);
      const FT sq_d_to_closest = geom_traits().compute_squared_distance_3_object()(seed_p, closest_pt);

      if(sq_d_to_closest < m_sq_offset)
      {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
        std::cerr << "Warning: seed " << seed_p << " is in the offset" << std::endl;
#endif
        continue;
      }

      // Mark the seeds and icosahedron vertices as "scaffolding" vertices such that the facets
      // incident to these vertices are always traversable regardless of their circumcenter.
      // This is done because otherwise some cavities can appear on the mesh: non-traversable facets
      // with two vertices on the offset, and the third being a deeper inside seed / ico_seed.
      // Making them traversable will cause more refinement than "alpha", but they will eventually
      // not appear anymore in the inside/outside boundary and the surface will look smoother.
      //
      // This problem only appears when the seed and icosahedron vertices are close to the offset surface,
      // which usually happens for large alpha values.

      Vertex_handle seed_v = m_tr.insert(seed_p);
      seed_v->type() = AW3i::Vertex_type:: SEED_VERTEX;
      seed_vs.push_back(seed_v);

      // Icosahedron vertices (see also BGL::make_icosahedron())
      const Point_3 center = seed_p;
      const FT radius = 1.65 * m_alpha;
      const FT phi = (FT(1) + approximate_sqrt(FT(5))) / FT(2);
      const FT t = radius / approximate_sqrt(1 + square(phi));
      const FT t_phi = t * phi;

      std::array<Point_3, 12> ico_ps =
      {
        Point_3(center.x(), center.y() + t, center.z() + t_phi),
        Point_3(center.x(), center.y() + t, center.z() - t_phi),
        Point_3(center.x(), center.y() - t, center.z() + t_phi),
        Point_3(center.x(), center.y() - t, center.z() - t_phi),

        Point_3(center.x() + t, center.y() + t_phi, center.z()),
        Point_3(center.x() + t, center.y() - t_phi, center.z()),
        Point_3(center.x() - t, center.y() + t_phi, center.z()),
        Point_3(center.x() - t, center.y() - t_phi, center.z()),

        Point_3(center.x() + t_phi, center.y(), center.z() + t),
        Point_3(center.x() + t_phi, center.y(), center.z() - t),
        Point_3(center.x() - t_phi, center.y(), center.z() + t),
        Point_3(center.x() - t_phi, center.y(), center.z() - t)
      };

      for(const Point_3& seed_neighbor_p : ico_ps)
      {
#ifdef CGAL_AW3_DEBUG_PP
        std::cout << seed_neighbor_p << std::endl;
#endif
        if(bbox.has_on_unbounded_side(seed_neighbor_p))
          continue;

        Vertex_handle ico_v = m_tr.insert(seed_neighbor_p, seed_v /*hint*/);
        ico_v->type() = AW3i::Vertex_type:: SEED_VERTEX;
      }
    }

    if(seed_vs.empty())
    {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
      std::cerr << "Error: no acceptable seed was provided" << std::endl;
#endif
      return false;
    }

#ifdef CGAL_AW3_DEBUG_INITIALIZATION
    std::cout << m_tr.number_of_vertices() - 8 /*bbox*/ << " vertice(s) due to seeds" << std::endl;
#endif

    for(Vertex_handle seed_v : seed_vs)
    {
      std::vector<Cell_handle> inc_cells;
      inc_cells.reserve(64);
      m_tr.incident_cells(seed_v, std::back_inserter(inc_cells));

      for(Cell_handle ch : inc_cells)
        ch->set_label(cavity_cell_label(ch));
    }

    // Should be cheap enough to go through the full triangulation as only seeds have been inserted
    for(Cell_handle ch : m_tr.all_cell_handles())
    {
      if(ch->is_inside())
        continue;

      // When the algorithm starts from a manually dug hole, infinite cells are initialized
      // as "inside" such that they do not appear on the boundary
      CGAL_assertion(!m_tr.is_infinite(ch));

      for(int i=0; i<4; ++i)
        push_facet(std::make_pair(ch, i));
    }

    if(m_queue.empty())
    {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
      std::cerr << "Could not initialize the algorithm with these seeds, and alpha|offset values" << std::endl;
#endif
      return false;
    }

    return true;
  }

  // tag all infinite cells "outside" and all finite cells "inside"
  // init queue with all convex hull facets
  bool initialize_from_infinity()
  {
    for(Cell_handle ch : m_tr.all_cell_handles())
    {
      if(m_tr.is_infinite(ch))
      {
        ch->set_label(Cell_label::OUTSIDE);
        const int inf_index = ch->index(m_tr.infinite_vertex());
        push_facet(std::make_pair(ch, inf_index));
      }
      else
      {
        ch->set_label(Cell_label::INSIDE);
      }
    }

    return true;
  }

  void reset_manifold_labels()
  {
    // No erase counter increment, or it will mess up with a possibly non-empty queue.
    for(Cell_handle ch : m_tr.all_cell_handles())
      if(ch->label() == Cell_label::MANIFOLD)
        ch->set_label(Cell_label::OUTSIDE);
  }

  // This function is used in the case of resumption of a previous run: m_tr is not cleared,
  // and we fill the queue with the new parameters.
  bool initialize_from_existing_triangulation()
  {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
    std::cout << "Restart from a DT of " << m_tr.number_of_cells() << " cells" << std::endl;
#endif

    for(Cell_handle ch : m_tr.all_cell_handles())
    {
      if(ch->is_inside())
        continue;

      for(int i=0; i<4; ++i)
      {
        if(ch->neighbor(i)->is_inside())
          push_facet(std::make_pair(ch, i));

      }
    }

    return true;
  }

private:
  // Manifoldness is tolerated while debugging and extracting at intermediate states
  // Not the preferred way because it uses 3*nv storage
  template <typename OutputMesh, typename OVPM>
  void extract_possibly_non_manifold_surface(OutputMesh& output_mesh,
                                             OVPM ovpm) const
  {
    namespace PMP = Polygon_mesh_processing;

#ifdef CGAL_AW3_DEBUG
    std::cout << "> Extract possibly non-manifold wrap... ()" << std::endl;
#endif

    remove_all_elements(output_mesh);

    CGAL_assertion_code(for(auto cit=m_tr.finite_cells_begin(), cend=m_tr.finite_cells_end(); cit!=cend; ++cit))
    CGAL_assertion(cit->tds_data().is_clear());

    for(auto cit=m_tr.finite_cells_begin(), cend=m_tr.finite_cells_end(); cit!=cend; ++cit)
    {
      Cell_handle seed = cit;
      if(seed->is_outside() || seed->tds_data().processed())
        continue;

      std::queue<Cell_handle> to_visit;
      to_visit.push(seed);

      std::vector<Point_3> points;
      std::vector<std::vector<size_t> > faces;
      std::size_t idx = 0;

      while(!to_visit.empty())
      {
        const Cell_handle cell = to_visit.front();
        CGAL_assertion(!cell->is_outside() && !m_tr.is_infinite(cell));

        to_visit.pop();

        if(cell->tds_data().processed())
          continue;
        cell->tds_data().mark_processed();

        for(int fid=0; fid<4; ++fid)
        {
          const Cell_handle neighbor = cell->neighbor(fid);
          if(neighbor->is_outside())
          {
            // There shouldn't be any artificial vertex on the inside/outside boundary
            // (past initialization)
//            CGAL_assertion(cell->vertex((fid + 1)&3)->type() == AW3i::Vertex_type:: DEFAULT);
//            CGAL_assertion(cell->vertex((fid + 2)&3)->type() == AW3i::Vertex_type:: DEFAULT);
//            CGAL_assertion(cell->vertex((fid + 3)&3)->type() == AW3i::Vertex_type:: DEFAULT);

            points.push_back(m_tr.point(cell, Triangulation::vertex_triple_index(fid, 0)));
            points.push_back(m_tr.point(cell, Triangulation::vertex_triple_index(fid, 1)));
            points.push_back(m_tr.point(cell, Triangulation::vertex_triple_index(fid, 2)));
            faces.push_back({idx, idx + 1, idx + 2});
            idx += 3;
          }
          else
          {
            to_visit.push(neighbor);
          }
        }
      }

      CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(faces));
      PMP::polygon_soup_to_polygon_mesh(points, faces, output_mesh,
                                        CGAL::parameters::default_values(),
                                        CGAL::parameters::vertex_point_map(ovpm));

      PMP::stitch_borders(output_mesh, CGAL::parameters::vertex_point_map(ovpm));
      CGAL_assertion(is_closed(output_mesh));
    }

    for(auto cit=m_tr.finite_cells_begin(), cend=m_tr.finite_cells_end(); cit!=cend; ++cit)
      cit->tds_data().clear();

    CGAL_postcondition(!is_empty(output_mesh));
    CGAL_postcondition(is_valid_polygon_mesh(output_mesh));
    CGAL_postcondition(is_closed(output_mesh));

    PMP::orient_to_bound_a_volume(output_mesh, CGAL::parameters::vertex_point_map(ovpm));

    collect_garbage(output_mesh);
  }

  template <typename OutputMesh, typename OVPM>
  void extract_manifold_surface(OutputMesh& output_mesh,
                                OVPM ovpm) const
  {
    namespace PMP = Polygon_mesh_processing;

#ifdef CGAL_AW3_DEBUG
    std::cout << "> Extract manifold wrap... ()" << std::endl;
#endif

    CGAL_assertion_code(for(Vertex_handle v : m_tr.finite_vertex_handles()))
    CGAL_assertion(!is_non_manifold(v));

    remove_all_elements(output_mesh);

    std::vector<Point_3> points;
    std::vector<std::array<std::size_t, 3> > faces;

    std::unordered_map<Vertex_handle, std::size_t> vertex_to_id;

    for(auto fit=m_tr.all_facets_begin(), fend=m_tr.all_facets_end(); fit!=fend; ++fit)
    {
      Facet f = *fit;
      if(!f.first->is_outside())
        f = m_tr.mirror_facet(f);

      const Cell_handle ch = f.first;
      const int s = f.second;
      const Cell_handle nh = ch->neighbor(s);
      if(!is_on_outside_boundary(ch, nh))
        continue;

      std::array<std::size_t, 3> ids;
      for(int pos=0; pos<3; ++pos)
      {
        Vertex_handle vh = ch->vertex(Triangulation::vertex_triple_index(s, pos));
        auto insertion_res = vertex_to_id.emplace(vh, vertex_to_id.size());
        if(insertion_res.second) // successful insertion, never-seen-before vertex
          points.push_back(m_tr.point(vh));

        ids[pos] = insertion_res.first->second;
      }

      faces.emplace_back(std::array<std::size_t, 3>{ids[0], ids[1], ids[2]});
    }

#ifdef CGAL_AW3_DEBUG
    std::cout << "\t" << points.size() << " points" << std::endl;
    std::cout << "\t" << faces.size() << " faces" << std::endl;
#endif

    if(faces.empty())
    {
#ifdef CGAL_AW3_DEBUG
      std::cerr << "Warning: empty wrap?..." << std::endl;
#endif
      return;
    }

    if(!PMP::is_polygon_soup_a_polygon_mesh(faces))
    {
      CGAL_warning_msg(false, "Failed to extract a manifold boundary!");
      return;
    }

    PMP::polygon_soup_to_polygon_mesh(points, faces, output_mesh,
                                      CGAL::parameters::default_values(),
                                      CGAL::parameters::vertex_point_map(ovpm));

    CGAL_postcondition(!is_empty(output_mesh));
    CGAL_postcondition(is_valid_polygon_mesh(output_mesh));
    CGAL_postcondition(is_closed(output_mesh));
    CGAL_postcondition(PMP::does_bound_a_volume(output_mesh, CGAL::parameters::vertex_point_map(ovpm)));

    PMP::orient_to_bound_a_volume(output_mesh, CGAL::parameters::vertex_point_map(ovpm));
  }

public:
  template <typename OutputMesh, typename OVPM>
  void extract_surface(OutputMesh& output_mesh,
                       OVPM ovpm,
                       const bool tolerate_non_manifoldness = false) const
  {
    if(tolerate_non_manifoldness)
      extract_possibly_non_manifold_surface(output_mesh, ovpm);
    else
      extract_manifold_surface(output_mesh, ovpm);
  }

private:
  // ooh I think I just need to change this a bit :)

  // bool is_traversable(const Facet& f) const { return less_squared_radius_of_min_empty_sphere(m_sq_alpha, f, m_tr); }

bool is_traversable(const Facet& f) const
{
  // A) Default CGAL behavior: when neither MAT nor Octree is used
  if ((!m_use_mat || !m_fs_mat_ptr) && (!m_use_octree || !m_octree_ptr))
  {
    // ---- 1. Alpha condition (compute ONCE) ----
    bool traversable =
        less_squared_radius_of_min_empty_sphere(m_sq_alpha, f, m_tr);

    if (traversable) {
      // std::cout << "face is traversible!" << std::endl;
      return traversable;
    }
    // else {
    //   // std::cout << "face is not traversible!" << std::endl;
    //   // ---- 2. Compute face midpoint ----
    //   const Cell_handle ch = f.first;
    //   const int i = f.second;
    //
    //   const Point_3& p0 = ch->vertex((i + 1) & 3)->point();
    //   const Point_3& p1 = ch->vertex((i + 2) & 3)->point();
    //   const Point_3& p2 = ch->vertex((i + 3) & 3)->point();
    //   // typename Geom_traits::Line_3 l12(p1, p2);
    //   // typename Geom_traits::Line_3 l02(p0, p2);
    //   // typename Geom_traits::Line_3 l01(p0, p1);
    //   //
    //   // double min_size_tr = 0.5*m_offset;
    //   //
    //   // double d0 = std::sqrt(CGAL::to_double(CGAL::squared_distance(p0, l12)));
    //   // double d1 = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, l02)));
    //   // double d2 = std::sqrt(CGAL::to_double(CGAL::squared_distance(p2, l01)));
    //   // if (d0 < min_size_tr || d1 < min_size_tr || d2 < min_size_tr) {
    //   //   return false;
    //   // }
    //
    //   // std::cout << "distances to opposite lines are:" << d0 << ", " << d1 << ", " << d2 << std::endl;
    //
    //   const Point_3 f_mid(
    //       (p0.x() + p1.x() + p2.x()) / 3.0,
    //       (p0.y() + p1.y() + p2.y()) / 3.0,
    //       (p0.z() + p1.z() + p2.z()) / 3.0
    //   );
    //
    //   // ---- 3. Distance to input surface via Oracle ----
    //   const Point_3 closest_pt = m_oracle.closest_point(f_mid);
    //
    //   const FT sq_d =
    //       geom_traits().compute_squared_distance_3_object()(f_mid, closest_pt);
    //
    //   const double dist_input =
    //       std::sqrt(CGAL::to_double(sq_d));
    //
    //   // ---- 4. Distance to offset surface (optional) ----
    //   const double dist_offset =
    //       dist_input - CGAL::to_double(m_offset);
    //
    //   const double tolerance = 1.2 * m_offset;
    //
    //   if (dist_input > tolerance) {
    //     // std::cout << "tolerance = " << tolerance << std::endl;
    //     // std::cout << "Distance midpoint -> input surface: "<< dist_input << std::endl;
    //     // std::cout << "Distance midpoint -> offset surface: " << dist_offset << std::endl;
    //     // std::cout << "Oof distance of triangle far from offset!" << std::endl;
    //     return true;
    //   }
    // }
    return traversable;
  }

    // B) MAT-based behavior: when MAT is used
    if ((m_use_mat || m_fs_mat_ptr) && (!m_use_octree || !m_octree_ptr)) {
        const Cell_handle ch = f.first;
        const int i = f.second;

        const Point_3& p0 = ch->vertex((i + 1) & 3)->point();
        const Point_3& p1 = ch->vertex((i + 2) & 3)->point();
        const Point_3& p2 = ch->vertex((i + 3) & 3)->point();

        const Point_3 f_mid(
            (p0.x() + p1.x() + p2.x()) / 3.0,
            (p0.y() + p1.y() + p2.y()) / 3.0,
            (p0.z() + p1.z() + p2.z()) / 3.0
        );

        const double fs = m_fs_mat_ptr->nearest_feature_size(f_mid);
        const double fs_safe = (fs > 0.0) ? fs : 1e-12;

        const double local_alpha = CGAL::to_double(m_alpha) * std::pow(fs_safe, 0.85); // alpha*3 for ref_depth=1, alpha*0.3 for ref_depth=9
        const FT local_sq_alpha = FT(local_alpha * local_alpha);
      // std:: cout << "does function come here?" << std::endl;

        return less_squared_radius_of_min_empty_sphere(local_sq_alpha, f, m_tr);
    }

    // C) Octree-based behavior: when Octree is used but no MAT
    if ((!m_use_mat || !m_fs_mat_ptr) && (m_use_octree || m_octree_ptr)) {
        const Cell_handle ch = f.first;
        const int i = f.second;

        const Point_3& p0 = ch->vertex((i + 1) & 3)->point();
        const Point_3& p1 = ch->vertex((i + 2) & 3)->point();
        const Point_3& p2 = ch->vertex((i + 3) & 3)->point();

        const Point_3 f_mid(
            (p0.x() + p1.x() + p2.x()) / 3.0,
            (p0.y() + p1.y() + p2.y()) / 3.0,
            (p0.z() + p1.z() + p2.z()) / 3.0
        );
      // std::cout << "f_mid: " << f_mid.x() << ", " << f_mid.y() << ", " << f_mid.z() << std::endl;
      double alpha_val = CGAL::to_double(m_alpha);
      double sq_dist = m_octree_ptr->squared_dist_to_finest(f_mid);
      double dist = std::sqrt(sq_dist);

      double edge_len = m_octree_ptr->finest_edge_length();
      double d_min = 0.87 * edge_len;
      double d_max = 15*alpha_val + d_min; // set length from concave area till a alpha ball of size 10 * m_alpha can be used

      double local_alpha;

      if (dist <= d_min) {
        // Region 1: Inside or very close to the concave feature
        local_alpha = alpha_val;
      }
      else if (dist >= d_max) {
        // Region 2: Far enough away to use full alpha
        local_alpha = 10 * alpha_val;
      }
      else {
        // Region 3: Linear transition from 0.1*alpha to 1.0*alpha
        // Interpolation factor (t) goes from 0.0 to 1.0
        double t = (dist - d_min) / (d_max - d_min);
        local_alpha = (alpha_val) + 10 * t * (alpha_val - 0.1 * alpha_val);
      }

      const FT local_sq_alpha = FT(local_alpha * local_alpha);
      return less_squared_radius_of_min_empty_sphere(local_sq_alpha, f, m_tr);
    }

    // D) Default case: if no MAT or Octree is provided, use the original logic.
    return less_squared_radius_of_min_empty_sphere(m_sq_alpha, f, m_tr);
}


//   bool compute_steiner_point(const Cell_handle ch,
//                              const Cell_handle neighbor,
//                              Point_3& steiner_point) const
//   {
//     CGAL_precondition(!m_tr.is_infinite(neighbor));
//
//     typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
//     typename Geom_traits::Construct_vector_3 vector = geom_traits().construct_vector_3_object();
//     typename Geom_traits::Construct_translated_point_3 translate = geom_traits().construct_translated_point_3_object();
//     typename Geom_traits::Construct_scaled_vector_3 scale = geom_traits().construct_scaled_vector_3_object();
//
//     const Point_3& neighbor_cc = circumcenter(neighbor);
//     const Ball_3 neighbor_cc_offset_ball = ball(neighbor_cc, m_sq_offset);
//     const bool is_neighbor_cc_in_offset = m_oracle.do_intersect(neighbor_cc_offset_ball);
//
// #ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
//     std::cout << "Compute_steiner_point(" << &*ch << ", " << &*neighbor << ")" << std::endl;
//
//     const Point_3& chc = circumcenter(ch);
//     std::cout << "CH" << std::endl;
//     std::cout << "\t" << ch->vertex(0)->point() << std::endl;
//     std::cout << "\t" << ch->vertex(1)->point() << std::endl;
//     std::cout << "\t" << ch->vertex(2)->point() << std::endl;
//     std::cout << "\t" << ch->vertex(3)->point() << std::endl;
//     std::cout << "cc: " << chc << std::endl;
//     std::cout << "CC Distance to input: " << CGAL::squared_distance(chc, m_oracle.closest_point(chc)) << std::endl;
//
//     std::cout << "NCH" << std::endl;
//     std::cout << "\t" << neighbor->vertex(0)->point() << std::endl;
//     std::cout << "\t" << neighbor->vertex(1)->point() << std::endl;
//     std::cout << "\t" << neighbor->vertex(2)->point() << std::endl;
//     std::cout << "\t" << neighbor->vertex(3)->point() << std::endl;
//     std::cout << "ncc: " << neighbor_cc << std::endl;
//     std::cout << "NC Distance to input: " << CGAL::squared_distance(neighbor_cc, m_oracle.closest_point(neighbor_cc)) << std::endl;
//     std::cout << "is_neighbor_cc_in_offset = " << is_neighbor_cc_in_offset << std::endl;
//
//     std::cout << "squared offset " << m_sq_offset << std::endl;
// #endif
//
//     // ch's circumcenter should not be within the offset volume
//     CGAL_assertion_code(const Point_3& ch_cc = circumcenter(ch);)
//     CGAL_assertion_code(const Ball_3 ch_cc_offset_ball = ball(ch_cc, m_sq_offset);)
//     CGAL_assertion(!m_oracle.do_intersect(ch_cc_offset_ball));
//
//     if(is_neighbor_cc_in_offset)
//     {
//       const Point_3& ch_cc = circumcenter(ch);
//
//       // std::cout << "ch_cc: " << ch_cc << std::endl;
//       // std::cout << "neighbor_cc: " << neighbor_cc << std::endl;
//       // If the voronoi edge intersects the offset, the steiner point is the first intersection
// //       if(m_oracle.first_intersection(ch_cc, neighbor_cc, steiner_point, m_offset))
// //       {
// // #ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
// //         std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
// // #endif
// //         return true;
// //       }
//     }
//
//     Tetrahedron_with_outside_info<Geom_traits> tet(neighbor, geom_traits());
//     if(m_oracle.do_intersect(tet))
//     {
//       const Point_3& ch_cc = circumcenter(ch);
//       // steiner point is the closest point on input from cell centroid with offset
//       const Point_3 closest_pt = m_oracle.closest_point(ch_cc);
//       CGAL_assertion(closest_pt != ch_cc);
//
//
//       Vector_3 unit = vector(closest_pt, ch_cc); //       Vector_3 unit = vector(closest_pt, ch_cc); produces less points, but outpuit is invalid :(
//
//       // PMP::internal::normalize() requires sqrt()
//       unit = scale(unit, m_offset / approximate_sqrt(geom_traits().compute_squared_length_3_object()(unit)));
//       steiner_point = translate(closest_pt, unit);
//
// #ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
//       std::cout << "Steiner found through neighboring tet intersecting the input: " << steiner_point << std::endl;
//       std::cout << "Closest point: " << closest_pt << std::endl;
//       std::cout << "Direction: " << vector(closest_pt, neighbor_cc) << std::endl;
// #endif
//
//       return true;
//     }
//
// #ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
//     std::cout << "No Steiner point" << std::endl;
// #endif
//
//     return false;
//   }

//     bool compute_steiner_point(const Cell_handle ch,
//                              const Cell_handle neighbor,
//                              Point_3& steiner_point) const
//   {
//     CGAL_precondition(!m_tr.is_infinite(neighbor));
//
//     typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
//     typename Geom_traits::Construct_vector_3 vector = geom_traits().construct_vector_3_object();
//     typename Geom_traits::Construct_translated_point_3 translate = geom_traits().construct_translated_point_3_object();
//     typename Geom_traits::Construct_scaled_vector_3 scale = geom_traits().construct_scaled_vector_3_object();
//     typename Geom_traits::Compute_squared_length_3 sq_length = geom_traits().compute_squared_length_3_object();
//     typename Geom_traits::Compute_scalar_product_3 dot = geom_traits().compute_scalar_product_3_object();
//
//     const Point_3& neighbor_cc = circumcenter(neighbor);
//     const Ball_3 neighbor_cc_offset_ball = ball(neighbor_cc, m_sq_offset);
//     const bool is_neighbor_cc_in_offset = m_oracle.do_intersect(neighbor_cc_offset_ball);
//
// #ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
//     std::cout << "Compute_steiner_point(" << &*ch << ", " << &*neighbor << ")" << std::endl;
//
//     const Point_3& chc = circumcenter(ch);
//     std::cout << "CH" << std::endl;
//     std::cout << "\t" << ch->vertex(0)->point() << std::endl;
//     std::cout << "\t" << ch->vertex(1)->point() << std::endl;
//     std::cout << "\t" << ch->vertex(2)->point() << std::endl;
//     std::cout << "\t" << ch->vertex(3)->point() << std::endl;
//     std::cout << "cc: " << chc << std::endl;
//     std::cout << "CC Distance to input: " << CGAL::squared_distance(chc, m_oracle.closest_point(chc)) << std::endl;
//
//     std::cout << "NCH" << std::endl;
//     std::cout << "\t" << neighbor->vertex(0)->point() << std::endl;
//     std::cout << "\t" << neighbor->vertex(1)->point() << std::endl;
//     std::cout << "\t" << neighbor->vertex(2)->point() << std::endl;
//     std::cout << "\t" << neighbor->vertex(3)->point() << std::endl;
//     std::cout << "ncc: " << neighbor_cc << std::endl;
//     std::cout << "NC Distance to input: " << CGAL::squared_distance(neighbor_cc, m_oracle.closest_point(neighbor_cc)) << std::endl;
//     std::cout << "is_neighbor_cc_in_offset = " << is_neighbor_cc_in_offset << std::endl;
//
//     std::cout << "squared offset " << m_sq_offset << std::endl;
// #endif
//
//     // ch's circumcenter should not be within the offset volume
//     CGAL_assertion_code(const Point_3& ch_cc = circumcenter(ch);)
//     CGAL_assertion_code(const Ball_3 ch_cc_offset_ball = ball(ch_cc, m_sq_offset);)
//     CGAL_assertion(!m_oracle.do_intersect(ch_cc_offset_ball));
//
//     if(is_neighbor_cc_in_offset)
//     {
//       const Point_3 closest_pt = m_oracle.closest_point(neighbor_cc);
//       CGAL_assertion(closest_pt != neighbor_cc);
//       const Point_3& ch_cc = circumcenter(ch);
//       Point_3 pt1, pt2, pt3;
//       int fi = ch->index(neighbor);
//       int k = 0;
//
//       for(int j = 0; j < 4; ++j)
//       {
//         if(j != fi)
//         {
//           const Point_3& p = ch->vertex(j)->point();
//           if(k == 0) pt1 = p;
//           else if(k == 1) pt2 = p;
//           else pt3 = p;
//           ++k;
//         }
//       }
//
//       const Point_3 closest_pt1 = m_oracle.closest_point(pt1);
//       const Point_3 closest_pt2 = m_oracle.closest_point(pt2);
//       const Point_3 closest_pt3 = m_oracle.closest_point(pt3);
//       Vector_3 v1 = vector(closest_pt1, pt1);
//       Vector_3 v2 = vector(closest_pt2, pt2);
//       Vector_3 v3 = vector(closest_pt3, pt3);
//       v1 = scale(v1, FT(1) / approximate_sqrt(sq_length(v1)));
//       v2 = scale(v2, FT(1) / approximate_sqrt(sq_length(v2)));
//       v3 = scale(v3, FT(1) / approximate_sqrt(sq_length(v3)));
//       std::cout << "=======================================" << std::endl;
//       std::cout << "v1 = " << v1 << std::endl;
//       std::cout << "v2 = " << v2 << std::endl;
//       std::cout << "v3 = " << v3 << std::endl;
//       Vector_3 offset_dir;
//
//       FT d12 = dot(v1, v2);
//       FT d13 = dot(v1, v3);
//       FT d23 = dot(v2, v3);
//
//       // helper lambdas
//       const FT same_threshold = FT(0.98);
//       auto same_dir = [&](const Vector_3& a, const Vector_3& b) -> bool {
//         return dot(a, b) > same_threshold;
//       };
//
//       auto normalize = [&](Vector_3 u) -> Vector_3 {
//         const FT len2 = sq_length(u);
//         CGAL_assertion(len2 > FT(0));
//         return scale(u, FT(1) / approximate_sqrt(len2));
//       };
//
//       // classify
//       const bool v12_same = same_dir(v1, v2);
//       const bool v13_same = same_dir(v1, v3);
//       const bool v23_same = same_dir(v2, v3);
//       Vector_3 vc, vd;
//
//       if (v12_same && v13_same && v23_same)
//       {
// #ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
//         std::cout << "all vectors are the same😯" << std::endl;
// #endif
//         // Move from closest_pt in direction v1 by offset length
//         steiner_point = translate(closest_pt, scale(v1, m_offset));
//         return true;
//       }
//       else if(!v12_same && !v13_same && !v23_same)
//       {
//         if(m_oracle.first_intersection(ch_cc, neighbor_cc, steiner_point, m_offset))
//         {
// #ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
//           std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
// #endif
//           return true;
//         }
//       }
//         // "two vectors are the same" branch (your requested behavior)
//         // Find the repeated direction vc and the different one vd
//
//         else if (v12_same && !v13_same)      { vc = v1; vd = v3; }
//         else if (v13_same && !v12_same) { vc = v1; vd = v2; }
//         else if (v23_same && !v12_same) { vc = v2; vd = v1; }
//         else
//         {
//           // If it's "all different" (or ambiguous), fall back to average direction of all three.
//           // (You can replace this with your own desired behavior.)
// #ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
//           std::cout << "all vectors are different/ambiguous -> fallback average" << std::endl;
// #endif
//           Vector_3 avg3 = normalize(v1 + v2 + v3);
//           steiner_point = translate(closest_pt, scale(avg3, m_offset));
//           return true;
//         }
//       //       if(m_oracle.first_intersection(ch_cc, closest_pt, steiner_point, m_offset))
//       //       {
//       // #ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
//       //         std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
//       // #endif
//       //         return true;
//       //       }
//       std::cout << "offset_dir = " << offset_dir << std::endl;
//       offset_dir = scale(offset_dir, FT(1) / approximate_sqrt(sq_length(offset_dir)));
//       std::cout << "unit length offset_dir = " << offset_dir << std::endl;
//       const Point_3 end_pt = translate(closest_pt, scale(offset_dir, 2*m_offset));
//       if(m_oracle.first_intersection(closest_pt, end_pt, steiner_point, m_offset)) {
//         std::cout << "omg steiner point constructed inside concave edge!" << std::endl;
//         return true;
//       }
//       // const Point_3& ch_cc = circumcenter(ch);
//
//       // If the voronoi edge intersects the offset, the steiner point is the first intersection
//       if(m_oracle.first_intersection(ch_cc, neighbor_cc, steiner_point, m_offset))
//       {
// #ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
//         std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
// #endif
//         return true;
//       }
//     }
//
//     Tetrahedron_with_outside_info<Geom_traits> tet(neighbor, geom_traits());
//     if(m_oracle.do_intersect(tet))
//     {
//       // steiner point is the closest point on input from cell centroid with offset
//       const Point_3 closest_pt = m_oracle.closest_point(neighbor_cc);
//       CGAL_assertion(closest_pt != neighbor_cc);
//       // const Ball_3 closest_pt_offset_ball = ball(neighbor_cc, m_sq_offset*1.1);
//
//       Vector_3 unit = vector(closest_pt, neighbor_cc);
//
//       // PMP::internal::normalize() requires sqrt()
//       unit = scale(unit, m_offset / approximate_sqrt(geom_traits().compute_squared_length_3_object()(unit)));
//       steiner_point = translate(closest_pt, unit);
//
// #ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
//       std::cout << "Steiner found through neighboring tet intersecting the input: " << steiner_point << std::endl;
//       std::cout << "Closest point: " << closest_pt << std::endl;
//       std::cout << "Direction: " << vector(closest_pt, neighbor_cc) << std::endl;
// #endif
//
//       return true;
//     }
//
// #ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
//     std::cout << "No Steiner point" << std::endl;
// #endif
//
//     return false;
//   }

    bool compute_steiner_point3(const Cell_handle ch,
                             const Cell_handle neighbor,
                             Point_3& steiner_point) const
  {
    CGAL_precondition(!m_tr.is_infinite(neighbor));

    typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
    typename Geom_traits::Construct_vector_3 vector = geom_traits().construct_vector_3_object();
    typename Geom_traits::Construct_translated_point_3 translate = geom_traits().construct_translated_point_3_object();
    typename Geom_traits::Construct_scaled_vector_3 scale = geom_traits().construct_scaled_vector_3_object();
    typename Geom_traits::Compute_squared_length_3 sq_length = geom_traits().compute_squared_length_3_object();
    typename Geom_traits::Compute_scalar_product_3 dot = geom_traits().compute_scalar_product_3_object();

    const Point_3& neighbor_cc = circumcenter(neighbor);
    const Ball_3 neighbor_cc_offset_ball = ball(neighbor_cc, m_sq_offset);
    const bool is_neighbor_cc_in_offset = m_oracle.do_intersect(neighbor_cc_offset_ball);

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "Compute_steiner_point(" << &*ch << ", " << &*neighbor << ")" << std::endl;

    const Point_3& chc = circumcenter(ch);
    std::cout << "CH" << std::endl;
    std::cout << "\t" << ch->vertex(0)->point() << std::endl;
    std::cout << "\t" << ch->vertex(1)->point() << std::endl;
    std::cout << "\t" << ch->vertex(2)->point() << std::endl;
    std::cout << "\t" << ch->vertex(3)->point() << std::endl;
    std::cout << "cc: " << chc << std::endl;
    std::cout << "CC Distance to input: " << CGAL::squared_distance(chc, m_oracle.closest_point(chc)) << std::endl;

    std::cout << "NCH" << std::endl;
    std::cout << "\t" << neighbor->vertex(0)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(1)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(2)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(3)->point() << std::endl;
    std::cout << "ncc: " << neighbor_cc << std::endl;
    std::cout << "NC Distance to input: " << CGAL::squared_distance(neighbor_cc, m_oracle.closest_point(neighbor_cc)) << std::endl;
    std::cout << "is_neighbor_cc_in_offset = " << is_neighbor_cc_in_offset << std::endl;

    std::cout << "squared offset " << m_sq_offset << std::endl;
#endif

    // ch's circumcenter should not be within the offset volume
    CGAL_assertion_code(const Point_3& ch_cc = circumcenter(ch);)
    CGAL_assertion_code(const Ball_3 ch_cc_offset_ball = ball(ch_cc, m_sq_offset);)
    CGAL_assertion(!m_oracle.do_intersect(ch_cc_offset_ball));

    if(is_neighbor_cc_in_offset)
    {
      const Point_3& ch_cc = circumcenter(ch);
      Point_3 pt1, pt2, pt3;
      int fi = ch->index(neighbor);
      int k = 0;

      for(int j = 0; j < 4; ++j)
      {
        if(j != fi)
        {
          const Point_3& p = ch->vertex(j)->point();
          if(k == 0) pt1 = p;
          else if(k == 1) pt2 = p;
          else pt3 = p;
          ++k;
        }
      }
      Facet f = Facet(ch, fi);
      bool traversable = less_squared_radius_of_min_empty_sphere(m_sq_alpha, f, m_tr);
      std::cout << "Checking traversable: " << traversable << std::endl; // Debug line before the condition

      if (!traversable) {
        std::cout << "🤮🤮🤮!!!!!!face is not traversible!!!!!!🤮🤮🤮" << std::endl;


        const Point_3 closest_pt1 = m_oracle.closest_point(pt1);
        const Point_3 closest_pt2 = m_oracle.closest_point(pt2);
        const Point_3 closest_pt3 = m_oracle.closest_point(pt3);
        Vector_3 v1 = vector(closest_pt1, pt1);
        Vector_3 v2 = vector(closest_pt2, pt2);
        Vector_3 v3 = vector(closest_pt3, pt3);
        v1 = scale(v1, FT(1) / approximate_sqrt(sq_length(v1)));
        v2 = scale(v2, FT(1) / approximate_sqrt(sq_length(v2)));
        v3 = scale(v3, FT(1) / approximate_sqrt(sq_length(v3)));
        // std::cout << "=======================================" << std::endl;
        // std::cout << "v1 = " << v1 << std::endl;
        // std::cout << "v2 = " << v2 << std::endl;
        // std::cout << "v3 = " << v3 << std::endl;
        Vector_3 offset_dir;

        FT d12 = dot(v1, v2);
        FT d13 = dot(v1, v3);
        FT d23 = dot(v2, v3);

        // helper lambdas
        const FT same_threshold = FT(0.95);
        auto same_dir = [&](const Vector_3& a, const Vector_3& b) -> bool {
          return dot(a, b) > same_threshold;
        };

        auto normalize = [&](Vector_3 u) -> Vector_3 {
          const FT len2 = sq_length(u);
          CGAL_assertion(len2 > FT(0));
          return scale(u, FT(1) / approximate_sqrt(len2));
        };

        // classify
        const bool v12_same = same_dir(v1, v2);
        const bool v13_same = same_dir(v1, v3);
        const bool v23_same = same_dir(v2, v3);
        Vector_3 vc(0, 0, 0);  // Initialize to zero vector
        Vector_3 vd(0, 0, 0);  // Initialize to zero vector
        Point_3 pc,pd;
        if (v12_same && !v13_same) {
          vc = v1;
          pc = closest_pt1;
          vd = v3;
          pd = closest_pt3;
        }
        if (v13_same && !v23_same) {
          vc = v1;
          pc = closest_pt1;
          vd = v2;
          pd = closest_pt2;
        }
        if (v23_same && !v12_same) {
          vc = v2;
          pc = closest_pt2;
          vd = v1;
          pd = closest_pt1;
        }
        const FT epsilon = FT(1e-10);  // Small threshold to account for floating-point precision

        // Function to check if a vector is non-zero
        auto is_non_zero = [&](const Vector_3& v) -> bool {
          return sq_length(v) > epsilon;  // Check if the squared length is greater than epsilon
        };
        if (is_non_zero(vc) && is_non_zero(vd)) {
          // std::cout << "✈️plane✈️ 1 = " << vc << std::endl;
          // std::cout << "✈️plane✈️ 2 = " << vd << std::endl;
          Vector_3 intersection_direction(
            vc.y() * vd.z() - vc.z() * vd.y(),  // x component
            vc.z() * vd.x() - vc.x() * vd.z(),  // y component
            vc.x() * vd.y() - vc.y() * vd.x()   // z component
          );
          // planes should not be parrallel
          if (sq_length(intersection_direction) > FT(1e-10)) {
            // std::cout << "intersection direction = " << intersection_direction << std::endl;
            Vector_3 vec_to_neighbor = vector(pc, ch_cc);
            FT t = dot(vec_to_neighbor, intersection_direction) / dot(intersection_direction, intersection_direction);
            Point_3 closest_point_on_line = pc + t * intersection_direction;
            // std::cout << "closest_point on_line = " << closest_point_on_line << ", with neighbor_cc = " << neighbor_cc << std::endl;
            if(m_oracle.first_intersection(closest_point_on_line, neighbor_cc, steiner_point, m_offset))
            {
              // std::cout << "steiner point = " << steiner_point << std::endl;
              // std::cout << "🤩steiner point found from closesest point on line --> neighbor_cc🤩" << std::endl;
              return true;
            }
            if(m_oracle.first_intersection(ch_cc, closest_point_on_line, steiner_point, m_offset))
            {
              // std::cout << "steiner point = " << steiner_point << std::endl;
              // std::cout << "🥰steiner point found from ch_cc --> closesest point on line🥰" << std::endl;
              return true;
            }
          }
        }
      }

      // If the voronoi edge intersects the offset, the steiner point is the first intersection
      if(m_oracle.first_intersection(ch_cc, neighbor_cc, steiner_point, m_offset))
      {
        // std::cout << "steiner point = " << steiner_point << std::endl;
        // std::cout << "🧠steiner point made with first intersection with offset of voronoi edge🧠" << std::endl;
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
        std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
#endif
        return true;
      }
    }

    Tetrahedron_with_outside_info<Geom_traits> tet(neighbor, geom_traits());
    if(m_oracle.do_intersect(tet))
    {
      // steiner point is the closest point on input from cell centroid with offset
      const Point_3 closest_pt = m_oracle.closest_point(neighbor_cc);
      CGAL_assertion(closest_pt != neighbor_cc);

      Vector_3 unit = vector(closest_pt, neighbor_cc);

      // PMP::internal::normalize() requires sqrt()
      unit = scale(unit, m_offset / approximate_sqrt(geom_traits().compute_squared_length_3_object()(unit)));
      steiner_point = translate(closest_pt, unit);
      // std::cout << "=======================================" << std::endl;
      // std::cout << "👀steiner point based on closest point to neighbor_cc👀" << std::endl;

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
      std::cout << "Steiner found through neighboring tet intersecting the input: " << steiner_point << std::endl;
      std::cout << "Closest point: " << closest_pt << std::endl;
      std::cout << "Direction: " << vector(closest_pt, neighbor_cc) << std::endl;
#endif

      return true;
    }

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "No Steiner point" << std::endl;
#endif

    return false;
  }

bool compute_steiner_point4(const Cell_handle ch,
                           const Cell_handle neighbor,
                           Point_3& steiner_point) const
{
  CGAL_precondition(!m_tr.is_infinite(neighbor));
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Vector_3 Vector_3;
  typedef typename Geom_traits::Plane_3 Plane_3;
  typedef typename Geom_traits::Line_3 Line_3;

  // --- 1. Traits and Setup ---
  typename Geom_traits::Construct_vector_3 vector = geom_traits().construct_vector_3_object();
  typename Geom_traits::Construct_translated_point_3 translate = geom_traits().construct_translated_point_3_object();
  typename Geom_traits::Construct_scaled_vector_3 scale = geom_traits().construct_scaled_vector_3_object();
  typename Geom_traits::Compute_squared_length_3 sq_length = geom_traits().compute_squared_length_3_object();
  typename Geom_traits::Compute_scalar_product_3 dot = geom_traits().compute_scalar_product_3_object();
  typename Geom_traits::Construct_plane_3 construct_plane = geom_traits().construct_plane_3_object();
  typename Geom_traits::Construct_projected_point_3 project_pt = geom_traits().construct_projected_point_3_object();

  const Point_3& ch_cc = circumcenter(ch);
  const Point_3& neighbor_cc = circumcenter(neighbor);

  Point_3 pts[3];
  int fi = ch->index(neighbor);
  for(int j=0, k=0; j<4; ++j) {
    if(j != fi) pts[k++] = ch->vertex(j)->point();
  }

  // --- 2. Check Traversability ---
  Facet f = Facet(ch, fi);
  if (less_squared_radius_of_min_empty_sphere(m_sq_alpha, f, m_tr)) {
    return compute_steiner_point2(ch, neighbor, steiner_point);
  }

  // --- 3. Concave Edge Logic (Ray-casting toward intersection line) ---
  const Point_3 cpts[3] = { m_oracle.closest_point(pts[0]),
                            m_oracle.closest_point(pts[1]),
                            m_oracle.closest_point(pts[2]) };

  Vector_3 vs[3] = { vector(cpts[0], pts[0]),
                     vector(cpts[1], pts[1]),
                     vector(cpts[2], pts[2]) };

  auto normalize = [&](Vector_3 v) {
    FT len2 = sq_length(v);
    return (len2 > 1e-16) ? scale(v, FT(1) / approximate_sqrt(len2)) : Vector_3(0,0,0);
  };
  Vector_3 n1 = normalize(vs[0]), n2 = normalize(vs[1]), n3 = normalize(vs[2]);

  // Select the two normals that are most different (closest to orthogonal)
  FT d12 = CGAL::abs(dot(n1, n2));
  FT d13 = CGAL::abs(dot(n1, n3));
  FT d23 = CGAL::abs(dot(n2, n3));

  Plane_3 planeA, planeB;
  bool planes_found = false;
  const FT diff_limit = FT(0.999); // Ignore if normals are within ~2 degrees

  if (d12 <= d13 && d12 <= d23 && d12 < diff_limit) {
    planeA = construct_plane(pts[0], n1);
    planeB = construct_plane(pts[1], n2);
    planes_found = true;
  } else if (d13 <= d12 && d13 <= d23 && d13 < diff_limit) {
    planeA = construct_plane(pts[0], n1);
    planeB = construct_plane(pts[2], n3);
    planes_found = true;
  } else if (d23 < diff_limit) {
    planeA = construct_plane(pts[1], n2);
    planeB = construct_plane(pts[2], n3);
    planes_found = true;
  }

  if (planes_found) {
    auto result = CGAL::intersection(planeA, planeB);
    Line_3 intersect_line;
    // Using CGAL::assign to avoid boost/std versioning issues
    if (result && CGAL::assign(intersect_line, *result)) {

      Point_3 intersection_point = project_pt(intersect_line, ch_cc);

      // PERFORMANCE GUARD: Only ray-cast if the line is relatively close to the facet.
      // If the line is extremely far, the planes are too parallel to be useful.
      FT facet_size_sq = CGAL::squared_distance(pts[0], pts[1]);
      if (CGAL::squared_distance(intersection_point, pts[0]) < 100.0 * facet_size_sq) {

        // SEARCH 1: From ch_cc toward the intersection line
        if (m_oracle.first_intersection(ch_cc, intersection_point, steiner_point, m_offset)) {
          return true;
        }

        // SEARCH 2: From intersection line toward neighbor_cc
        if (m_oracle.first_intersection(intersection_point, neighbor_cc, steiner_point, m_offset)) {
          return true;
        }
      }
    }
  }

  // --- 4. Fallback ---
  // If no intersection was found on the segments, use the standard logic.
  return compute_steiner_point2(ch, neighbor, steiner_point);
}

  // normal steiner computation
  bool compute_steiner_point(const Cell_handle ch,
                             const Cell_handle neighbor,
                             Point_3& steiner_point) const
  {
    CGAL_precondition(!m_tr.is_infinite(neighbor));

    typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
    typename Geom_traits::Construct_vector_3 vector = geom_traits().construct_vector_3_object();
    typename Geom_traits::Construct_translated_point_3 translate = geom_traits().construct_translated_point_3_object();
    typename Geom_traits::Construct_scaled_vector_3 scale = geom_traits().construct_scaled_vector_3_object();

    const Point_3& neighbor_cc = circumcenter(neighbor);
    const Ball_3 neighbor_cc_offset_ball = ball(neighbor_cc, m_sq_offset);
    const bool is_neighbor_cc_in_offset = m_oracle.do_intersect(neighbor_cc_offset_ball);

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "Compute_steiner_point(" << &*ch << ", " << &*neighbor << ")" << std::endl;

    const Point_3& chc = circumcenter(ch);
    std::cout << "CH" << std::endl;
    std::cout << "\t" << ch->vertex(0)->point() << std::endl;
    std::cout << "\t" << ch->vertex(1)->point() << std::endl;
    std::cout << "\t" << ch->vertex(2)->point() << std::endl;
    std::cout << "\t" << ch->vertex(3)->point() << std::endl;
    std::cout << "cc: " << chc << std::endl;
    std::cout << "CC Distance to input: " << CGAL::squared_distance(chc, m_oracle.closest_point(chc)) << std::endl;

    std::cout << "NCH" << std::endl;
    std::cout << "\t" << neighbor->vertex(0)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(1)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(2)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(3)->point() << std::endl;
    std::cout << "ncc: " << neighbor_cc << std::endl;
    std::cout << "NC Distance to input: " << CGAL::squared_distance(neighbor_cc, m_oracle.closest_point(neighbor_cc)) << std::endl;
    std::cout << "is_neighbor_cc_in_offset = " << is_neighbor_cc_in_offset << std::endl;

    std::cout << "squared offset " << m_sq_offset << std::endl;
#endif

    // ch's circumcenter should not be within the offset volume
    CGAL_assertion_code(const Point_3& ch_cc = circumcenter(ch);)
    CGAL_assertion_code(const Ball_3 ch_cc_offset_ball = ball(ch_cc, m_sq_offset);)
    CGAL_assertion(!m_oracle.do_intersect(ch_cc_offset_ball));

    if(is_neighbor_cc_in_offset)
    {
      const Point_3& ch_cc = circumcenter(ch);

      // If the voronoi edge intersects the offset, the steiner point is the first intersection
      if(m_oracle.first_intersection(ch_cc, neighbor_cc, steiner_point, m_offset))
      {
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
        std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
#endif
        return true;
      }
    }

    Tetrahedron_with_outside_info<Geom_traits> tet(neighbor, geom_traits());
    if(m_oracle.do_intersect(tet))
    {
      // steiner point is the closest point on input from cell centroid with offset
      const Point_3 closest_pt = m_oracle.closest_point(neighbor_cc);
      CGAL_assertion(closest_pt != neighbor_cc);

      Vector_3 unit = vector(closest_pt, neighbor_cc);

      // PMP::internal::normalize() requires sqrt()
      unit = scale(unit, m_offset / approximate_sqrt(geom_traits().compute_squared_length_3_object()(unit)));
      steiner_point = translate(closest_pt, unit);

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
      std::cout << "Steiner found through neighboring tet intersecting the input: " << steiner_point << std::endl;
      std::cout << "Closest point: " << closest_pt << std::endl;
      std::cout << "Direction: " << vector(closest_pt, neighbor_cc) << std::endl;
#endif

      return true;
    }

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "No Steiner point" << std::endl;
#endif

    return false;
  }

  bool compute_edge_optimum(const Point_3& a,
                            const Point_3& b,
                            const Point_3& edge_mid,
                            Point_3& p_opt) const
{
  typename Geom_traits::Construct_vector_3 vector =
      geom_traits().construct_vector_3_object();
  auto cross = geom_traits().construct_cross_product_vector_3_object();
  auto sq_length = geom_traits().compute_squared_length_3_object();

  const Point_3 closest_a = m_oracle.closest_point(a);
  const Point_3 closest_b = m_oracle.closest_point(b);

  const Vector_3 e  = vector(a, b);
  const Vector_3 na = vector(a, closest_a);
  const Vector_3 nb = vector(b, closest_b);

  const double len_e  = std::sqrt(CGAL::to_double(sq_length(e)));
  const double len_na = std::sqrt(CGAL::to_double(sq_length(na)));
  const double len_nb = std::sqrt(CGAL::to_double(sq_length(nb)));

  const double eps = 1e-12;
  if (len_e < eps || len_na < eps || len_nb < eps) {
    p_opt = edge_mid;
    return false;
  }

  double sin_a =
      std::sqrt(CGAL::to_double(sq_length(cross(na, e)))) / (len_na * len_e);
  double sin_b =
      std::sqrt(CGAL::to_double(sq_length(cross(nb, e)))) / (len_nb * len_e);

  sin_a = std::max(0.0, std::min(1.0, sin_a));
  sin_b = std::max(0.0, std::min(1.0, sin_b));

  const double denom = sin_a + sin_b;
  if (denom < eps) {
    p_opt = edge_mid;
    return false;
  }

  const double t = sin_b / denom;
  p_opt = a + t * e;
  return true;
}

    // modified steiner computation
  bool compute_steiner_point8(const Cell_handle ch,
                             const Cell_handle neighbor,
                             Point_3& steiner_point) const
  {
    CGAL_precondition(!m_tr.is_infinite(neighbor));

    typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
    typename Geom_traits::Construct_vector_3 vector = geom_traits().construct_vector_3_object();
    typename Geom_traits::Construct_translated_point_3 translate = geom_traits().construct_translated_point_3_object();
    typename Geom_traits::Construct_scaled_vector_3 scale = geom_traits().construct_scaled_vector_3_object();
    typename Geom_traits::Construct_midpoint_3 construct_midpoint = geom_traits().construct_midpoint_3_object();
    typename Geom_traits::Compute_squared_distance_3 sq_dist = geom_traits().compute_squared_distance_3_object();
    auto cross = geom_traits().construct_cross_product_vector_3_object();
    auto sq_length = geom_traits().compute_squared_length_3_object();

    const Point_3& neighbor_cc = circumcenter(neighbor);
    const Ball_3 neighbor_cc_offset_ball = ball(neighbor_cc, m_sq_offset);
    const bool is_neighbor_cc_in_offset = m_oracle.do_intersect(neighbor_cc_offset_ball);


#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "Compute_steiner_point(" << &*ch << ", " << &*neighbor << ")" << std::endl;

    const Point_3& chc = circumcenter(ch);
    std::cout << "CH" << std::endl;
    std::cout << "\t" << ch->vertex(0)->point() << std::endl;
    std::cout << "\t" << ch->vertex(1)->point() << std::endl;
    std::cout << "\t" << ch->vertex(2)->point() << std::endl;
    std::cout << "\t" << ch->vertex(3)->point() << std::endl;
    std::cout << "cc: " << chc << std::endl;
    std::cout << "CC Distance to input: " << CGAL::squared_distance(chc, m_oracle.closest_point(chc)) << std::endl;

    std::cout << "NCH" << std::endl;
    std::cout << "\t" << neighbor->vertex(0)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(1)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(2)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(3)->point() << std::endl;
    std::cout << "ncc: " << neighbor_cc << std::endl;
    std::cout << "NC Distance to input: " << CGAL::squared_distance(neighbor_cc, m_oracle.closest_point(neighbor_cc)) << std::endl;
    std::cout << "is_neighbor_cc_in_offset = " << is_neighbor_cc_in_offset << std::endl;

    std::cout << "squared offset " << m_sq_offset << std::endl;
#endif

    // ch's circumcenter should not be within the offset volume
    CGAL_assertion_code(const Point_3& ch_cc = circumcenter(ch);)
    CGAL_assertion_code(const Ball_3 ch_cc_offset_ball = ball(ch_cc, m_sq_offset);)
    CGAL_assertion(!m_oracle.do_intersect(ch_cc_offset_ball));

  Point_3 pts[3];
  int fi = ch->index(neighbor);
  for(int j=0, k=0; j<4; ++j) {
    if(j != fi) pts[k++] = ch->vertex(j)->point();
  }

  // --- 2. Check Traversability ---
  Facet f = Facet(ch, fi);


  const bool traversible = less_squared_radius_of_min_empty_sphere(m_sq_alpha, f, m_tr);
  if (traversible) {
    return compute_steiner_point_normal(ch, neighbor, steiner_point);

  }


  if (!traversible) {
    Tetrahedron_with_outside_info<Geom_traits> tet(neighbor, geom_traits());
    if(m_oracle.do_intersect(tet))
  {
  // steiner point is the closest point on input from cell centroid with offset
  const Point_3 closest_pt = m_oracle.closest_point(neighbor_cc);
  CGAL_assertion(closest_pt != neighbor_cc);
      Point_3 face_mid(
        (pts[0].x() + pts[1].x() + pts[2].x()) / 3.0,
        (pts[0].y() + pts[1].y() + pts[2].y()) / 3.0,
        (pts[0].z() + pts[1].z() + pts[2].z()) / 3.0);

  // compute midpoints on triangle edges
  Point_3 edge_mid01 = construct_midpoint(pts[0], pts[1]);
  Point_3 edge_mid12 = construct_midpoint(pts[1], pts[2]);
  Point_3 edge_mid20 = construct_midpoint(pts[2], pts[0]);

  // // compute optimal points on each edge
  // Point_3 p_opt01 = edge_mid01;
  // Point_3 p_opt12 = edge_mid12;
  // Point_3 p_opt20 = edge_mid20;
  //
  // compute_edge_optimum(pts[0], pts[1], edge_mid01, p_opt01);
  // compute_edge_optimum(pts[1], pts[2], edge_mid12, p_opt12);
  // compute_edge_optimum(pts[2], pts[0], edge_mid20, p_opt20);

  // closest points on input for the 3 edge midpoints
  const Point_3 closest_pt01 = m_oracle.closest_point(edge_mid01);
  const Point_3 closest_pt12 = m_oracle.closest_point(edge_mid12);
  const Point_3 closest_pt20 = m_oracle.closest_point(edge_mid20);

  // // closest points on input for the 3 optimal edge points
  // const Point_3 closest_pt_opt01 = m_oracle.closest_point(p_opt01);
  // const Point_3 closest_pt_opt12 = m_oracle.closest_point(p_opt12);
  // const Point_3 closest_pt_opt20 = m_oracle.closest_point(p_opt20);

  // squared distances to input for all 6 candidates
  FT d_mid01_sq = sq_dist(edge_mid01, closest_pt01);
  FT d_mid12_sq = sq_dist(edge_mid12, closest_pt12);
  FT d_mid20_sq = sq_dist(edge_mid20, closest_pt20);

  // FT d_opt01_sq = sq_dist(p_opt01, closest_pt_opt01);
  // FT d_opt12_sq = sq_dist(p_opt12, closest_pt_opt12);
  // FT d_opt20_sq = sq_dist(p_opt20, closest_pt_opt20);

  // choose the point farthest from the input
  Point_3 edge_mid = edge_mid01;
  FT max_d_sq = d_mid01_sq;

  if (d_mid12_sq > max_d_sq) {
    max_d_sq = d_mid12_sq;
    edge_mid = edge_mid12;
  }
  if (d_mid20_sq > max_d_sq) {
    max_d_sq = d_mid20_sq;
    edge_mid = edge_mid20;
  }

  // if (d_opt01_sq > max_d_sq) {
  //   max_d_sq = d_opt01_sq;
  //   face_mid = p_opt01;
  // }
  // if (d_opt12_sq > max_d_sq) {
  //   max_d_sq = d_opt12_sq;
  //   face_mid = p_opt12;
  // }
  // if (d_opt20_sq > max_d_sq) {
  //   max_d_sq = d_opt20_sq;
  //   face_mid = p_opt20;
  // }
      Vector_3 from_face_mid = vector(face_mid,edge_mid);
      // std::cout << "=========================================================" << std::endl;
      // std::cout << "opt point should lay between: face_mid = " << face_mid << " and edge_mid = " << edge_mid << std::endl;
      face_mid = face_mid + 0.95*from_face_mid;
      // std::cout << "opt point = " << face_mid << std::endl;

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
  std::cout << "edge_mid01 = " << edge_mid01 << ", d^2 = " << d_mid01_sq << std::endl;
  std::cout << "edge_mid12 = " << edge_mid12 << ", d^2 = " << d_mid12_sq << std::endl;
  std::cout << "edge_mid20 = " << edge_mid20 << ", d^2 = " << d_mid20_sq << std::endl;

  std::cout << "p_opt01    = " << p_opt01    << ", d^2 = " << d_opt01_sq << std::endl;
  std::cout << "p_opt12    = " << p_opt12    << ", d^2 = " << d_opt12_sq << std::endl;
  std::cout << "p_opt20    = " << p_opt20    << ", d^2 = " << d_opt20_sq << std::endl;

  std::cout << "selected face_mid = " << face_mid << ", max d^2 = " << max_d_sq << std::endl;
#endif

  if(m_oracle.first_intersection(face_mid, closest_pt, steiner_point, m_offset))
  {
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
#endif
    return true;
  }
  }
  }
  return compute_steiner_point_normal(ch, neighbor, steiner_point);
  }

    // normal steiner computation
  bool compute_steiner_point_testtest(const Cell_handle ch,
                             const Cell_handle neighbor,
                             Point_3& steiner_point) const
  {
    CGAL_precondition(!m_tr.is_infinite(neighbor));

    typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
    typename Geom_traits::Construct_vector_3 vector = geom_traits().construct_vector_3_object();
    typename Geom_traits::Construct_translated_point_3 translate = geom_traits().construct_translated_point_3_object();
    typename Geom_traits::Construct_scaled_vector_3 scale = geom_traits().construct_scaled_vector_3_object();
    typename Geom_traits::Construct_midpoint_3 construct_midpoint = geom_traits().construct_midpoint_3_object();
    typename Geom_traits::Compute_squared_distance_3 sq_dist = geom_traits().compute_squared_distance_3_object();
    auto cross = geom_traits().construct_cross_product_vector_3_object();
    auto sq_length = geom_traits().compute_squared_length_3_object();

    const Point_3& neighbor_cc = circumcenter(neighbor);
    const Ball_3 neighbor_cc_offset_ball = ball(neighbor_cc, m_sq_offset);
    const bool is_neighbor_cc_in_offset = m_oracle.do_intersect(neighbor_cc_offset_ball);


#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "Compute_steiner_point(" << &*ch << ", " << &*neighbor << ")" << std::endl;

    const Point_3& chc = circumcenter(ch);
    std::cout << "CH" << std::endl;
    std::cout << "\t" << ch->vertex(0)->point() << std::endl;
    std::cout << "\t" << ch->vertex(1)->point() << std::endl;
    std::cout << "\t" << ch->vertex(2)->point() << std::endl;
    std::cout << "\t" << ch->vertex(3)->point() << std::endl;
    std::cout << "cc: " << chc << std::endl;
    std::cout << "CC Distance to input: " << CGAL::squared_distance(chc, m_oracle.closest_point(chc)) << std::endl;

    std::cout << "NCH" << std::endl;
    std::cout << "\t" << neighbor->vertex(0)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(1)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(2)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(3)->point() << std::endl;
    std::cout << "ncc: " << neighbor_cc << std::endl;
    std::cout << "NC Distance to input: " << CGAL::squared_distance(neighbor_cc, m_oracle.closest_point(neighbor_cc)) << std::endl;
    std::cout << "is_neighbor_cc_in_offset = " << is_neighbor_cc_in_offset << std::endl;

    std::cout << "squared offset " << m_sq_offset << std::endl;
#endif

    // ch's circumcenter should not be within the offset volume
    CGAL_assertion_code(const Point_3& ch_cc = circumcenter(ch);)
    CGAL_assertion_code(const Ball_3 ch_cc_offset_ball = ball(ch_cc, m_sq_offset);)
    CGAL_assertion(!m_oracle.do_intersect(ch_cc_offset_ball));

  Point_3 pts[3];
  int fi = ch->index(neighbor);
  for(int j=0, k=0; j<4; ++j) {
    if(j != fi) pts[k++] = ch->vertex(j)->point();
  }

  // --- 2. Check Traversability ---
  Facet f = Facet(ch, fi);
  const bool traversible = less_squared_radius_of_min_empty_sphere(m_sq_alpha, f, m_tr);
  if (traversible) {
    return compute_steiner_point_normal(ch, neighbor, steiner_point);
    // std::cout << "face is traversible!" << std::endl;
  }
  if (!traversible) {
    // std::cout << "face is NOT traversible!" << std::endl;
    Tetrahedron_with_outside_info<Geom_traits> tet(neighbor, geom_traits());
    if(m_oracle.do_intersect(tet))
{
  // steiner point is the closest point on input from cell centroid with offset
  const Point_3 closest_pt = m_oracle.closest_point(neighbor_cc);
  CGAL_assertion(closest_pt != neighbor_cc);

  // compute midpoints on triangle edges
  Point_3 edge_mid01 = construct_midpoint(pts[0], pts[1]);
  Point_3 edge_mid12 = construct_midpoint(pts[1], pts[2]);
  Point_3 edge_mid20 = construct_midpoint(pts[2], pts[0]);

  // compute optimal points on each edge
  Point_3 p_opt01 = edge_mid01;
  Point_3 p_opt12 = edge_mid12;
  Point_3 p_opt20 = edge_mid20;

  compute_edge_optimum(pts[0], pts[1], edge_mid01, p_opt01);
  compute_edge_optimum(pts[1], pts[2], edge_mid12, p_opt12);
  compute_edge_optimum(pts[2], pts[0], edge_mid20, p_opt20);

  // closest points on input for the 3 edge midpoints
  const Point_3 closest_pt01 = m_oracle.closest_point(edge_mid01);
  const Point_3 closest_pt12 = m_oracle.closest_point(edge_mid12);
  const Point_3 closest_pt20 = m_oracle.closest_point(edge_mid20);

  // closest points on input for the 3 optimal edge points
  const Point_3 closest_pt_opt01 = m_oracle.closest_point(p_opt01);
  const Point_3 closest_pt_opt12 = m_oracle.closest_point(p_opt12);
  const Point_3 closest_pt_opt20 = m_oracle.closest_point(p_opt20);

  // squared distances to input for all 6 candidates
  FT d_mid01_sq = sq_dist(edge_mid01, closest_pt01);
  FT d_mid12_sq = sq_dist(edge_mid12, closest_pt12);
  FT d_mid20_sq = sq_dist(edge_mid20, closest_pt20);

  FT d_opt01_sq = sq_dist(p_opt01, closest_pt_opt01);
  FT d_opt12_sq = sq_dist(p_opt12, closest_pt_opt12);
  FT d_opt20_sq = sq_dist(p_opt20, closest_pt_opt20);

  // choose the point farthest from the input
  Point_3 face_mid = edge_mid01;
  FT max_d_sq = d_mid01_sq;

  if (d_mid12_sq > max_d_sq) {
    max_d_sq = d_mid12_sq;
    face_mid = edge_mid12;
  }
  if (d_mid20_sq > max_d_sq) {
    max_d_sq = d_mid20_sq;
    face_mid = edge_mid20;
  }
  if (d_opt01_sq > max_d_sq) {
    max_d_sq = d_opt01_sq;
    face_mid = p_opt01;
  }
  if (d_opt12_sq > max_d_sq) {
    max_d_sq = d_opt12_sq;
    face_mid = p_opt12;
  }
  if (d_opt20_sq > max_d_sq) {
    max_d_sq = d_opt20_sq;
    face_mid = p_opt20;
  }

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
  std::cout << "edge_mid01 = " << edge_mid01 << ", d^2 = " << d_mid01_sq << std::endl;
  std::cout << "edge_mid12 = " << edge_mid12 << ", d^2 = " << d_mid12_sq << std::endl;
  std::cout << "edge_mid20 = " << edge_mid20 << ", d^2 = " << d_mid20_sq << std::endl;

  std::cout << "p_opt01    = " << p_opt01    << ", d^2 = " << d_opt01_sq << std::endl;
  std::cout << "p_opt12    = " << p_opt12    << ", d^2 = " << d_opt12_sq << std::endl;
  std::cout << "p_opt20    = " << p_opt20    << ", d^2 = " << d_opt20_sq << std::endl;

  std::cout << "selected face_mid = " << face_mid << ", max d^2 = " << max_d_sq << std::endl;
#endif

  if(m_oracle.first_intersection(face_mid, closest_pt, steiner_point, m_offset))
  {
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
#endif
    return true;
  }

  Vector_3 unit = vector(closest_pt, neighbor_cc);
  unit = scale(unit, m_offset / approximate_sqrt(geom_traits().compute_squared_length_3_object()(unit)));
  steiner_point = translate(closest_pt, unit);

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
  std::cout << "Steiner found through neighboring tet intersecting the input: " << steiner_point << std::endl;
  std::cout << "Closest point: " << closest_pt << std::endl;
  std::cout << "Direction: " << vector(closest_pt, neighbor_cc) << std::endl;
#endif

  return true;
}
  }
    if(is_neighbor_cc_in_offset)
    {
      const Point_3& ch_cc = circumcenter(ch);

      // If the voronoi edge intersects the offset, the steiner point is the first intersection
      if(m_oracle.first_intersection(ch_cc, neighbor_cc, steiner_point, m_offset))
      {
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
        std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
#endif
        return true;
      }
    }

    Tetrahedron_with_outside_info<Geom_traits> tet(neighbor, geom_traits());
    if(m_oracle.do_intersect(tet))
    {
      // steiner point is the closest point on input from cell centroid with offset
      const Point_3 closest_pt = m_oracle.closest_point(neighbor_cc);
      CGAL_assertion(closest_pt != neighbor_cc);

      Vector_3 unit = vector(closest_pt, neighbor_cc);

      // PMP::internal::normalize() requires sqrt()
      unit = scale(unit, m_offset / approximate_sqrt(geom_traits().compute_squared_length_3_object()(unit)));
      steiner_point = translate(closest_pt, unit);

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
      std::cout << "Steiner found through neighboring tet intersecting the input: " << steiner_point << std::endl;
      std::cout << "Closest point: " << closest_pt << std::endl;
      std::cout << "Direction: " << vector(closest_pt, neighbor_cc) << std::endl;
#endif

      return true;
    }

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "No Steiner point" << std::endl;
#endif

    return false;
  }

  // normal steiner computation
  bool compute_steiner_point_test(const Cell_handle ch,
                             const Cell_handle neighbor,
                             Point_3& steiner_point) const
  {
    CGAL_precondition(!m_tr.is_infinite(neighbor));

    typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
    typename Geom_traits::Construct_vector_3 vector = geom_traits().construct_vector_3_object();
    typename Geom_traits::Construct_translated_point_3 translate = geom_traits().construct_translated_point_3_object();
    typename Geom_traits::Construct_scaled_vector_3 scale = geom_traits().construct_scaled_vector_3_object();
    typename Geom_traits::Construct_midpoint_3 construct_midpoint = geom_traits().construct_midpoint_3_object();
    typename Geom_traits::Compute_squared_distance_3 sq_dist = geom_traits().compute_squared_distance_3_object();
    auto cross = geom_traits().construct_cross_product_vector_3_object();
    auto sq_length = geom_traits().compute_squared_length_3_object();

    const Point_3& neighbor_cc = circumcenter(neighbor);
    const Ball_3 neighbor_cc_offset_ball = ball(neighbor_cc, m_sq_offset);
    const bool is_neighbor_cc_in_offset = m_oracle.do_intersect(neighbor_cc_offset_ball);


#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "Compute_steiner_point(" << &*ch << ", " << &*neighbor << ")" << std::endl;

    const Point_3& chc = circumcenter(ch);
    std::cout << "CH" << std::endl;
    std::cout << "\t" << ch->vertex(0)->point() << std::endl;
    std::cout << "\t" << ch->vertex(1)->point() << std::endl;
    std::cout << "\t" << ch->vertex(2)->point() << std::endl;
    std::cout << "\t" << ch->vertex(3)->point() << std::endl;
    std::cout << "cc: " << chc << std::endl;
    std::cout << "CC Distance to input: " << CGAL::squared_distance(chc, m_oracle.closest_point(chc)) << std::endl;

    std::cout << "NCH" << std::endl;
    std::cout << "\t" << neighbor->vertex(0)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(1)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(2)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(3)->point() << std::endl;
    std::cout << "ncc: " << neighbor_cc << std::endl;
    std::cout << "NC Distance to input: " << CGAL::squared_distance(neighbor_cc, m_oracle.closest_point(neighbor_cc)) << std::endl;
    std::cout << "is_neighbor_cc_in_offset = " << is_neighbor_cc_in_offset << std::endl;

    std::cout << "squared offset " << m_sq_offset << std::endl;
#endif

    // ch's circumcenter should not be within the offset volume
    CGAL_assertion_code(const Point_3& ch_cc = circumcenter(ch);)
    CGAL_assertion_code(const Ball_3 ch_cc_offset_ball = ball(ch_cc, m_sq_offset);)
    CGAL_assertion(!m_oracle.do_intersect(ch_cc_offset_ball));

  Point_3 pts[3];
  int fi = ch->index(neighbor);
  for(int j=0, k=0; j<4; ++j) {
    if(j != fi) pts[k++] = ch->vertex(j)->point();
  }

  // --- 2. Check Traversability ---
  Facet f = Facet(ch, fi);
  const bool traversible = less_squared_radius_of_min_empty_sphere(m_sq_alpha, f, m_tr);
  if (traversible) {
    return compute_steiner_point_normal(ch, neighbor, steiner_point);
    // std::cout << "face is traversible!" << std::endl;
  }
  if (!traversible) {
    // std::cout << "face is NOT traversible!" << std::endl;
    Tetrahedron_with_outside_info<Geom_traits> tet(neighbor, geom_traits());
    if(m_oracle.do_intersect(tet))
    {
      // steiner point is the closest point on input from cell centroid with offset
      const Point_3 closest_pt = m_oracle.closest_point(neighbor_cc);
      CGAL_assertion(closest_pt != neighbor_cc);

      // compute midpoints on triangle edges
      Point_3 edge_mid01 = construct_midpoint(pts[0], pts[1]);
      Point_3 edge_mid12 = construct_midpoint(pts[1], pts[2]);
      Point_3 edge_mid20 = construct_midpoint(pts[2], pts[0]);
      // compute closest points on input from the triangle edges midpoints
      const Point_3 closest_pt01 = m_oracle.closest_point(edge_mid01);
      const Point_3 closest_pt12 = m_oracle.closest_point(edge_mid12);
      const Point_3 closest_pt20 = m_oracle.closest_point(edge_mid20);

      // 2. Compute the distances from triangle edges midpoints to input
      FT d01_sq = sq_dist(edge_mid01, closest_pt01);
      FT d12_sq = sq_dist(edge_mid12, closest_pt12);
      FT d20_sq = sq_dist(edge_mid20, closest_pt20);

      // compute furthest point on triangle from midpoint
      Point_3 face_mid(
          (pts[0].x() + pts[1].x() + pts[2].x()) / 3.0,
          (pts[0].y() + pts[1].y() + pts[2].y()) / 3.0,
          (pts[0].z() + pts[1].z() + pts[2].z()) / 3.0
      );
      std::cout << "face_mid = " << face_mid << std::endl;

      if (d01_sq >= d12_sq && d01_sq >= d20_sq)
      {
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
        std::cout << "===================== EDGE 0 --> 1 ======================\n";
#endif
        std::cout << "================= EDCGE 0 --> 1 ===============" << std::endl;
        face_mid = edge_mid01;
        std::cout << "face_mid: " << face_mid << std::endl;
        compute_edge_optimum(pts[0], pts[1], edge_mid01, face_mid);
        const Point_3 closest_pt_opt01 = m_oracle.closest_point(face_mid);
        // 2. Compute the distances from triangle edges midpoints to input
        FT dopt_01= sq_dist(face_mid, closest_pt_opt01);
        if (d01_sq >= dopt_01) {
          face_mid = edge_mid01;
        }
        std::cout << "pt_0: " << pts[0] << std::endl;
        std::cout << "face_mid: " << face_mid << std::endl;
        std::cout << "pt_1: " << pts[1] << std::endl;
      }
      else if (d12_sq >= d20_sq)
      {
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
        std::cout << "===================== EDGE 1 --> 2 ======================\n";
#endif
        std::cout << "================= EDCGE 1 --> 2 ===============" << std::endl;
        face_mid = edge_mid12;
        std::cout << "face_mid: " << face_mid << std::endl;
        compute_edge_optimum(pts[1], pts[2], edge_mid12, face_mid);

        const Point_3 closest_pt_opt12 = m_oracle.closest_point(face_mid);
        // 2. Compute the distances from triangle edges midpoints to input
        FT dopt_12 = sq_dist(face_mid, closest_pt_opt12);
        if (d12_sq >= dopt_12) {
          face_mid = edge_mid12;
        }
        std::cout << "pt_1: " << pts[1] << std::endl;
        std::cout << "face_mid: " << face_mid << std::endl;
        std::cout << "pt_2: " << pts[2] << std::endl;
      }
      else
      {
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
        std::cout << "===================== EDGE 2 --> 0 ======================\n";
#endif
        std::cout << "================= EDGE 2 --> 0 ===============" << std::endl;
        face_mid = edge_mid20;
        std::cout << "face_mid: " << face_mid << std::endl;
        compute_edge_optimum(pts[2], pts[0], edge_mid20, face_mid);
        const Point_3 closest_pt_opt20 = m_oracle.closest_point(face_mid);
        // 2. Compute the distances from triangle edges midpoints to input
        FT dopt_20 = sq_dist(face_mid, closest_pt_opt20);
        if (d20_sq >= dopt_20) {
          face_mid = edge_mid20;
        }
        std::cout << "pt_2: " << pts[2] << std::endl;
        std::cout << "face_mid: " << face_mid << std::endl;
        std::cout << "pt_0: " << pts[0] << std::endl;
      }

      // Point_3 face_mid = CGAL::centroid(pts[0], pts[1], pts[2]);
      if(m_oracle.first_intersection(face_mid, closest_pt, steiner_point, m_offset))
      // if(m_oracle.first_intersection(closest_pt, face_mid, steiner_point, m_offset))
      {
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
        std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
#endif
        return true;
      }

      // normal logic
      std::cout << "uhhhh, function comes here???" << std::endl;
      Vector_3 unit = vector(closest_pt, neighbor_cc);
      // PMP::internal::normalize() requires sqrt()
      unit = scale(unit, m_offset / approximate_sqrt(geom_traits().compute_squared_length_3_object()(unit)));
      steiner_point = translate(closest_pt, unit);

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
      std::cout << "Steiner found through neighboring tet intersecting the input: " << steiner_point << std::endl;
      std::cout << "Closest point: " << closest_pt << std::endl;
      std::cout << "Direction: " << vector(closest_pt, neighbor_cc) << std::endl;
#endif

      return true;
    }
  }
    if(is_neighbor_cc_in_offset)
    {
      const Point_3& ch_cc = circumcenter(ch);

      // If the voronoi edge intersects the offset, the steiner point is the first intersection
      if(m_oracle.first_intersection(ch_cc, neighbor_cc, steiner_point, m_offset))
      {
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
        std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
#endif
        return true;
      }
    }

    Tetrahedron_with_outside_info<Geom_traits> tet(neighbor, geom_traits());
    if(m_oracle.do_intersect(tet))
    {
      // steiner point is the closest point on input from cell centroid with offset
      const Point_3 closest_pt = m_oracle.closest_point(neighbor_cc);
      CGAL_assertion(closest_pt != neighbor_cc);

      Vector_3 unit = vector(closest_pt, neighbor_cc);

      // PMP::internal::normalize() requires sqrt()
      unit = scale(unit, m_offset / approximate_sqrt(geom_traits().compute_squared_length_3_object()(unit)));
      steiner_point = translate(closest_pt, unit);

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
      std::cout << "Steiner found through neighboring tet intersecting the input: " << steiner_point << std::endl;
      std::cout << "Closest point: " << closest_pt << std::endl;
      std::cout << "Direction: " << vector(closest_pt, neighbor_cc) << std::endl;
#endif

      return true;
    }

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "No Steiner point" << std::endl;
#endif

    return false;
  }

private:
  // A permissive gate is a gate that we can traverse without checking its circumradius
  enum class Facet_status
  {
    IRRELEVANT = 0,
    HAS_INFINITE_NEIGHBOR, // the cell incident to the mirrored facet is infinite (permissive)
    SCAFFOLDING, // incident to a SEED or BBOX vertex (permissive)
    TRAVERSABLE
  };

  inline const char* get_status_message(const Facet_status status)
  {
    constexpr std::size_t status_count = 4;

    // Messages corresponding to Error_code list above. Must be kept in sync!
    static const char* message[status_count] = {
      "Irrelevant facet",
      "Facet incident to infinite neighbor",
      "Facet with a bbox/seed vertex",
      "Traversable facet"
    };

    const std::size_t status_id = static_cast<std::size_t>(status);
    if(status_id > status_count || status_id < 0)
      return "Unknown status";
    else
      return message[status_id];
  }

public:
  // @speed some decent time may be spent constructing Facet (pairs) for no reason as it's always
  // just to grab the .first and .second as soon as it's constructed, and not due to API requirements
  Facet_status facet_status(const Facet& f) const
  {
    CGAL_precondition(!m_tr.is_infinite(f));

#ifdef CGAL_AW3_DEBUG_FACET_STATUS
    std::cout << "facet status: "
              << m_tr.point(f.first, Triangulation::vertex_triple_index(f.second, 0)) << " "
              << m_tr.point(f.first, Triangulation::vertex_triple_index(f.second, 1)) << " "
              << m_tr.point(f.first, Triangulation::vertex_triple_index(f.second, 2)) << std::endl;
#endif

    // skip if neighbor is "outside" or infinite
    const Cell_handle ch = f.first;
    const int id = f.second;
    CGAL_precondition(ch->label() == Cell_label::INSIDE || ch->label() == Cell_label::OUTSIDE);

    if(!ch->is_outside()) // "inside" or "manifold"
    {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
      std::cout << "Facet is inside" << std::endl;
#endif
      return Facet_status::IRRELEVANT;
    }

    const Cell_handle nh = ch->neighbor(id);
    CGAL_precondition(ch->label() == Cell_label::INSIDE || ch->label() == Cell_label::OUTSIDE);

    if(m_tr.is_infinite(nh))
      return Facet_status::HAS_INFINITE_NEIGHBOR;

    if(nh->is_outside())
    {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
      std::cout << "Neighbor already outside" << std::endl;
#endif
      return Facet_status::IRRELEVANT;
    }

    // push if facet is connected to scaffolding vertices
    for(int i=0; i<3; ++i)
    {
      const Vertex_handle vh = ch->vertex(Triangulation::vertex_triple_index(id, i));
      if(vh->type() == AW3i::Vertex_type:: BBOX_VERTEX ||
         vh->type() == AW3i::Vertex_type:: SEED_VERTEX)
      {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
        std::cout << "Scaffolding facet due to vertex #" << i << std::endl;
#endif
        return Facet_status::SCAFFOLDING;
      }
    }

    // skip if f min empty sphere radius is smaller than alpha
    if(is_traversable(f))
    {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
      std::cout << "traversable" << std::endl;
#endif
      return Facet_status::TRAVERSABLE;
    }

#ifdef CGAL_AW3_DEBUG_FACET_STATUS
    std::cout << "not traversable" << std::endl;
#endif
    return Facet_status::IRRELEVANT;
  }

private:
  bool push_facet(const Facet& f)
  {
    CGAL_precondition(f.first->is_outside());

#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
    // skip if f is already in queue
    if(m_queue.contains_with_bounds_check(Gate(f)))
      return false;
#endif

    // @todo could avoid useless facet_status() calls by doing it after the zombie check
    // for the unsorted priority queue, but AFAIR, it doesn't save noticeable time (and that
    // increases the queue size).
    const Facet_status status = facet_status(f);
    if(status == Facet_status::IRRELEVANT)
      return false;

#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
    const FT sqr = smallest_squared_radius_3(f, m_tr);
    const bool is_permissive = (status == Facet_status::HAS_INFINITE_NEIGHBOR ||
                                status == Facet_status::SCAFFOLDING);
    m_queue.resize_and_push(Gate(f, sqr, is_permissive));
#else
    m_queue.push(Gate(f, m_tr));
#endif

#ifdef CGAL_AW3_DEBUG_QUEUE
    const Cell_handle ch = f.first;
    const int s = f.second;
    const Point_3& p0 = m_tr.point(ch, Triangulation::vertex_triple_index(s, 0));
    const Point_3& p1 = m_tr.point(ch, Triangulation::vertex_triple_index(s, 1));
    const Point_3& p2 = m_tr.point(ch, Triangulation::vertex_triple_index(s, 2));

    static int gid = 0;
    std::cout << "Queue insertion #" << gid++ << "\n"
              << "  ch = " << &*ch << " (" << m_tr.is_infinite(ch) << ") " << "\n"
              << "\t" << p0 << "\n\t" << p1 << "\n\t" << p2 << std::endl;
    std::cout << "  Status: " << get_status_message(status) << std::endl;
 #ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
    std::cout << "  SQR: " << sqr << std::endl;
    std::cout << "  Permissiveness: " << is_permissive << std::endl;

    CGAL_assertion(is_permissive || sqr >= m_sq_alpha);
 #endif // CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
#endif // CGAL_AW3_DEBUG_QUEUE

    return true;
  }

private:
  bool initialize(const double alpha,
                  const double offset,
                  const bool refining)
  {
#ifdef CGAL_AW3_DEBUG
    std::cout << "> Initialize..." << std::endl;
#endif

    const bool resuming = refining && (alpha == m_alpha) && (offset == m_offset);

#ifdef CGAL_AW3_DEBUG
    std::cout << "\tAlpha: " << alpha << std::endl;
    std::cout << "\tOffset: " << offset << std::endl;
    std::cout << "\tRefining? " << refining << std::endl;
    std::cout << "\tResuming? " << resuming << std::endl;
#endif

    if(!is_positive(alpha) || !is_positive(offset))
    {
#ifdef CGAL_AW3_DEBUG
      std::cerr << "Error: invalid input parameters: " << alpha << " and " << offset << std::endl;
#endif
      return false;
    }

#ifdef CGAL_AW3_DEBUG
    if(refining && alpha > m_alpha)
      std::cerr << "Warning: refining with an alpha greater than the last iteration's!" << std::endl;
    if(refining && offset != m_offset)
      std::cerr << "Warning: refining with a different offset value!" << std::endl;
#endif

    m_alpha = FT(alpha);
    m_sq_alpha = square(m_alpha);
    m_offset = FT(offset);
    m_sq_offset = square(m_offset);

    const Bbox_3 ib = m_oracle.bbox();
    m_xmid = 0.5 * (ib.xmin() + ib.xmax());
    m_ymid = 0.5 * (ib.ymin() + ib.ymax());
    m_zmid = 0.5 * (ib.zmin() + ib.zmax());

    // Resuming means that we do not need to re-initialize the queue: either we have finished
    // and there is nothing to do, or the interruption was due to a user callback in the visitor,
    // and we can resume with the current queue
    if(resuming)
    {
#ifdef CGAL_AW3_DEBUG
      std::cout << "Resuming with a queue of size: " << m_queue.size() << std::endl;
#endif

      reset_manifold_labels();
      return true;
    }

#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
    m_queue.clear();
#else
    m_queue = { };
#endif

    if(refining)
    {
      // If we are reusing the triangulation, change the label of the extra elements
      // that we have added to ensure a manifold result back to external ("manifold" -> "outside")
      reset_manifold_labels();

      return initialize_from_existing_triangulation();
    }
    else
    {
      m_tr.clear();

      insert_bbox_corners();

      if(m_seeds.empty())
        return initialize_from_infinity();
      else
        return initialize_with_cavities();
    }
  }

  template <typename Visitor>

// ooh found the alpha flood fill function :)
  // this function starts after initialisation
  //    -> builds/updates 3D triangulation
  //    -> computes offset surface
  //    -> inserts points
  //    -> prepares data structures
  bool alpha_flood_fill(Visitor& visitor)
  {
#ifdef CGAL_AW3_DEBUG
    std::cout << "> Flood fill..." << std::endl;
#endif

    visitor.on_flood_fill_begin(*this);

    // Explore all finite cells that are reachable from one of the initial outside cells.
    // 1 iteration per gate (gate = facet that has on one side a OUTSIDE and on the other side a INSIDE tetrahedral cell)
    // loop trhough these gates till this 'gate list' is empty
    while(!m_queue.empty())
    {
      if(!visitor.go_further(*this))
        return false;

#ifdef CGAL_AW3_DEBUG_QUEUE_PP
      check_queue_sanity();
#endif

      // const& to something that will be popped, but safe as `ch` && `id` are extracted before the pop
      // take highest priority gate (safest first -> based on alpha, istance, offset, etc.?)
      const Gate& gate = m_queue.top();

#ifndef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
      // face is no longer a valid gate; e.g. now both sides outside, or things changed after Steiner Point insertion
      if(gate.is_zombie())
      {
        m_queue.pop();
        continue;
      }
#endif

      // extract the facet and its two incident cells
      const Facet& f = gate.facet();
      CGAL_precondition(!m_tr.is_infinite(f));

      // 'ch' must be OUTSIDE
      const Cell_handle ch = f.first; // 'ch' = OUTSIDE cell on one side of the facet
      const int s = f.second; // 'nh' = cell on other side
      CGAL_precondition(ch->is_outside());

      const Cell_handle nh = ch->neighbor(s); // 'nh' is neighbor of 'ch' (they share a triangular face)
      // 'nh' must be labeled inside or outside (zombie) -> later popped
      CGAL_precondition(nh->label() == Cell_label::INSIDE || nh->label() == Cell_label::OUTSIDE);

#ifdef CGAL_AW3_DEBUG_QUEUE
      static int fid = 0;
      std::cout << m_tr.number_of_vertices() << " DT vertices" << std::endl;
      std::cout << m_queue.size() << " facets in the queue" << std::endl;
      std::cout << "Face " << fid++ << "\n"
                << "c = " << &*ch << " (" << m_tr.is_infinite(ch) << "), n = " << &*nh << " (" << m_tr.is_infinite(nh) << ")" << "\n"
                << m_tr.point(ch, Triangulation::vertex_triple_index(s, 0)) << "\n"
                << m_tr.point(ch, Triangulation::vertex_triple_index(s, 1)) << "\n"
                << m_tr.point(ch, Triangulation::vertex_triple_index(s, 2)) << std::endl;
      std::cout << "Priority: " << gate.priority() << " (sq alpha: " << m_sq_alpha << ")" << std::endl;
      std::cout << "Permissiveness: " << gate.is_permissive_facet() << std::endl;
#endif

#ifdef CGAL_AW3_DEBUG_DUMP_EVERY_STEP
      static int i = 0;
      std::string step_name = "results/steps/step_" + std::to_string(static_cast<int>(i)) + ".off";
      dump_triangulation_faces(step_name, true /*only_boundary_faces*/);

      std::string face_name = "results/steps/face_" + std::to_string(static_cast<int>(i++)) + ".xyz";
      std::ofstream face_out(face_name);
      face_out.precision(17);
      face_out << "3\n" << m_tr.point(ch, Triangulation::vertex_triple_index(s, 0)) << "\n"
                        << m_tr.point(ch, Triangulation::vertex_triple_index(s, 1)) << "\n"
                        << m_tr.point(ch, Triangulation::vertex_triple_index(s, 2)) << std::endl;
      face_out.close();
#endif

      // gate can be procesed and popped from the queue
      visitor.before_facet_treatment(*this, gate);

      m_queue.pop();

      // Infinite cells per definition OUTSIDE
      if(m_tr.is_infinite(nh))
      {
        nh->set_label(Cell_label::OUTSIDE);
#ifndef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
        nh->increment_erase_counter();
#endif
        continue;
      }

      // if compute_steiner_point returns true:
      //    -> non traversal cell: insert a steiner point
      Point_3 steiner_point;
      if(compute_steiner_point(ch, nh, steiner_point))
      {
//        std::cout << CGAL::abs(CGAL::approximate_sqrt(m_oracle.squared_distance(steiner_point)) - m_offset)
//                  << " vs " << 1e-2 * m_offset << std::endl;
        // for steiner point insertion, looks if steiner point is 'near enough' the offset surface
        // so actually just a sanity check if they inserted steiner point good enough
        //    -> e.g. floating point error
        CGAL_assertion(CGAL::abs(CGAL::approximate_sqrt(m_oracle.squared_distance(steiner_point)) - m_offset) <= 1e-2 * m_offset);

        // first destroy conflicting cells caused by the steiner point insertion
        // locate cells that are going to be destroyed and remove their facet from the queue
        int li, lj = 0;
        // multiple options for lt:
        //    -> lt == CELL -> steiner point lies inside cell
        //        -> insertion splits tetrahedron into 4 new tetrahedra
        //    -> lt == FACET -> steiner point lies on triangular face, li used (face index (0-3))
        //        -> point lies on face li of this tetrahedron
        //    -> lt == EDGE -> steiner point lies on edge, li and lj used
        //        -> point lies on the edge between vertices li and lj
        //    -> lt == VERTEX -> FORBIDDEN
        Locate_type lt;
        const Cell_handle conflict_cell = m_tr.locate(steiner_point, lt, li, lj, nh);
        CGAL_assertion(lt != Triangulation::VERTEX);

        // Using small vectors like in Triangulation_3 does not bring any runtime improvement
        std::vector<Facet> boundary_facets; // = boundary between conflict zone and the rest
        std::vector<Cell_handle> conflict_zone; // = tetrahedra whose circumspheres contain the new point -> will be removed
        boundary_facets.reserve(32);
        conflict_zone.reserve(32);

        m_tr.find_conflicts(steiner_point, conflict_cell,
                            std::back_inserter(boundary_facets),
                            std::back_inserter(conflict_zone));

#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
        // Purge the queue of facets that will be deleted/modified by the Steiner point insertion,
        // and which might have been gates
        for(const Cell_handle& cch : conflict_zone)
        {
          for(int i=0; i<4; ++i)
          {
            const Facet cf = std::make_pair(cch, i);
            if(m_queue.contains_with_bounds_check(Gate(cf)))
              m_queue.erase(Gate(cf));
          }
        }

        for(const Facet& f : boundary_facets)
        {
          const Facet mf = m_tr.mirror_facet(f); // boundary facets have incident cells in the CZ
          if(m_queue.contains_with_bounds_check(Gate(mf)))
            m_queue.erase(Gate(mf));
        }
#endif

        visitor.before_Steiner_point_insertion(*this, steiner_point);

        // Actual insertion of the Steiner point :)
        // We could use TDS functions to avoid recomputing the conflict zone, but in practice
        // it does not bring any runtime improvements
        // also passes lt, li, lj
        Vertex_handle vh = m_tr.insert(steiner_point, lt, conflict_cell, li, lj);
        vh->type() = AW3i::Vertex_type:: DEFAULT;

        visitor.after_Steiner_point_insertion(*this, vh);

        std::vector<Cell_handle> new_cells;
        new_cells.reserve(32);
        m_tr.incident_cells(vh, std::back_inserter(new_cells));
        // label new inserted cells: infinite OUTSIDE, otherwise INSIDE
        for(const Cell_handle& new_ch : new_cells)
        {
          // std::cout << "new cell has time stamp " << new_ch->time_stamp() << std::endl;
          new_ch->set_label(m_tr.is_infinite(new_ch) ? Cell_label::OUTSIDE : Cell_label::INSIDE);
        }

        // Push all new boundary facets to the queue.
        // It is not performed by looking at the facets on the boundary of the conflict zones
        // because we need to handle internal facets, infinite facets, and also more subtle changes
        // such as a new cell being marked inside which now creates a boundary
        // with its incident "outside" flagged cell.
        for(Cell_handle new_ch : new_cells)
        {
          for(int i=0; i<4; ++i)
          {
            // ensure ch is outside cell
            if(m_tr.is_infinite(new_ch, i))
              continue;

            const Cell_handle new_nh = new_ch->neighbor(i);
            if(new_nh->label() == new_ch->label()) // not on a boundary
              continue;

            const Facet boundary_f = std::make_pair(new_ch, i);
            if(new_ch->is_outside())
              push_facet(boundary_f);
            else
              push_facet(m_tr.mirror_facet(boundary_f));
          }
        }
      }
      else // no need for a Steiner point, carve through and continue
      {
        // gate is traversal -> so carving operation is performed
        nh->set_label(Cell_label::OUTSIDE);
#ifndef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
        nh->increment_erase_counter();
#endif

        // for each finite facet of neighbor, push it to the queue
        // carved into nh via gate face, so nh now outside, other 3 faces might now be boundary facets, so push to queue
        const int mi = m_tr.mirror_index(ch, s);
        for(int i=1; i<4; ++i)
        {
          const Facet neighbor_f = std::make_pair(nh, (mi+i)&3);
          push_facet(neighbor_f);
        }
      }
    } // while(!queue.empty())

    // if queue is emtpy carving operations end
    visitor.on_flood_fill_end(*this);

    // Check that no useful facet has been ignored
    CGAL_postcondition_code(for(auto fit=m_tr.finite_facets_begin(), fend=m_tr.finite_facets_end(); fit!=fend; ++fit) {)
    CGAL_postcondition_code(  Cell_handle ch = fit->first; Cell_handle nh = fit->first->neighbor(fit->second); )
    CGAL_postcondition_code(  if(ch->label() == nh->label()) continue;)
    CGAL_postcondition_code(  Facet f = *fit;)
    CGAL_postcondition_code(  if(ch->is_inside()) f = m_tr.mirror_facet(f);)
    CGAL_postcondition(       facet_status(f) == Facet_status::IRRELEVANT);
    CGAL_postcondition_code(})

    return true;
  }

  // Any outside cell that isn't reachable from infinity is a cavity that can be discarded.
  std::size_t purge_inner_connected_components()
  {
#ifdef CGAL_AW3_DEBUG
    std::cout << "> Purge inner islands..." << std::endl;

    std::size_t pre_counter = 0;
    for(Cell_handle ch : m_tr.all_cell_handles())
      if(ch->is_outside())
        ++pre_counter;
    std::cout << pre_counter << " / " << m_tr.all_cell_handles().size() << " (pre purge)" << std::endl;
#endif

    std::size_t label_change_counter = 0;

    std::stack<Cell_handle> cells_to_visit;

    if(!m_seeds.empty())
    {
      for(const Point_3& seed : m_seeds)
      {
        Locate_type lt;
        int li, lj;
        Cell_handle ch = m_tr.locate(seed, lt, li, lj);

        if(!ch->is_outside())
        {
          std::cerr << "Warning: cell containing seed is not outside?!" << std::endl;
          continue;
        }

        cells_to_visit.push(ch);
      }
    }
    else // typical flooding from outside
    {
      cells_to_visit.push(m_tr.infinite_vertex()->cell());
    }

    while(!cells_to_visit.empty())
    {
      Cell_handle curr_c = cells_to_visit.top();
      cells_to_visit.pop();

      CGAL_assertion(curr_c->is_outside());

      if(curr_c->tds_data().processed())
        continue;
      curr_c->tds_data().mark_processed();

      for(int j=0; j<4; ++j)
      {
        Cell_handle neigh_c = curr_c->neighbor(j);
        if(neigh_c->tds_data().processed() || !neigh_c->is_outside())
          continue;

        cells_to_visit.push(neigh_c);
      }
    }

    for(Cell_handle ch : m_tr.all_cell_handles())
    {
      if(ch->tds_data().is_clear() && ch->is_outside())
      {
        ch->set_label(Cell_label::INSIDE);
#ifndef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
        ch->increment_erase_counter();
#endif
        ++label_change_counter;
      }
    }

    // reset the conflict flags
    for(Cell_handle ch : m_tr.all_cell_handles())
      ch->tds_data().clear();

#ifdef CGAL_AW3_DEBUG
    std::size_t post_counter = 0;
    for(Cell_handle ch : m_tr.all_cell_handles())
      if(ch->is_outside())
        ++post_counter;
    std::cout << post_counter << " / " << m_tr.all_cell_handles().size() << " (pre purge)" << std::endl;

    std::cout << label_change_counter << " label changes" << std::endl;
#endif

    return label_change_counter;
  }

private:
  bool is_non_manifold(Vertex_handle v) const
  {
    CGAL_precondition(!m_tr.is_infinite(v));

    bool is_non_manifold = false;

    std::vector<Cell_handle> inc_cells;
    inc_cells.reserve(64);
    m_tr.incident_cells(v, std::back_inserter(inc_cells));

    // Flood one inside and outside CC within the cell umbrella of the vertex.
    // Process both an inside and an outside CC to also detect edge pinching.
    // If there are still unprocessed afterwards, there is a non-manifoldness issue.
    //
    // Squat the conflict cell data to mark visits
    Cell_handle inside_start = Cell_handle();
    Cell_handle outside_start = Cell_handle();

    for(Cell_handle ic : inc_cells)
    {
      ic->tds_data().clear();
      if(ic->is_outside())
        outside_start = ic;
      else if(inside_start == Cell_handle())
        inside_start = ic; // can be "inside" or "manifold"
    }

    // fully inside / outside
    if(inside_start == Cell_handle() || outside_start == Cell_handle())
      return false;

    std::stack<Cell_handle> cells_to_visit;
    cells_to_visit.push(inside_start);
    cells_to_visit.push(outside_start);
    while(!cells_to_visit.empty())
    {
      Cell_handle curr_c = cells_to_visit.top();
      cells_to_visit.pop();

      if(curr_c->tds_data().processed())
        continue;
      curr_c->tds_data().mark_processed();

      int v_pos = -1;
      CGAL_assertion_code(bool res = ) curr_c->has_vertex(v, v_pos);
      CGAL_assertion(res);

      for(int j=0; j<4; ++j)
      {
        if(j == v_pos)
          continue;

        Cell_handle neigh_c = curr_c->neighbor(j);
        CGAL_assertion(neigh_c->has_vertex(v));

        if(neigh_c->tds_data().processed() ||
           is_on_outside_boundary(neigh_c, curr_c)) // do not cross the boundary
        {
          continue;
        }

        cells_to_visit.push(neigh_c);
      }
    }

    for(Cell_handle ic : inc_cells)
    {
      if(ic->tds_data().is_clear()) // <=> unvisited cell
      {
        is_non_manifold = true;
        break;
      }
    }

    // reset the conflict flags
    for(Cell_handle ic : inc_cells)
      ic->tds_data().clear();

    return is_non_manifold;
  }

  bool is_non_manifold(Cell_handle c) const
  {
    CGAL_precondition(!m_tr.is_infinite(c));

    for(int i=0; i<4; ++i)
    {
      Vertex_handle v = c->vertex(i);
      if(is_non_manifold(v))
        return true;
    }

    return false;
  }

  bool is_manifold() const
  {
    // Not the best complexity, but it's not important: this function is purely for information
    // Better complexity --> see PMP::non_manifold_vertices + throw
    for(const Vertex_handle v : m_tr.finite_vertex_handles())
      if(is_non_manifold(v))
        return true;

    return false;
  }

  // Remove bbox vertices, if they are not necessary (i.e., no "inside" incident cell)
  // This is to try and avoid having long tets with bbox vertices being tagged "inside" as part
  // of the manifold re-tagging
  bool remove_bbox_vertices()
  {
    bool do_remove = true;
    auto vit = m_tr.finite_vertices_begin();
    for(std::size_t i=0; i<8; ++i)
    {
      Vertex_handle v = vit++;

      std::vector<Cell_handle> inc_cells;
      inc_cells.reserve(64);
      m_tr.finite_incident_cells(v, std::back_inserter(inc_cells));

      for(Cell_handle c : inc_cells)
      {
        if(!c->is_outside())
        {
          do_remove = false;
          break;
        }
      }

      if(!do_remove)
        break;
    }

    std::cout << "Removing bbox vertices? " << do_remove << std::endl;
    if(!do_remove)
      return false;

    vit = m_tr.finite_vertices_begin();
    for(std::size_t i=0; i<8; ++i)
    {
      Vertex_handle v = vit++;
      m_tr.remove(v);
    }

    return true;
  }

public:
  // Not the best complexity, but it's very cheap compared to the rest of the algorithm.
  void make_manifold()
  {
#ifdef CGAL_AW3_DEBUG
    std::cout << "> Make manifold..." << std::endl;

    auto wrap_volume = [&]()
    {
      FT vol = 0;
      for(Cell_handle ch : m_tr.finite_cell_handles())
        if(!ch->is_outside())
          vol += volume(m_tr.point(ch, 0), m_tr.point(ch, 1), m_tr.point(ch, 2), m_tr.point(ch, 3));

      return vol;
    };

 #ifdef CGAL_AW3_DEBUG_DUMP_INTERMEDIATE_WRAPS
    dump_triangulation_faces("carved_tr.off", true /*only_boundary_faces*/);
 #endif

    FT base_vol = wrap_volume();
    if(!is_positive(base_vol))
      std::cerr << "Warning: wrap with non-positive volume?" << std::endl;
#endif

    // This ends up more harmful than useful after the priority queue has been introduced since
    // it usually results in a lot of flat cells into the triangulation, which then get added
    // to the mesh during manifoldness fixing.
//    remove_bbox_vertices();

    std::stack<Vertex_handle> non_manifold_vertices; // @todo sort somehow?
    for(Vertex_handle v : m_tr.finite_vertex_handles())
    {
      if(is_non_manifold(v))
      {
#ifdef CGAL_AW3_DEBUG_MANIFOLDNESS_PP
        std::cout << v->point() << " is non-manifold" << std::endl;
#endif
        non_manifold_vertices.push(v);
      }
    }

    // Some lambdas for the comparer
    auto has_scaffolding_vertex = [](Cell_handle c) -> bool
    {
      for(int i=0; i<4; ++i)
      {
        if(c->vertex(i)->type() == AW3i::Vertex_type:: BBOX_VERTEX ||
           c->vertex(i)->type() == AW3i::Vertex_type:: SEED_VERTEX)
        {
          return true;
        }
      }

      return false;
    };

    // This originally seemed like a good idea, but in the end it can have strong cascading issues,
    // whereas some cells with much smaller volume could have solved the non-manifoldness.
//    auto is_on_boundary = [](Cell_handle c, int i) -> bool
//    {
//      return is_on_outside_boundary(c, c->neighbor(i));
//    };
//
//    auto count_boundary_facets = [&](Cell_handle c, Vertex_handle v) -> int
//    {
//      const int v_index_in_c = c->index(v);
//      int boundary_facets = 0;
//      for(int i=0; i<3; ++i) // also take into account the opposite facet?
//      {
//        if(i == v_index_in_c)
//          continue;
//
//        boundary_facets += is_on_boundary(c, i);
//      }
//
//      return boundary_facets;
//    };

    // Experimentally, longest edge works better
//    auto sq_circumradius = [&](Cell_handle c) -> FT
//    {
//      const Point_3& cc = circumcenter(c);
//      return geom_traits().compute_squared_distance_3_object()(m_tr.point(c, 0), cc);
//    };

    // The reasoning behind using longest edge rather than volume is that we want to avoid
    // spikes (which would have a small volume), and can often happen since we do not spend
    // any care on the quality of tetrahedra.
    auto sq_longest_edge = [&](Cell_handle c) -> FT
    {
      return (std::max)({ squared_distance(m_tr.point(c, 0), m_tr.point(c, 1)),
                          squared_distance(m_tr.point(c, 0), m_tr.point(c, 2)),
                          squared_distance(m_tr.point(c, 0), m_tr.point(c, 3)),
                          squared_distance(m_tr.point(c, 1), m_tr.point(c, 2)),
                          squared_distance(m_tr.point(c, 1), m_tr.point(c, 3)),
                          squared_distance(m_tr.point(c, 2), m_tr.point(c, 3)) });
    };

#ifdef CGAL_AW3_DEBUG_MANIFOLDNESS_PP
    std::cout << non_manifold_vertices.size() << " initial NMV" << std::endl;
#endif

    while(!non_manifold_vertices.empty())
    {
#ifdef CGAL_AW3_DEBUG_MANIFOLDNESS_PP
      std::cout << non_manifold_vertices.size() << " NMV in queue" << std::endl;
#endif

      Vertex_handle v = non_manifold_vertices.top();
      non_manifold_vertices.pop();

      if(!is_non_manifold(v))
        continue;

      // Prioritize:
      // - cells without bbox vertices
      // - small cells when equal number of boundary facets
      //
      // Note that these are properties that do not depend on the cell labels, and so we only need
      // to sort once. However, if a criterion such as the number of incident inside cells were added,
      // we would need to sort after each modification of "inside"/"outside" labels.
      auto comparer = [&](Cell_handle l, Cell_handle r) -> bool
      {
        CGAL_precondition(!m_tr.is_infinite(l) && !m_tr.is_infinite(r));

        if(has_scaffolding_vertex(l) != has_scaffolding_vertex(r))
          return has_scaffolding_vertex(r);

        return sq_longest_edge(l) < sq_longest_edge(r);
      };

      std::vector<Cell_handle> inc_cells;
      inc_cells.reserve(64);
      m_tr.finite_incident_cells(v, std::back_inserter(inc_cells));

      std::vector<Cell_handle> finite_outside_inc_cells;
      finite_outside_inc_cells.reserve(64);
      std::copy_if(inc_cells.begin(), inc_cells.end(), std::back_inserter(finite_outside_inc_cells),
                   [&](Cell_handle c) -> bool { return !m_tr.is_infinite(c) && c->is_outside(); });

      // 'std::stable_sort' to have determinism without having to write something like:
      //     if(longest_edge(l) == longest_edge(r)) return ...
      // in the comparer. It's almost always a small range, so the extra cost does not matter.
      std::stable_sort(finite_outside_inc_cells.begin(), finite_outside_inc_cells.end(), comparer);

      for(Cell_handle ic : finite_outside_inc_cells)
      {
        CGAL_assertion(!m_tr.is_infinite(ic) && ic->is_outside());

        // This is where new material is added
        ic->set_label(Cell_label::MANIFOLD);

#ifdef CGAL_AW3_DEBUG_DUMP_EVERY_STEP
        static int i = 0;
        std::string step_name = "results/steps_manifold/step" + std::to_string(static_cast<int>(i++)) + ".off";
        dump_triangulation_faces(step_name, true /*only_boundary_faces*/);
#endif

        // @speed could update the manifold status while tagging
        if(!is_non_manifold(v))
          break;
      }

      CGAL_assertion(!is_non_manifold(v));

      // Check if the new material has not created a non-manifold configuration.
      // @speed this could be done on only the vertices of cells whose labels have changed.
      std::vector<Vertex_handle> adj_vertices;
      adj_vertices.reserve(64);
      m_tr.finite_adjacent_vertices(v, std::back_inserter(adj_vertices));

      for(Vertex_handle nv : adj_vertices)
        if(is_non_manifold(nv))
          non_manifold_vertices.push(nv);
    }

    CGAL_assertion_code(for(Vertex_handle v : m_tr.finite_vertex_handles()))
    CGAL_assertion(!is_non_manifold(v));

#ifdef CGAL_AW3_DEBUG
    std::size_t nm_cells_counter = 0;
    for(Cell_handle ch : m_tr.all_cell_handles())
      if(ch->label() == Cell_label::MANIFOLD)
        ++nm_cells_counter;
    std::cout << "Number of added cells: " << nm_cells_counter << std::endl;

    if(!is_zero(base_vol))
    {
      const FT manifold_vol = wrap_volume();
      const FT ratio = manifold_vol / base_vol;

      std::cout << "Volumes post-manifoldness fix:\n"
                << "before: " << base_vol << "\n"
                << "after:  " << manifold_vol << "\n"
                << "ratio:  " << ratio << std::endl;
      if(ratio > 1.1) // more than 10% extra volume
        std::cerr << "Warning: large increase of volume after manifoldness resolution" << std::endl;
    }
#endif
  }

private:
  void check_queue_sanity()
  {
    std::cout << "\t~~~ Check queue sanity ~~~" << std::endl;

    std::vector<Gate> queue_gates;
    Gate previous_top_gate = m_queue.top();
    while(!m_queue.empty())
    {
      const Gate& current_gate = m_queue.top();
      queue_gates.push_back(current_gate);
      const Facet& current_f = current_gate.facet();
      const Cell_handle ch = current_f.first;
      const int id = current_f.second;
      const Point_3& p0 = m_tr.point(ch, Triangulation::vertex_triple_index(id, 0));
      const Point_3& p1 = m_tr.point(ch, Triangulation::vertex_triple_index(id, 1));
      const Point_3& p2 = m_tr.point(ch, Triangulation::vertex_triple_index(id, 2));

      std::cout << "Facet with VID " << get(Gate_ID_PM<Triangulation>(), current_gate) << "\n";
      std::cout << "\t" << p0 << "\n\t" << p1 << "\n\t" << p2 << "\n";

#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
      std::cout << "  Permissiveness: " << current_gate.is_permissive_facet() << "\n";
      std::cout << "  SQR: " << geom_traits().compute_squared_radius_3_object()(p0, p1, p2) << "\n";
      std::cout << "  Priority " << current_gate.priority() << std::endl;

      if(Less_gate()(current_gate, previous_top_gate))
        std::cerr << "Error: current gate has higher priority than the previous top" << std::endl;

      previous_top_gate = current_gate;
#else
      if(current_gate.is_zombie())
        std::cout << "Gate is a zombie!" << std::endl;
#endif

      m_queue.pop();
    }

    std::cout << "\t~~~ End queue sanity check ~~~" << std::endl;

    // Rebuild
    CGAL_assertion(m_queue.empty());
    for(const auto& g : queue_gates)
      m_queue.push(g); // no need for a resize here since the vector capacity is unchanged
  }

  void dump_triangulation_faces(const std::string filename,
                                bool only_boundary_faces = false)
  {
    std::stringstream vertices_ss;
    vertices_ss.precision(17);

    std::stringstream facets_ss;
    facets_ss.precision(17);

    std::unordered_map<Vertex_handle, std::size_t> vertex_to_id;
    std::size_t nv = 0;
    std::size_t nf = 0;

    for(auto fit=m_tr.finite_facets_begin(), fend=m_tr.finite_facets_end(); fit!=fend; ++fit)
    {
      Cell_handle ch = fit->first;
      int s = fit->second;

      Cell_handle nh = ch->neighbor(s);
      if(only_boundary_faces && ch->label() == nh->label())
        continue;

      std::array<std::size_t, 3> ids;
      for(std::size_t pos=0; pos<3; ++pos)
      {
        Vertex_handle v = ch->vertex((s+pos+1)&3);
        auto insertion_res = vertex_to_id.emplace(v, nv);
        if(insertion_res.second)
        {
          vertices_ss << m_tr.point(v) << "\n";
          ++nv;
        }

        ids[pos] = insertion_res.first->second;
      }

      facets_ss << "3 " << ids[0] << " " << ids[1] << " " << ids[2] << "\n";
      ++nf;
    }

    std::ofstream out(filename.c_str());
    out << "OFF\n" << nv << " " << nf << " 0\n";
    out << vertices_ss.str() << "\n" << facets_ss.str() << std::endl;
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_ALPHA_WRAP_3_H
