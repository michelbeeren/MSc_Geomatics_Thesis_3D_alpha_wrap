// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL: https://github.com/CGAL/cgal/blob/v6.1/Alpha_wrap_3/include/CGAL/Alpha_wrap_3/internal/Triangle_mesh_oracle.h $
// $Id: include/CGAL/Alpha_wrap_3/internal/Triangle_mesh_oracle.h b26b07a1242 $
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_TRIANGLE_MESH_ORACLE_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_TRIANGLE_MESH_ORACLE_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_AABB_geom_traits.h>
#include <CGAL/Alpha_wrap_3/internal/Oracle_base.h>
#include <CGAL/Alpha_wrap_3/internal/splitting_helper.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

#include <algorithm>
#include <iostream>
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

// Just some typedefs for readability in the main oracle class
template <typename GT_>
struct TM_oracle_traits
{
  using Geom_traits = Alpha_wrap_AABB_geom_traits<GT_>; // Wrap the kernel to add Ball_3 + custom Do_intersect_3

  using Point_3 = typename Geom_traits::Point_3;
  using AABB_traits = typename AABB_tree_splitter_traits<Point_3, Geom_traits>::AABB_traits;
  using AABB_tree = typename AABB_tree_splitter_traits<Point_3, Geom_traits>::AABB_tree;
};

// @speed could do a partial specialization 'subdivide = false' with simpler code for speed?
template <typename GT_,
          typename BaseOracle = int,
          bool subdivide = true>
class Triangle_mesh_oracle
  : // this is the base that handles calls to the AABB tree
    public AABB_tree_oracle<typename TM_oracle_traits<GT_>::Geom_traits,
                            typename TM_oracle_traits<GT_>::AABB_tree,
                            typename std::conditional<
                                      /*condition*/subdivide,
                                      /*true*/Splitter_traversal_traits<typename TM_oracle_traits<GT_>::AABB_traits>,
                                      /*false*/Default_traversal_traits<typename TM_oracle_traits<GT_>::AABB_traits> >::type,
                            BaseOracle>,
    // this is the base that handles splitting input faces and inserting them into the AABB tree
    public AABB_tree_oracle_splitter<subdivide,
                                     typename TM_oracle_traits<GT_>::Point_3,
                                     typename TM_oracle_traits<GT_>::Geom_traits>
{
  using TMOT = TM_oracle_traits<GT_>;
  using Base_GT = GT_;

public:
  using Geom_traits = typename TMOT::Geom_traits;

private:
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Segment_3 = typename Geom_traits::Segment_3;
  using Triangle_3 = typename Geom_traits::Triangle_3;

  using AABB_traits = typename TMOT::AABB_traits;
  using AABB_tree = typename TMOT::AABB_tree;
  using AABB_traversal_traits = typename std::conditional<
                                  /*condition*/subdivide,
                                  /*true*/Splitter_traversal_traits<AABB_traits>,
                                  /*false*/Default_traversal_traits<AABB_traits> >::type;

  using Oracle_base = AABB_tree_oracle<Geom_traits, AABB_tree, AABB_traversal_traits, BaseOracle>;
  using Splitter_base = AABB_tree_oracle_splitter<subdivide, Point_3, Geom_traits>;
  using Primitive_id = typename Oracle_base::Primitive_id;

  struct Face_info
  {
    Point_3 a, b, c;

    bool e0_has_neighbor = false; // edge (a, b)
    Point_3 e0_n0, e0_n1, e0_n2;

    bool e1_has_neighbor = false; // edge (b, c)
    Point_3 e1_n0, e1_n1, e1_n2;

    bool e2_has_neighbor = false; // edge (c, a)
    Point_3 e2_n0, e2_n1, e2_n2;
  };

  std::vector<Face_info> m_face_info_by_fid;

public:
  // Constructors
  //
  // When using this constructor (and thus doing actual splitting), note that the oracle
  // will be adapted to this particular 'alpha', and so when calling again AW3(other_alpha)
  // the oracle might not have performed a split that is adapted to this other alpha value.
  Triangle_mesh_oracle(const double alpha,
                       const BaseOracle& base_oracle = BaseOracle(),
                       const Base_GT& gt = Base_GT())
    : Oracle_base(base_oracle, gt), Splitter_base(alpha)
  {
    Splitter_base::initialize_tree_property_maps(this->tree());
  }

  Triangle_mesh_oracle(const double alpha,
                       const Base_GT& gt,
                       const BaseOracle& base_oracle = BaseOracle())
    : Triangle_mesh_oracle(alpha, base_oracle, gt)
  { }

 Triangle_mesh_oracle(const BaseOracle& base_oracle,
                      const Base_GT& gt = Base_GT())
   : Triangle_mesh_oracle(0. /*alpha*/, base_oracle, gt)
 { }

 Triangle_mesh_oracle(const Base_GT& gt,
                      const BaseOracle& base_oracle = BaseOracle())
   : Triangle_mesh_oracle(0. /*alpha*/, base_oracle, gt)
 { }

 Triangle_mesh_oracle()
   : Triangle_mesh_oracle(0. /*alpha*/, BaseOracle(), Base_GT())
 { }

  template <typename Result>
  bool analyze_closest_point_and_primitive(const Point_3& closest_pt,
                                           const Primitive_id& primitive_id,
                                           Result& out) const
  {
    if(primitive_id.second >= m_face_info_by_fid.size())
      return false;

    const Face_info& f = m_face_info_by_fid[primitive_id.second];

    out.closest_point = closest_pt;
    out.t0_a = f.a;
    out.t0_b = f.b;
    out.t0_c = f.c;
    out.has_second_face = false;
    out.has_shared_edge_points = false;

    const FT eps = FT(1e-12);
    const auto sqd = this->geom_traits().compute_squared_distance_3_object();

    const auto same_point =
      [&](const Point_3& p, const Point_3& q) -> bool
      {
        return sqd(p, q) <= eps;
      };

    const auto set_location =
      [&](const int location_id)
      {
        out.location = static_cast<decltype(out.location)>(location_id);
      };

    const auto set_shared_edge_points =
      [&](const Point_3& a,
          const Point_3& b,
          const Point_3& c,
          const Point_3& n0,
          const Point_3& n1,
          const Point_3& n2)
      {
        out.shared_edge_a = a;      // A
        out.shared_edge_b = b;      // B
        out.triangle_1_opposite = c; // C

        // D = neighbor-triangle vertex that is not on shared edge A-B.
        if(!same_point(n0, a) && !same_point(n0, b))
        {
          out.triangle_2_opposite = n0;
          out.has_shared_edge_points = true;
          return;
        }
        if(!same_point(n1, a) && !same_point(n1, b))
        {
          out.triangle_2_opposite = n1;
          out.has_shared_edge_points = true;
          return;
        }
        if(!same_point(n2, a) && !same_point(n2, b))
        {
          out.triangle_2_opposite = n2;
          out.has_shared_edge_points = true;
          return;
        }

        out.has_shared_edge_points = false;
      };

    if(same_point(closest_pt, f.a) || same_point(closest_pt, f.b) || same_point(closest_pt, f.c))
    {
      set_location(2); // ON_VERTEX
      return true;
    }

    const Segment_3 s0(f.a, f.b);
    const Segment_3 s1(f.b, f.c);
    const Segment_3 s2(f.c, f.a);

    int edge_idx = -1;
    if(s0.has_on(closest_pt))
      edge_idx = 0;
    else if(s1.has_on(closest_pt))
      edge_idx = 1;
    else if(s2.has_on(closest_pt))
      edge_idx = 2;

    if(edge_idx == -1)
    {
      set_location(0); // ON_FACE
      return true;
    }

    set_location(1); // ON_EDGE

    if(edge_idx == 0 && f.e0_has_neighbor)
    {
      out.has_second_face = true;
      out.t1_a = f.e0_n0;
      out.t1_b = f.e0_n1;
      out.t1_c = f.e0_n2;
      set_shared_edge_points(f.a, f.b, f.c, f.e0_n0, f.e0_n1, f.e0_n2);
      return true;
    }

    if(edge_idx == 1 && f.e1_has_neighbor)
    {
      out.has_second_face = true;
      out.t1_a = f.e1_n0;
      out.t1_b = f.e1_n1;
      out.t1_c = f.e1_n2;
      set_shared_edge_points(f.b, f.c, f.a, f.e1_n0, f.e1_n1, f.e1_n2);
      return true;
    }

    if(edge_idx == 2 && f.e2_has_neighbor)
    {
      out.has_second_face = true;
      out.t1_a = f.e2_n0;
      out.t1_b = f.e2_n1;
      out.t1_c = f.e2_n2;
      set_shared_edge_points(f.c, f.a, f.b, f.e2_n0, f.e2_n1, f.e2_n2);
      return true;
    }

    return true;
  }

public:
  template <typename TriangleMesh,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_triangle_mesh(const TriangleMesh& tmesh,
                         const CGAL_NP_CLASS& np = CGAL::parameters::default_values())
  {
    using parameters::get_parameter;
    using parameters::choose_parameter;

    using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

    using VPM = typename GetVertexPointMap<TriangleMesh>::const_type;
    using Point_ref = typename boost::property_traits<VPM>::reference;

    CGAL_precondition(CGAL::is_triangle_mesh(tmesh));

    if(is_empty(tmesh))
    {
#ifdef CGAL_AW3_DEBUG
      std::cout << "Warning: Input is empty (TM)" << std::endl;
#endif
      return;
    }

#ifdef CGAL_AW3_DEBUG
    std::cout << "Insert into AABB tree (faces)..." << std::endl;
#endif

    VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tmesh));
    static_assert(std::is_same<typename boost::property_traits<VPM>::value_type, Point_3>::value);

    Splitter_base::reserve(num_faces(tmesh));

    for(face_descriptor f : faces(tmesh))
    {
      if(Polygon_mesh_processing::is_degenerate_triangle_face(f, tmesh, np))
      {
#ifdef CGAL_AW3_DEBUG
        std::cerr << "Warning: ignoring degenerate face " << f << std::endl;
#endif
        continue;
      }

      const Point_ref p0 = get(vpm, source(halfedge(f, tmesh), tmesh));
      const Point_ref p1 = get(vpm, target(halfedge(f, tmesh), tmesh));
      const Point_ref p2 = get(vpm, target(next(halfedge(f, tmesh), tmesh), tmesh));

      const Triangle_3 tr = this->geom_traits().construct_triangle_3_object()(p0, p1, p2);

      const std::size_t fid = Splitter_base::fid;
      CGAL_assertion(fid == m_face_info_by_fid.size());
      m_face_info_by_fid.push_back(make_face_info(f, tmesh, vpm));

      Splitter_base::split_and_insert_datum(tr, this->tree(), this->geom_traits());
    }

    // Manually constructing it here purely for profiling reasons: if we keep the lazy approach,
    // it will be done at the first treatment of a facet that needs a Steiner point.
    // So if one wanted to bench the flood fill runtime, it would be skewed by the time it takes
    // to accelerate the tree.
    this->tree().accelerate_distance_queries();

#ifdef CGAL_AW3_DEBUG
    std::cout << "Tree: " << this->tree().size() << " primitives (" << num_faces(tmesh) << " faces in input)" << std::endl;
#endif
  }

private:
  template <typename TriangleMesh, typename VPM>
  Face_info make_face_info(const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                           const TriangleMesh& tmesh,
                           VPM vpm) const
  {
    using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
    using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

    Face_info out;

    const halfedge_descriptor h0 = halfedge(f, tmesh);
    const halfedge_descriptor h1 = next(h0, tmesh);
    const halfedge_descriptor h2 = next(h1, tmesh);

    out.a = get(vpm, source(h0, tmesh));
    out.b = get(vpm, target(h0, tmesh));
    out.c = get(vpm, target(h1, tmesh));

    const auto set_neighbor_triangle =
      [&](const halfedge_descriptor edge_h,
          bool& has_neighbor,
          Point_3& n0,
          Point_3& n1,
          Point_3& n2)
      {
        const halfedge_descriptor opp = opposite(edge_h, tmesh);
        const face_descriptor nfd = face(opp, tmesh);

        if(nfd == boost::graph_traits<TriangleMesh>::null_face())
        {
          has_neighbor = false;
          return;
        }

        has_neighbor = true;
        const halfedge_descriptor nh0 = halfedge(nfd, tmesh);
        const halfedge_descriptor nh1 = next(nh0, tmesh);

        n0 = get(vpm, source(nh0, tmesh));
        n1 = get(vpm, target(nh0, tmesh));
        n2 = get(vpm, target(nh1, tmesh));
      };

    set_neighbor_triangle(h0, out.e0_has_neighbor, out.e0_n0, out.e0_n1, out.e0_n2);
    set_neighbor_triangle(h1, out.e1_has_neighbor, out.e1_n0, out.e1_n1, out.e1_n2);
    set_neighbor_triangle(h2, out.e2_has_neighbor, out.e2_n0, out.e2_n1, out.e2_n2);

    return out;
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_TRIANGLE_MESH_ORACLE_H
