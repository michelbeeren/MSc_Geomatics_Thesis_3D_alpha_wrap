// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL: https://github.com/CGAL/cgal/blob/v6.1/Alpha_wrap_3/include/CGAL/Alpha_wrap_3/internal/Oracle_base.h $
// $Id: include/CGAL/Alpha_wrap_3/internal/Oracle_base.h b26b07a1242 $
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_ORACLE_BASE_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_ORACLE_BASE_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Alpha_wrap_3/internal/offset_intersection.h>

#include <CGAL/AABB_tree/internal/AABB_traversal_traits.h>
#include <CGAL/AABB_tree/internal/Primitive_helper.h>
#include <CGAL/Default.h>

#include <algorithm>
#include <memory>
#include <type_traits>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

template <typename AABBTraits>
struct Default_traversal_traits
{
  using Projection_traits = CGAL::internal::AABB_tree::Projection_traits<AABBTraits>;

  template <typename Query>
  using Do_intersect_traits = CGAL::internal::AABB_tree::Do_intersect_traits<AABBTraits, Query>;

  template <typename Query>
  using First_intersection_traits = CGAL::internal::AABB_tree::First_intersection_traits<AABBTraits, Query>;
};

// Factorize the implementation of the functions calling the AABB tree
template <typename AABBTree,
          typename AABBTraversalTraits>
struct AABB_tree_oracle_helper
{
  using Self = AABB_tree_oracle_helper<AABBTree, AABBTraversalTraits>;

  using AABB_traits = typename AABBTree::AABB_traits;
  using GT = typename AABB_traits::Geom_traits;

  using FT = typename AABB_traits::FT;
  using Point_3 = typename AABB_traits::Point;

  template <typename Query>
  static bool do_intersect(const Query& query,
                           const AABBTree& tree)
  {
    CGAL_precondition(!tree.empty());

    using Do_intersect_traits = typename AABBTraversalTraits::template Do_intersect_traits<Query>;

    Do_intersect_traits traversal_traits(tree.traits());
    tree.traversal(query, traversal_traits);
    return traversal_traits.is_intersection_found();
  }

  static Point_3 closest_point(const Point_3& p,
                               const AABBTree& tree)
  {
    CGAL_precondition(!tree.empty());

    using Projection_traits = typename AABBTraversalTraits::Projection_traits;

    const auto& hint = tree.best_hint(p);

    Projection_traits projection_traits(hint.first, hint.second, tree.traits());
    tree.traversal(p, projection_traits);
    return projection_traits.closest_point();
  }

  static Point_3 closest_point_edge_vertex(const Point_3& p,
                                           const AABBTree& tree)
  {
    CGAL_precondition(!tree.empty());

    using Primitive = typename AABB_traits::Primitive;
    using Primitive_id = typename Primitive::Id;
    using Node = ::CGAL::AABB_node<AABB_traits>;

    class Edge_vertex_projection_traits
    {
    public:
      Edge_vertex_projection_traits(const Point_3& hint,
                                    const Primitive_id& hint_primitive,
                                    const AABB_traits& traits)
        : m_closest_point(hint),
          m_closest_primitive(hint_primitive),
          m_traits(traits)
      { }

      constexpr bool go_further() const { return true; }

      void intersection(const Point_3& query, const Primitive& primitive)
      {
        const Point_3 new_closest_point =
          Self::closest_point_on_edges_vertices(query, primitive, m_closest_point, m_traits);

        if(!m_traits.equal_object()(new_closest_point, m_closest_point))
        {
          m_closest_primitive = primitive.id();
          m_closest_point = new_closest_point;
        }
      }

      bool do_intersect(const Point_3& query, const Node& node) const
      {
        return m_traits.compare_distance_object()
          (query, node.bbox(), m_closest_point) == CGAL::SMALLER;
      }

      Point_3 closest_point() const { return m_closest_point; }

    private:
      Point_3 m_closest_point;
      Primitive_id m_closest_primitive;
      const AABB_traits& m_traits;
    };

    const auto& hint = tree.best_hint(p);
    const Primitive hint_primitive(hint.second);
    const Point_3 edge_vertex_hint =
      Self::closest_point_on_edges_vertices_unbounded(p, hint_primitive, tree.traits());

    Edge_vertex_projection_traits projection_traits(edge_vertex_hint, hint.second, tree.traits());
    tree.traversal(p, projection_traits);
    return projection_traits.closest_point();
  }

  static FT squared_distance(const Point_3& p,
                             const AABBTree& tree)
  {
    CGAL_precondition(!tree.empty());

    const Point_3 closest = Self::closest_point(p, tree);
    return tree.traits().squared_distance_object()(p, closest);
  }

private:
  template <typename Primitive>
  static Point_3 closest_point_on_edges_vertices_unbounded(const Point_3& query,
                                                           const Primitive& primitive,
                                                           const AABB_traits& traits)
  {
    using Datum = typename Primitive::Datum;
    using Triangle_3 = typename GT::Triangle_3;

    const GT gt;
    const auto datum = CGAL::internal::Primitive_helper<AABB_traits>::get_datum(primitive, traits);

    if constexpr(std::is_same<Datum, Triangle_3>::value)
    {
      const auto vertex = gt.construct_vertex_3_object();
      const auto segment = gt.construct_segment_3_object();
      const auto project = gt.construct_projected_point_3_object();
      const auto compare_distance = gt.compare_distance_3_object();

      const Point_3& p0 = vertex(datum, 0);
      const Point_3& p1 = vertex(datum, 1);
      const Point_3& p2 = vertex(datum, 2);

      const Point_3 c01 = project(segment(p0, p1), query);
      const Point_3 c12 = project(segment(p1, p2), query);
      const Point_3 c20 = project(segment(p2, p0), query);

      Point_3 candidate = c01;
      if(compare_distance(query, c12, candidate) == CGAL::SMALLER)
        candidate = c12;
      if(compare_distance(query, c20, candidate) == CGAL::SMALLER)
        candidate = c20;

      return candidate;
    }

    return gt.construct_projected_point_3_object()(datum, query);
  }

  template <typename Primitive>
  static Point_3 closest_point_on_edges_vertices(const Point_3& query,
                                                 const Primitive& primitive,
                                                 const Point_3& bound,
                                                 const AABB_traits& traits)
  {
    const GT gt;
    const auto compare_distance = gt.compare_distance_3_object();
    const Point_3 candidate = Self::closest_point_on_edges_vertices_unbounded(query, primitive, traits);
    return (compare_distance(query, candidate, bound) == CGAL::SMALLER) ? candidate : bound;
  }

public:

  static bool first_intersection(const Point_3& p, const Point_3& q, Point_3& o,
                                 const FT offset_size,
                                 const FT intersection_precision,
                                 const AABBTree& tree)
  {
    CGAL_precondition(!tree.empty());

    using AABB_distance_oracle = internal::AABB_distance_oracle<AABBTree, AABBTraversalTraits>;
    using Offset_intersection = internal::Offset_intersection<GT, AABB_distance_oracle>;

    AABB_distance_oracle dist_oracle(tree);
    Offset_intersection offset_intersection(dist_oracle, offset_size, intersection_precision, 1 /*lip*/);
    return offset_intersection.first_intersection(p, q, o);
  }
};

template <typename GT,
          typename AABBTree,
          typename AABBTraversalTraits = CGAL::Default,
          typename BaseOracle = int> // base oracle
class AABB_tree_oracle
  : public BaseOracle
{
protected:
  using Geom_traits = GT;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;

  // When building oracle stacks, there are copies of (empty) trees, which isn't possible, thus pointer
  using AABB_tree = AABBTree;
  using AABB_tree_ptr = std::shared_ptr<AABB_tree>;
  using AABB_traits = typename AABB_tree::AABB_traits;
  using AABB_traversal_traits = typename Default::Get<AABBTraversalTraits,
                                                      Default_traversal_traits<AABB_traits> >::type;

  using AABB_helper = AABB_tree_oracle_helper<AABB_tree, AABB_traversal_traits>;

public:
  AABB_tree_oracle(const BaseOracle& base,
                   const GT& gt)
    : BaseOracle(base),
      m_gt(gt),
      m_tree_ptr(std::make_shared<AABB_tree>())
  { }

  AABB_tree_oracle(const AABB_tree_oracle&) = default;

public:
  const Geom_traits& geom_traits() const { return m_gt; }

  AABB_tree& tree() { return *m_tree_ptr; }
  const AABB_tree& tree() const { return *m_tree_ptr; }
  BaseOracle& base() { return static_cast<BaseOracle&>(*this); }
  const BaseOracle& base() const { return static_cast<const BaseOracle&>(*this); }

  bool empty() const { return m_tree_ptr->empty(); }
  bool do_call() const { return (!empty() || base().do_call()); }

  void clear() { m_tree_ptr->clear() && base().clear(); }

public:
  typename AABB_tree::Bounding_box bbox() const
  {
    CGAL_precondition(do_call());

    typename AABB_tree::Bounding_box bb;

    if(base().do_call())
      bb += base().bbox();

    if(!empty())
      bb += tree().bbox();

    return bb;
  }

  template <typename T>
  bool do_intersect(const T& t) const
  {
    CGAL_precondition(do_call());

    if(base().do_call() && base().do_intersect(t))
      return true;

    if(!empty())
      return AABB_helper::do_intersect(t, tree());

    return false;
  }

  FT squared_distance(const Point_3& p) const
  {
    CGAL_precondition(do_call());

    if(base().do_call())
    {
      if(!empty()) // both non empty
      {
        const FT base_sqd = base().squared_distance(p);
        // @speed could do a smarter traversal, no need to search deeper than the current best
        const FT this_sqd = AABB_helper::squared_distance(p, tree());
        return (std::min)(base_sqd, this_sqd);
      }
      else // this level is empty
      {
        return base().squared_distance(p);
      }
    }
    else // empty base
    {
      return AABB_helper::squared_distance(p, tree());
    }
  }

  Point_3 closest_point(const Point_3& p) const
  {
    CGAL_precondition(do_call());

    if(base().do_call())
    {
      if(!empty()) // both non empty
      {
        const Point_3 base_c = base().closest_point(p);
        // @speed could do a smarter traversal, no need to search deeper than the current best
        const Point_3 this_c = AABB_helper::closest_point(p, tree());
        return (compare_distance_to_point(p, base_c, this_c) == CGAL::SMALLER) ? base_c : this_c;
      }
      else // this level is empty
      {
        return base().closest_point(p);
      }
    }
    else // empty base
    {
      return AABB_helper::closest_point(p, tree());
    }
  }

  Point_3 closest_point_edge_vertex(const Point_3& p) const
  {
    CGAL_precondition(do_call());

    if(base().do_call())
    {
      if(!empty()) // both non empty
      {
        const Point_3 base_c = base().closest_point_edge_vertex(p);
        // @speed could do a smarter traversal, no need to search deeper than the current best
        const Point_3 this_c = AABB_helper::closest_point_edge_vertex(p, tree());
        return (compare_distance_to_point(p, base_c, this_c) == CGAL::SMALLER) ? base_c : this_c;
      }
      else // this level is empty
      {
        return base().closest_point_edge_vertex(p);
      }
    }
    else // empty base
    {
      return AABB_helper::closest_point_edge_vertex(p, tree());
    }
  }

  bool first_intersection(const Point_3& p, const Point_3& q,
                          Point_3& o,
                          const FT offset_size,
                          const FT intersection_precision) const
  {
    CGAL_precondition(do_call());

    if(base().do_call())
    {
      if(!empty()) // both non empty
      {
        Point_3 base_o;
        bool base_b = base().first_intersection(p, q, base_o, offset_size, intersection_precision);

        if(base_b) // intersection found in base
        {
          // @speed could do a smarter traversal, no need to search deeper than the current best
          Point_3 this_o;
          bool this_b = AABB_helper::first_intersection(p, q, this_o, offset_size, intersection_precision, tree());
          if(this_b)
            o = (compare_distance_to_point(p, base_o, this_o) == SMALLER) ? base_o : this_o;
          else
            o = base_o;

          return true;
        }
        else // no intersection found in non-empty base
        {
          return AABB_helper::first_intersection(p, q, o, offset_size, intersection_precision, tree());
        }
      }
      else // this level is empty
      {
        return base().first_intersection(p, q, o, offset_size, intersection_precision);
      }
    }
    else // empty base
    {
      return AABB_helper::first_intersection(p, q, o, offset_size, intersection_precision, tree());
    }
  }

  bool first_intersection(const Point_3& p, const Point_3& q,
                          Point_3& o,
                          const FT offset_size) const
  {
    return first_intersection(p, q, o, offset_size, 1e-2 * offset_size);
  }

protected:
  Geom_traits m_gt;
  AABB_tree_ptr m_tree_ptr;
};

// partial specialization, when there is no further oracle beneath in the stack.
//
// `int` is used to denote the absence of base rather than `void`,
// as to use the same constructor for both versions (thus requires a default construction)
template <typename GT,
          typename AABBTree,
          typename AABBTraversalTraits>
class AABB_tree_oracle<GT, AABBTree, AABBTraversalTraits, int>
{
protected:
  using Geom_traits = GT;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;

  using AABB_tree = AABBTree;
  using AABB_tree_ptr = std::shared_ptr<AABB_tree>;
  using AABB_traits = typename AABB_tree::AABB_traits;
  using AABB_traversal_traits = typename Default::Get<AABBTraversalTraits,
                                                      Default_traversal_traits<AABB_traits> >::type;

  using AABB_helper = AABB_tree_oracle_helper<AABB_tree, AABB_traversal_traits>;

public:
  AABB_tree_oracle(const int, // to have a common constructor API between both versions
                   const GT& gt)
    : m_gt(gt), m_tree_ptr(std::make_shared<AABB_tree>())
  { }

public:
  const Geom_traits& geom_traits() const { return m_gt; }
  AABB_tree& tree() { return *m_tree_ptr; }
  const AABB_tree& tree() const { return *m_tree_ptr; }

  bool empty() const { return m_tree_ptr->empty(); }
  bool do_call() const { return !empty(); }

  void clear() { m_tree_ptr->clear(); }

public:
  typename AABB_tree::Bounding_box bbox() const
  {
    CGAL_precondition(!empty());
    return tree().bbox();
  }

  template <typename T>
  bool do_intersect(const T& t) const
  {
    CGAL_precondition(!empty());
    return AABB_helper::do_intersect(t, tree());
  }

  FT squared_distance(const Point_3& p) const
  {
    CGAL_precondition(!empty());
    return AABB_helper::squared_distance(p, tree());
  }

  Point_3 closest_point(const Point_3& p) const
  {
    CGAL_precondition(!empty());
    return AABB_helper::closest_point(p, tree());
  }

  Point_3 closest_point_edge_vertex(const Point_3& p) const
  {
    CGAL_precondition(!empty());
    return AABB_helper::closest_point_edge_vertex(p, tree());
  }

  bool first_intersection(const Point_3& p, const Point_3& q, Point_3& o,
                          const FT offset_size, const FT intersection_precision) const
  {
    CGAL_precondition(!empty());
    return AABB_helper::first_intersection(p, q, o, offset_size, intersection_precision, tree());
  }

  bool first_intersection(const Point_3& p, const Point_3& q, Point_3& o, const FT offset_size) const
  {
    CGAL_precondition(!empty());
    return AABB_helper::first_intersection(p, q, o, offset_size, 1e-2 * offset_size, tree());
  }

private:
  Geom_traits m_gt;
  AABB_tree_ptr m_tree_ptr;
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_ORACLE_BASE_H
