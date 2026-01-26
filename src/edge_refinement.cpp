//
// Created by Michel Beeren on 26/01/2026.
//
#include <string>
#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <vector>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
using Vector_3 = K::Vector_3;
namespace PMP = CGAL::Polygon_mesh_processing;
using Mesh = CGAL::Surface_mesh<Point_3>;
using face_descriptor = Mesh::Face_index;
using Ray_3   = K::Ray_3;
using Segment_3 = K::Segment_3;
using Primitive   = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using AABB_traits = CGAL::AABB_traits<K, Primitive>;
using Tree        = CGAL::AABB_tree<AABB_traits>;

#include "edge_refinement.h"

// ------------------------------------REFINE ROUND EDGES-----------------------------------
Mesh offset_mesh(const Mesh& mesh, const double offset)
{
    const double d = 2.0 * offset;   // your factor

    using vd = Mesh::Vertex_index;

    Mesh out = mesh;                 // copy

    // vertex normals on the *output* mesh
    typedef Mesh::Property_map<vd, Vector_3> VNormalMap;
    VNormalMap vnormals;
    bool created;
    boost::tie(vnormals, created) =
        out.add_property_map<vd, Vector_3>("v:normal", CGAL::NULL_VECTOR);

    PMP::compute_vertex_normals(out, vnormals);

    for (vd v : out.vertices())
    {
        Point_3  p = out.point(v);
        Vector_3 n = vnormals[v];

        if (n.squared_length() == 0.0)
            continue;

        n = n / std::sqrt(n.squared_length()); // normalize
        out.point(v) = p + d * n;              // use -d*n for inward offset
    }

    return out;
}

Mesh offset_mesh_by_5cm(Mesh& mesh, const double offset_)
{
    const double d = 2*offset_; // 5 cm if units are meters

    // correct vertex descriptor type
    using vd = boost::graph_traits<Mesh>::vertex_descriptor;

    // new mesh
    Mesh offset_mesh = mesh;

    // 1) create / get a vertex normal property map
    typedef Mesh::Property_map<vd, Vector_3> VNormalMap;
    VNormalMap vnormals;
    bool created;
    boost::tie(vnormals, created) =
        offset_mesh.add_property_map<vd, Vector_3>("v:normal", CGAL::NULL_VECTOR);

    // 2) compute vertex normals
    PMP::compute_vertex_normals(offset_mesh, vnormals);

    // 3) move each vertex along its normal
    for (vd v : offset_mesh.vertices())
    {
        Point_3 p = offset_mesh.point(v);
        Vector_3  n = vnormals[v];

        if (n.squared_length() == 0.0)
            continue;

        n = n / std::sqrt(n.squared_length()); // normalize
        offset_mesh.point(v) = p + d * n;             // use -d * n to offset inward
    }
    return offset_mesh;
} // probably not necessary anymore

Mesh add_midpoint_distance_tag(const Mesh& wrapped_mesh, const Mesh& original_mesh, const double filter)
{
    using fd = Mesh::Face_index;
    using vd = Mesh::Vertex_index;
    using Color = CGAL::IO::Color;
    using hd = Mesh::Halfedge_index;

    Mesh mesh_to_tag = wrapped_mesh;

    // --- AABB tree on original mesh ---
    typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
    typedef CGAL::AABB_traits<K, Primitive>                AABB_traits;
    typedef CGAL::AABB_tree<AABB_traits>                   Tree;

    Tree tree(faces(original_mesh).begin(),
              faces(original_mesh).end(),
              original_mesh);
    tree.accelerate_distance_queries();

    // distance property
    typedef Mesh::Property_map<fd, double> FaceWeightMap;
    FaceWeightMap fweight;
    bool created;
    boost::tie(fweight, created) =
        mesh_to_tag.add_property_map<fd, double>("f:dist_to_original", 0.0);

    double dmin = std::numeric_limits<double>::max();
    double dmax = 0.0;

    // 1) compute distances + min/max
    for (fd f : mesh_to_tag.faces())
    {
        auto h  = mesh_to_tag.halfedge(f);
        vd v0   = mesh_to_tag.source(h);
        vd v1   = mesh_to_tag.target(h);
        vd v2   = mesh_to_tag.target(mesh_to_tag.next(h));

        const Point_3& p0 = mesh_to_tag.point(v0);
        const Point_3& p1 = mesh_to_tag.point(v1);
        const Point_3& p2 = mesh_to_tag.point(v2);

        Point_3 c( (p0.x()+p1.x()+p2.x())/3.0,
                   (p0.y()+p1.y()+p2.y())/3.0,
                   (p0.z()+p1.z()+p2.z())/3.0 );

        double dist = std::sqrt(tree.squared_distance(c));
        fweight[f] = dist;

        dmin = std::min(dmin, dist);
        dmax = std::max(dmax, dist);
    }

    // ------------------------------------create face color map--------------------------------------------
    Mesh::Property_map<fd, Color> fcolor;
    bool created_color;
    boost::tie(fcolor, created_color) =
        mesh_to_tag.add_property_map<fd, Color>("f:color", Color(255,255,255));

    // -----------------------create face bool map "f:refine" (default = false)-----------------------------
    Mesh::Property_map<fd, unsigned char> frefine;
    bool created_refine;
    boost::tie(frefine, created_refine) =
        mesh_to_tag.add_property_map<fd, unsigned char>("f:refine", 0);

    double range = (dmax > dmin) ? (dmax - dmin) : 1.0;

    // 3) map distance -> white..red
    for (fd f : mesh_to_tag.faces())
    {
        double d = fweight[f];
        double t = (d - dmin) / range;      // 0 = closest, 1 = farthest

        // white (255,255,255) -> red (255,0,0)
        unsigned char r = 255;
        unsigned char g = static_cast<unsigned char>(255 * (1.0 - t));
        unsigned char b = static_cast<unsigned char>(255 * (1.0 - t));
        fcolor[f] = Color(r, g, b);
        frefine[f] = (g < filter) ? 1 : 0;

    }

    // -------------------------Property if face is neighbor of too far face-------------------------

    // neighbor_refine flag: 1 if any neighbor has refine == 1
    Mesh::Property_map<fd, unsigned char> fneighbor_refine;
    bool created_neighbor;
    boost::tie(fneighbor_refine, created_neighbor) =
        mesh_to_tag.add_property_map<fd, unsigned char>("f:neighbor_refine", 0);

    // 2) compute neighbor_refine
    for (fd f : mesh_to_tag.faces())
    {
        bool has_refined_neighbor = false;

        // loop over halfedges around the face
        for (hd h : CGAL::halfedges_around_face(mesh_to_tag.halfedge(f), mesh_to_tag))
        {
            hd ho = mesh_to_tag.opposite(h);
            fd fn = mesh_to_tag.face(ho);

            if (fn != Mesh::null_face() && frefine[fn] == 1)
            {
                has_refined_neighbor = true;
                break;
            }
        }

        fneighbor_refine[f] = has_refined_neighbor ? 1 : 0;
    }



    return mesh_to_tag;
}

Mesh refine_round_edges(const Mesh& tagged_mesh, const Mesh& original_mesh)
{
    using fd    = Mesh::Face_index;
    using vd    = Mesh::Vertex_index;
    using Color = CGAL::IO::Color;

    Mesh refined;

    // ---- 0. refine flag on faces (tagged_mesh) ----
    auto opt_refine = tagged_mesh.property_map<fd, unsigned char>("f:refine");
    if (!opt_refine) {
        std::cerr << "Error: no face property 'f:refine' on input mesh\n";
        return refined;
    }
    auto frefine = *opt_refine;

    // ---- 0b. per-face color on tagged_mesh (optional) ----
    auto opt_fcolor_in = tagged_mesh.property_map<fd, Color>("f:color");
    bool has_color = static_cast<bool>(opt_fcolor_in);
    Mesh::Property_map<fd, Color> fcolor_in;
    if (has_color)
        fcolor_in = *opt_fcolor_in;

    // ---- 1. copy vertices to refined mesh ----
    std::vector<vd> vmap(num_vertices(tagged_mesh));
    for (vd v : tagged_mesh.vertices())
    {
        vd nv = refined.add_vertex(tagged_mesh.point(v));
        vmap[v.idx()] = nv;
    }

    // ---- 1b. create color map on refined mesh (if needed) ----
    Mesh::Property_map<fd, Color> fcolor_out;
    if (has_color) {
        bool created;
        boost::tie(fcolor_out, created) =
            refined.add_property_map<fd, Color>("f:color", Color(255,255,255));
    }

    // ---- 2. AABB tree on original_mesh for ray intersections ----
    // using Primitive   = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
    // using AABB_traits = CGAL::AABB_traits<K, Primitive>;
    // using Tree        = CGAL::AABB_tree<AABB_traits>;
    using edge_descriptor     = boost::graph_traits<Mesh>::edge_descriptor;
    using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;

    Tree tree(faces(original_mesh).begin(),
              faces(original_mesh).end(),
              original_mesh);
    tree.accelerate_distance_queries();

    // ---- helpers: normals computed on TAGGED mesh ----

    // face normal on tagged_mesh
    auto face_normal = [&](fd f) -> Vector_3 {
        halfedge_descriptor h = tagged_mesh.halfedge(f);
        vd v0   = tagged_mesh.source(h);
        vd v1   = tagged_mesh.target(h);
        vd v2   = tagged_mesh.target(tagged_mesh.next(h));

        const Point_3& p0 = tagged_mesh.point(v0);
        const Point_3& p1 = tagged_mesh.point(v1);
        const Point_3& p2 = tagged_mesh.point(v2);

        Vector_3 u = p1 - p0;
        Vector_3 v = p2 - p0;
        return CGAL::cross_product(u, v);  // not normalized
    };

    // average normal of the two faces incident to edge (a,b) in TAGGED mesh
    auto edge_avg_normal = [&](vd a, vd b) -> Vector_3 {
        edge_descriptor e;
        bool found = false;

        boost::tie(e, found) = edge(a, b, tagged_mesh);
        if (!found)
            boost::tie(e, found) = edge(b, a, tagged_mesh);
        if (!found)
            return Vector_3(0,0,0); // no such edge in tagged_mesh

        halfedge_descriptor h = halfedge(e, tagged_mesh);

        Vector_3 n(0,0,0);

        fd f1 = face(h, tagged_mesh);
        if (f1 != Mesh::null_face())
            n = n + face_normal(f1);

        fd f2 = face(opposite(h, tagged_mesh), tagged_mesh);
        if (f2 != Mesh::null_face() && f2 != f1)
            n = n + face_normal(f2);

        if (n.squared_length() > 0.0) {
            double len = std::sqrt(n.squared_length());
            n = n / len;
        }
        return n;
    };

    // ---- 3. mark edges to split (those touching a refined face) ----
    using EdgeKey = std::pair<std::size_t, std::size_t>;

    auto edge_key = [](vd a, vd b) -> EdgeKey {
        std::size_t ia = a.idx();
        std::size_t ib = b.idx();
        if (ia > ib) std::swap(ia, ib);
        return EdgeKey(ia, ib);
    };

    std::set<EdgeKey> edges_to_split;

    for (fd f : tagged_mesh.faces())
    {
        if (frefine[f] == 0) continue;

        auto h = tagged_mesh.halfedge(f);
        vd v0  = tagged_mesh.source(h);
        vd v1  = tagged_mesh.target(h);
        vd v2  = tagged_mesh.target(tagged_mesh.next(h));

        edges_to_split.insert(edge_key(v0, v1));
        edges_to_split.insert(edge_key(v1, v2));
        edges_to_split.insert(edge_key(v2, v0));
    }

    // ---- 4. split edges, create midpoint vertices (projected onto original_mesh) ----
    std::map<EdgeKey, vd> edge_midpoint;

    auto get_midpoint_vertex = [&](vd a, vd b) -> vd
    {
        EdgeKey key = edge_key(a, b);
        auto it = edge_midpoint.find(key);
        if (it != edge_midpoint.end())
            return it->second;

        // geometric midpoint from tagged_mesh
        const Point_3& p0 = tagged_mesh.point(a);
        const Point_3& p1 = tagged_mesh.point(b);

        Point_3 mid(
            (p0.x() + p1.x()) / 2.0,
            (p0.y() + p1.y()) / 2.0,
            (p0.z() + p1.z()) / 2.0
        );

        Point_3 final_pos = mid;

        // average normal along edge in TAGGED mesh
        Vector_3 n = edge_avg_normal(a, b);
        if (n.squared_length() > 0.0)
        {
            Vector_3 nn = n / std::sqrt(n.squared_length());

            // move in *opposite* normal direction
            Ray_3 ray(mid, mid - nn);   // flip to mid + nn if needed

            auto opt = tree.first_intersection(ray);
            if (opt) {
                // opt->first is a CGAL::Object
                const CGAL::Object& obj = opt->first;
                Point_3 ip;
                if (CGAL::assign(ip, obj)) {
                    final_pos = ip;
                }
            }
        }

        vd mv = refined.add_vertex(final_pos);
        edge_midpoint[key] = mv;
        return mv;
    };

    // ---- 5. build faces in refined mesh (0/1/2/3 split edges) ----
    for (fd f : tagged_mesh.faces())
    {
        auto h = tagged_mesh.halfedge(f);
        vd v0  = tagged_mesh.source(h);
        vd v1  = tagged_mesh.target(h);
        vd v2  = tagged_mesh.target(tagged_mesh.next(h));

        vd A = vmap[v0.idx()];
        vd B = vmap[v1.idx()];
        vd C = vmap[v2.idx()];

        bool sAB = edges_to_split.count(edge_key(v0, v1)) != 0;
        bool sBC = edges_to_split.count(edge_key(v1, v2)) != 0;
        bool sCA = edges_to_split.count(edge_key(v2, v0)) != 0;

        int n_split = int(sAB) + int(sBC) + int(sCA);

        // color of this original face
        Color c_face = has_color ? fcolor_in[f] : Color(255,255,255);

        auto add_face = [&](vd a, vd b, vd c) -> fd {
            fd nf = refined.add_face(a, b, c);
            if (has_color && nf != Mesh::null_face()) {
                fcolor_out[nf] = c_face;
            }
            return nf;
        };

        if (n_split == 0)
        {
            add_face(A, B, C);
        }
        else if (n_split == 3)
        {
            vd MAB = get_midpoint_vertex(v0, v1);
            vd MBC = get_midpoint_vertex(v1, v2);
            vd MCA = get_midpoint_vertex(v2, v0);

            add_face(A,   MAB, MCA);
            add_face(B,   MBC, MAB);
            add_face(C,   MCA, MBC);
            add_face(MAB, MBC, MCA);
        }
        else if (n_split == 1)
        {
            if (sAB)
            {
                vd MAB = get_midpoint_vertex(v0, v1);
                add_face(A,   MAB, C);
                add_face(MAB, B,   C);
            }
            else if (sBC)
            {
                vd MBC = get_midpoint_vertex(v1, v2);
                add_face(B,   MBC, A);
                add_face(MBC, C,   A);
            }
            else // sCA
            {
                vd MCA = get_midpoint_vertex(v2, v0);
                add_face(C,   MCA, B);
                add_face(MCA, A,   B);
            }
        }
        else // n_split == 2
        {
            if (sAB && sBC && !sCA)
            {
                vd MAB = get_midpoint_vertex(v0, v1);
                vd MBC = get_midpoint_vertex(v1, v2);

                add_face(A,   MAB, C);
                add_face(MAB, MBC, C);
                add_face(MAB, B,   MBC);
            }
            else if (sBC && sCA && !sAB)
            {
                vd MBC = get_midpoint_vertex(v1, v2);
                vd MCA = get_midpoint_vertex(v2, v0);

                add_face(B,   MBC, A);
                add_face(MBC, MCA, A);
                add_face(MBC, C,   MCA);
            }
            else if (sCA && sAB && !sBC)
            {
                vd MCA = get_midpoint_vertex(v2, v0);
                vd MAB = get_midpoint_vertex(v0, v1);

                add_face(C,   MCA, B);
                add_face(MCA, MAB, B);
                add_face(MCA, A,   MAB);
            }
        }
    }

    return refined;
}

Mesh refine_round_edges3(const Mesh& tagged_mesh, const Mesh& original_mesh)
{
    using fd    = Mesh::Face_index;
    using vd    = Mesh::Vertex_index;
    using Color = CGAL::IO::Color;

    Mesh refined;

    // ---- 0. refine flag on faces (tagged_mesh) ----
    auto opt_refine = tagged_mesh.property_map<fd, unsigned char>("f:refine");
    if (!opt_refine) {
        std::cerr << "Error: no face property 'f:refine' on input mesh\n";
        return refined;
    }
    auto frefine = *opt_refine;

    // ---- 0b. per-face color on tagged_mesh (optional) ----
    auto opt_fcolor_in = tagged_mesh.property_map<fd, Color>("f:color");
    bool has_color = static_cast<bool>(opt_fcolor_in);
    Mesh::Property_map<fd, Color> fcolor_in;
    if (has_color)
        fcolor_in = *opt_fcolor_in;

    // ---- 1. copy vertices to refined mesh ----
    std::vector<vd> vmap(num_vertices(tagged_mesh));
    for (vd v : tagged_mesh.vertices())
    {
        vd nv = refined.add_vertex(tagged_mesh.point(v));
        vmap[v.idx()] = nv;
    }

    // ---- 1b. create color map on refined mesh (if needed) ----
    Mesh::Property_map<fd, Color> fcolor_out;
    if (has_color) {
        bool created;
        boost::tie(fcolor_out, created) =
            refined.add_property_map<fd, Color>("f:color", Color(255,255,255));
    }

    // ---- 2. AABB tree on original_mesh for ray intersections ----
    using Primitive   = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
    using AABB_traits = CGAL::AABB_traits<K, Primitive>;
    using Tree        = CGAL::AABB_tree<AABB_traits>;
    using edge_descriptor     = boost::graph_traits<Mesh>::edge_descriptor;
    using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;

    Tree tree(faces(original_mesh).begin(),
              faces(original_mesh).end(),
              original_mesh);
    tree.accelerate_distance_queries();

    // helper: face normal on original_mesh
    auto face_normal = [&](fd f) -> Vector_3 {
        halfedge_descriptor h = original_mesh.halfedge(f);
        vd v0   = original_mesh.source(h);
        vd v1   = original_mesh.target(h);
        vd v2   = original_mesh.target(original_mesh.next(h));

        const Point_3& p0 = original_mesh.point(v0);
        const Point_3& p1 = original_mesh.point(v1);
        const Point_3& p2 = original_mesh.point(v2);

        Vector_3 u = p1 - p0;
        Vector_3 v = p2 - p0;
        return CGAL::cross_product(u, v);  // not normalized
    };

    // helper: average normal of the two faces incident to edge (a,b)
    auto edge_avg_normal = [&](vd a, vd b) -> Vector_3 {
        edge_descriptor e;
        bool found = false;

        boost::tie(e, found) = edge(a, b, original_mesh);
        if (!found)
            boost::tie(e, found) = edge(b, a, original_mesh);
        if (!found)
            return Vector_3(0,0,0); // edge not in original_mesh

        halfedge_descriptor h = halfedge(e, original_mesh);

        Vector_3 n(0,0,0);

        fd f1 = face(h, original_mesh);
        if (f1 != Mesh::null_face())
            n = n + face_normal(f1);

        fd f2 = face(opposite(h, original_mesh), original_mesh);
        if (f2 != Mesh::null_face() && f2 != f1)
            n = n + face_normal(f2);

        if (n.squared_length() > 0.0) {
            double len = std::sqrt(n.squared_length());
            n = n / len;
        }
        return n;
    };

    // ---- 3. mark edges to split (those touching a refined face) ----
    using EdgeKey = std::pair<std::size_t, std::size_t>;

    auto edge_key = [](vd a, vd b) -> EdgeKey {
        std::size_t ia = a.idx();
        std::size_t ib = b.idx();
        if (ia > ib) std::swap(ia, ib);
        return EdgeKey(ia, ib);
    };

    std::set<EdgeKey> edges_to_split;

    for (fd f : tagged_mesh.faces())
    {
        if (frefine[f] == 0) continue;

        auto h = tagged_mesh.halfedge(f);
        vd v0  = tagged_mesh.source(h);
        vd v1  = tagged_mesh.target(h);
        vd v2  = tagged_mesh.target(tagged_mesh.next(h));

        edges_to_split.insert(edge_key(v0, v1));
        edges_to_split.insert(edge_key(v1, v2));
        edges_to_split.insert(edge_key(v2, v0));
    }

    // ---- 4. split edges, create midpoint vertices (projected to original) ----
    std::map<EdgeKey, vd> edge_midpoint;

    // Return type of first_intersection for a Ray_3
    using Ray_intersection =
        Tree::Intersection_and_primitive_id<Ray_3>::Type;

    auto get_midpoint_vertex = [&](vd a, vd b) -> vd
    {
        EdgeKey key = edge_key(a, b);
        auto it = edge_midpoint.find(key);
        if (it != edge_midpoint.end())
            return it->second;

        // geometric midpoint from tagged_mesh
        const Point_3& p0 = tagged_mesh.point(a);
        const Point_3& p1 = tagged_mesh.point(b);

        Point_3 mid( (p0.x() + p1.x()) / 2.0,
                     (p0.y() + p1.y()) / 2.0,
                     (p0.z() + p1.z()) / 2.0 );

        Point_3 final_pos = mid;

        // average normal along edge in original_mesh
        Vector_3 n = edge_avg_normal(a, b);
        if (n.squared_length() > 0.0)
        {
            Vector_3 nn = n / std::sqrt(n.squared_length());

            // move in *opposite* normal direction
            Ray_3 ray(mid, mid - nn);   // Point, Point (mid - nn is a Point_3)

            // let the compiler deduce the exact optional type
            auto opt = tree.first_intersection(ray);

            if (opt) {

                // in current CGAL, *opt is a pair<CGAL::Object, Primitive_id>
                const CGAL::Object& obj = opt->first;

                Point_3 ip;
                if (CGAL::assign(ip, obj)) {
                    final_pos = ip;
                }
            }
        }

        vd mv = refined.add_vertex(final_pos);
        edge_midpoint[key] = mv;
        return mv;
    };

    // ---- 5. build faces in refined mesh (0/1/2/3 split edges) ----
    for (fd f : tagged_mesh.faces())
    {
        auto h = tagged_mesh.halfedge(f);
        vd v0  = tagged_mesh.source(h);
        vd v1  = tagged_mesh.target(h);
        vd v2  = tagged_mesh.target(tagged_mesh.next(h));

        vd A = vmap[v0.idx()];
        vd B = vmap[v1.idx()];
        vd C = vmap[v2.idx()];

        bool sAB = edges_to_split.count(edge_key(v0, v1)) != 0;
        bool sBC = edges_to_split.count(edge_key(v1, v2)) != 0;
        bool sCA = edges_to_split.count(edge_key(v2, v0)) != 0;

        int n_split = int(sAB) + int(sBC) + int(sCA);

        // color of this original face
        Color c_face = has_color ? fcolor_in[f] : Color(255,255,255);

        auto add_face = [&](vd a, vd b, vd c) -> fd {
            fd nf = refined.add_face(a, b, c);
            if (has_color && nf != Mesh::null_face()) {
                fcolor_out[nf] = c_face;
            }
            return nf;
        };

        if (n_split == 0)
        {
            add_face(A, B, C);
        }
        else if (n_split == 3)
        {
            vd MAB = get_midpoint_vertex(v0, v1);
            vd MBC = get_midpoint_vertex(v1, v2);
            vd MCA = get_midpoint_vertex(v2, v0);

            add_face(A,   MAB, MCA);
            add_face(B,   MBC, MAB);
            add_face(C,   MCA, MBC);
            add_face(MAB, MBC, MCA);
        }
        else if (n_split == 1)
        {
            if (sAB)
            {
                vd MAB = get_midpoint_vertex(v0, v1);
                add_face(A,   MAB, C);
                add_face(MAB, B,   C);
            }
            else if (sBC)
            {
                vd MBC = get_midpoint_vertex(v1, v2);
                add_face(B,   MBC, A);
                add_face(MBC, C,   A);
            }
            else // sCA
            {
                vd MCA = get_midpoint_vertex(v2, v0);
                add_face(C,   MCA, B);
                add_face(MCA, A,   B);
            }
        }
        else // n_split == 2
        {
            if (sAB && sBC && !sCA)
            {
                vd MAB = get_midpoint_vertex(v0, v1);
                vd MBC = get_midpoint_vertex(v1, v2);

                add_face(A,   MAB, C);
                add_face(MAB, MBC, C);
                add_face(MAB, B,   MBC);
            }
            else if (sBC && sCA && !sAB)
            {
                vd MBC = get_midpoint_vertex(v1, v2);
                vd MCA = get_midpoint_vertex(v2, v0);

                add_face(B,   MBC, A);
                add_face(MBC, MCA, A);
                add_face(MBC, C,   MCA);
            }
            else if (sCA && sAB && !sBC)
            {
                vd MCA = get_midpoint_vertex(v2, v0);   // <-- fixed: v0, not 0
                vd MAB = get_midpoint_vertex(v0, v1);

                add_face(C,   MCA, B);
                add_face(MCA, MAB, B);
                add_face(MCA, A,   MAB);
            }
        }
    }

    return refined;
}

Mesh refine_round_edges2(const Mesh& tagged_mesh)
{
    using fd    = Mesh::Face_index;
    using vd    = Mesh::Vertex_index;
    using Color = CGAL::IO::Color;

    Mesh refined;

    // --- 0. Get refine-flag per face on the original mesh ---
    auto opt_refine = tagged_mesh.property_map<fd, unsigned char>("f:refine");
    if (!opt_refine) {
        std::cerr << "Error: no face property 'f:refine' on input mesh\n";
        return refined;
    }
    auto frefine = *opt_refine;

    // --- 0b. Get per-face color map on the original mesh (if present) ---
    auto opt_fcolor_in = tagged_mesh.property_map<fd, Color>("f:color");
    bool has_color = static_cast<bool>(opt_fcolor_in);
    Mesh::Property_map<fd, Color> fcolor_in;
    if (has_color) {
        fcolor_in = *opt_fcolor_in;
    }

    // --- 1. Copy all original vertices to the new mesh ---
    std::vector<vd> vmap(num_vertices(tagged_mesh));
    for (vd v : tagged_mesh.vertices())
    {
        vd nv = refined.add_vertex(tagged_mesh.point(v));
        vmap[v.idx()] = nv;
    }

    // --- 1b. Create per-face color map on the refined mesh (if needed) ---
    Mesh::Property_map<fd, Color> fcolor_out;
    if (has_color) {
        bool created_color_out;
        boost::tie(fcolor_out, created_color_out) =
            refined.add_property_map<fd, Color>("f:color", Color(255,255,255));
    }

    // --- 2. Decide which edges must be split (any edge touching a refined face) ---
    using EdgeKey = std::pair<std::size_t, std::size_t>;

    auto edge_key = [](vd a, vd b) -> EdgeKey {
        std::size_t ia = a.idx();
        std::size_t ib = b.idx();
        if (ia > ib) std::swap(ia, ib);
        return EdgeKey(ia, ib);
    };

    std::set<EdgeKey> edges_to_split;

    for (fd f : tagged_mesh.faces())
    {
        if (frefine[f] == 0) continue;

        auto h = tagged_mesh.halfedge(f);
        vd v0  = tagged_mesh.source(h);
        vd v1  = tagged_mesh.target(h);
        vd v2  = tagged_mesh.target(tagged_mesh.next(h));

        edges_to_split.insert(edge_key(v0, v1));
        edges_to_split.insert(edge_key(v1, v2));
        edges_to_split.insert(edge_key(v2, v0));
    }

    // --- 3. Map each split edge to a single midpoint vertex in the new mesh ---
    std::map<EdgeKey, vd> edge_midpoint;

    auto get_midpoint_vertex = [&](vd a, vd b) -> vd
    {
        EdgeKey key = edge_key(a, b);
        auto it = edge_midpoint.find(key);
        if (it != edge_midpoint.end())
            return it->second;

        const Point_3& p0 = tagged_mesh.point(a);
        const Point_3& p1 = tagged_mesh.point(b);

        Point_3 mid( (p0.x() + p1.x()) / 2.0,
                     (p0.y() + p1.y()) / 2.0,
                     (p0.z() + p1.z()) / 2.0 );

        vd mv = refined.add_vertex(mid);
        edge_midpoint[key] = mv;
        return mv;
    };

    // --- 4. Build faces in the refined mesh, case by case ---
    for (fd f : tagged_mesh.faces())
    {
        auto h = tagged_mesh.halfedge(f);
        vd v0  = tagged_mesh.source(h);
        vd v1  = tagged_mesh.target(h);
        vd v2  = tagged_mesh.target(tagged_mesh.next(h));

        vd A = vmap[v0.idx()];
        vd B = vmap[v1.idx()];
        vd C = vmap[v2.idx()];

        bool sAB = edges_to_split.count(edge_key(v0, v1)) != 0;
        bool sBC = edges_to_split.count(edge_key(v1, v2)) != 0;
        bool sCA = edges_to_split.count(edge_key(v2, v0)) != 0;

        int n_split = int(sAB) + int(sBC) + int(sCA);

        // color of this original face
        Color c_face = has_color ? fcolor_in[f] : Color(255,255,255);

        // helper: add face & copy color
        auto add_face = [&](vd a, vd b, vd c) -> fd {
            fd nf = refined.add_face(a, b, c);
            if (has_color && nf != Mesh::null_face()) {
                fcolor_out[nf] = c_face;
            }
            return nf;
        };

        if (n_split == 0)
        {
            // no refinement on this face
            add_face(A, B, C);
        }
        else if (n_split == 3)
        {
            // all three edges split -> 4 triangles
            vd MAB = get_midpoint_vertex(v0, v1);
            vd MBC = get_midpoint_vertex(v1, v2);
            vd MCA = get_midpoint_vertex(v2, v0);

            add_face(A,   MAB, MCA); // 1 4 6
            add_face(B,   MBC, MAB); // 2 5 4
            add_face(C,   MCA, MBC); // 3 6 5
            add_face(MAB, MBC, MCA); // 4 5 6
        }
        else if (n_split == 1)
        {
            // exactly one split edge -> 2 triangles
            if (sAB)
            {
                vd MAB = get_midpoint_vertex(v0, v1);
                add_face(A,   MAB, C); // 1 4 3
                add_face(MAB, B,   C); // 4 2 3
            }
            else if (sBC)
            {
                vd MBC = get_midpoint_vertex(v1, v2);
                add_face(B,   MBC, A); // 2 5 1
                add_face(MBC, C,   A); // 5 3 1
            }
            else // sCA
            {
                vd MCA = get_midpoint_vertex(v2, v0);
                add_face(C,   MCA, B); // 3 6 2
                add_face(MCA, A,   B); // 6 1 2
            }
        }
        else // n_split == 2
        {
            // two split edges -> 3 triangles
            if (sAB && sBC && !sCA)
            {
                // edges around B: AB, BC
                vd MAB = get_midpoint_vertex(v0, v1);
                vd MBC = get_midpoint_vertex(v1, v2);

                add_face(A,   MAB, C);    // 1 4 3
                add_face(MAB, MBC, C);    // 4 5 3
                add_face(MAB, B,   MBC);  // 4 2 5
            }
            else if (sBC && sCA && !sAB)
            {
                // edges around C: BC, CA
                vd MBC = get_midpoint_vertex(v1, v2);
                vd MCA = get_midpoint_vertex(v2, v0);

                add_face(B,   MBC, A);    // 2 5 1
                add_face(MBC, MCA, A);    // 5 6 1
                add_face(MBC, C,   MCA);  // 5 3 6
            }
            else if (sCA && sAB && !sBC)
            {
                // edges around A: CA, AB
                vd MCA = get_midpoint_vertex(v2, v0);
                vd MAB = get_midpoint_vertex(v0, v1);

                add_face(C,   MCA, B);    // 3 6 2
                add_face(MCA, MAB, B);    // 6 4 2
                add_face(MCA, A,   MAB);  // 6 1 4
            }
        }
    }

    return refined;
}