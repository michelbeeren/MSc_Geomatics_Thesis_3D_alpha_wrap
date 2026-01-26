#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <map>
#include <set>

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Object.h>
#include <CGAL/boost/graph/helpers.h>

#include <boost/optional.hpp>
#include <boost/variant/get.hpp>

#include <iostream>
#include <string>
#include <filesystem>

#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Random.h>
#include <CGAL/IO/PLY.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <limits>
#include <CGAL/IO/Color.h>
#include <limits>
#include <CGAL/IO/Color.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>

#include <boost/variant/get.hpp>

// random

namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;

using Mesh = CGAL::Surface_mesh<Point_3>;
using face_descriptor = Mesh::Face_index;
using Ray_3   = K::Ray_3;
using Segment_3 = K::Segment_3;
using Primitive   = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using AABB_traits = CGAL::AABB_traits<K, Primitive>;
using Tree        = CGAL::AABB_tree<AABB_traits>;

#include <fstream>
#include <unordered_map>
#include <vector>
#include <iomanip>

// ------------------------------ VAL3DITY ------------------------------------------
// write CityJSON file from mesh to validate it using val3dity
void write_cityjson_from_mesh(const Mesh& mesh, const std::string& filename)
{
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error("Cannot open output CityJSON file");

    // Set high precision to prevent snapping vertices
    out << std::fixed << std::setprecision(17);

    using VI = Mesh::Vertex_index;
    using FI = Mesh::Face_index;

    // Building index
    std::vector<VI> index_to_vertex;
    index_to_vertex.reserve(num_vertices(mesh));
    std::unordered_map<VI,int> vidx;

    int idx = 0;
    for (VI v : mesh.vertices()) {
        index_to_vertex.push_back(v);
        vidx[v] = idx++;
    }

    // Collect vertices
    std::vector<Point_3> verts;
    verts.reserve(index_to_vertex.size());
    for (VI v : index_to_vertex)
        verts.push_back(mesh.point(v));

    // Collect faces
    std::vector<std::vector<int>> faces;
    faces.reserve(num_faces(mesh));

    for (FI f : mesh.faces()) {
        std::vector<int> face;
        for (VI v : CGAL::vertices_around_face(mesh.halfedge(f), mesh))
            face.push_back(vidx[v]);

        if (face.size() < 3) continue;
        faces.push_back(std::move(face));
    }

    // Write CityJSON for 1 building (1 solid)
    out << "{\n";
    out << "  \"type\": \"CityJSON\",\n";
    out << "  \"version\": \"1.1\",\n";
    out << "  \"CityObjects\": {\n";
    out << "    \"wrap\": {\n";
    out << "      \"type\": \"Building\",\n";
    out << "      \"geometry\": [\n";
    out << "        {\n";
    out << "          \"type\": \"Solid\",\n";
    out << "          \"lod\": 2,\n";
    out << "          \"boundaries\": [\n";
    out << "            [\n"; // Start Shell


    for (std::size_t i = 0; i < faces.size(); ++i) {
        // Each surface is a new polygon which starts with an exterior ring; structure [ [ v1, v2, v3 ] ]
        out << "              [ [";
        for (std::size_t j = 0; j < faces[i].size(); ++j) {
            out << faces[i][j];
            if (j + 1 < faces[i].size()) out << ", ";
        }
        out << "] ]"; // Close ring; close polygon;

        if (i + 1 < faces.size()) out << ",";
        out << "\n";
    }

    out << "            ]\n";     // End Shell
    out << "          ]\n";       // End boundaries
    out << "        }\n";        // End geometry object
    out << "      ]\n";          // End geometry array
    out << "    }\n";            // End CityObject wrap
    out << "  },\n";             // End CityObjects
    out << "  \"vertices\": [\n";

    for (std::size_t k = 0; k < verts.size(); ++k) {
        out << "    ["
            << verts[k].x() << ", "
            << verts[k].y() << ", "
            << verts[k].z() << "]";
        if (k + 1 < verts.size()) out << ",";
        out << "\n";
    }

    out << "  ]\n";              // vertices
    out << "}\n";                // root
}

// run val3dity and get error code back
int run_val3dity(const std::string& input)
{
    std::string exe = "/opt/homebrew/bin/val3dity";
    // default snap tolerance = 0.001, Now set lower for really small edges (in case of small alpha)
    std::string cmd = exe + " --snap_tol 1e-06 " + input;
    return std::system(cmd.c_str());
} // probably not necessary anymore

// if you do not want to print the report
int run_val3dity_silent(const std::string& input)
{
    std::string exe = "/opt/homebrew/bin/val3dity";

    // Add " > /dev/null 2>&1" at the end.
    // > /dev/null    : sends standard output to the void.
    // 2>&1           : sends standard error errors to the same place as standard output.
    std::string cmd = exe + " --snap_tol 1e-06 " + input + " > /dev/null 2>&1";

    return std::system(cmd.c_str());
} // probably not necessary anymore

// check val3dity report if mesh is valid
bool check_report_for_validity(const std::string& report_path) {
    std::ifstream file(report_path);
    if (!file.is_open()) return false;

    std::string content((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());

    // val3dity reports typically contain an "errors" array.
    // If the mesh is valid, the errors list is empty or the summary says valid.
    // NOTE: Adjust this logic based on the exact version of val3dity you use.
    // Newer versions (v2.6+) generally produce a structure where you look for empty error arrays.

    // A robust check (pseudo-logic):
    // if content contains "errors": [] -> Valid
    // if content contains "errors": [ ... ] -> Invalid

    // Simple heuristic: If the file mentions "validity": true or has no error codes.
    // For safety, assume invalid unless proven valid.

    // You might want to print the report to debug once to see the JSON structure.
    return content.find("\"errors\": []") != std::string::npos;
}

// returns and prints if input mesh is valid
bool valid_mesh_boolean(const Mesh& mesh_)
{
    std::string input_path = "../data/Test_3D_Alphawrap/Output/val3dity_check/check_me.json";
    std::string report_path = "../data/Test_3D_Alphawrap/Output/val3dity_check/report.json";

    write_cityjson_from_mesh(mesh_, input_path);

    // 1. Run val3dity silently, but force it to write a report file
    std::string exe = "/opt/homebrew/bin/val3dity";

    // We add "--report" to save the results to a file
    // We keep "> /dev/null" so it doesn't clutter your console
    std::string cmd = exe + " --snap_tol 1e-06 --report " + report_path + " " + input_path + " > /dev/null 2>&1";

    int system_code = std::system(cmd.c_str());

    // 2. Check if the PROGRAM ran successfully (system_code == 0)
    // AND if the REPORT says the mesh is valid.
    bool program_ran_ok = (system_code == 0);
    bool mesh_is_valid = false;

    if (program_ran_ok) {
        mesh_is_valid = check_report_for_validity(report_path);
    }

    if (mesh_is_valid) {
        std::cout << "Mesh is VALID :)" << std::endl;
        return true;
    } else {
        std::cout << "Mesh is INVALID :(" << std::endl;
        // Optional: Print the report location so you can inspect it manually if needed
        // std::cout << "See report at: " << report_path << std::endl;
        return false;
    }
}

// --------------------------------------ALPHA WRAP-------------------------------------
// generate output name
std::string generate_output_name(const std::string input_name_, const double rel_alpha_, const double rel_offset_)
{
    const int rel_alpha_for_txt = rel_alpha_;
    const int rel_offset_for_txt = rel_offset_;
    const std::string rel_alpha_txt = std::to_string(rel_alpha_for_txt);
    const std::string rel_offset_txt = std::to_string(rel_offset_for_txt);

    // Replace "Input" → "Output"
    std::string out_path = input_name_;
    std::string from = "Input";
    std::string to   = "Output";

    size_t pos = out_path.find(from);
    if (pos != std::string::npos) {
        out_path.replace(pos, from.length(), to);
    }

    // Remove last 4 characters (".off")
    if (out_path.size() > 4 && out_path.substr(out_path.size() - 4) == ".off") {
        out_path = out_path.substr(0, out_path.size() - 4);
    }

    // Add cuurent alpha and offset to the output name
    std::string out_name = out_path + "_a=" + rel_alpha_txt + "_offset=" + rel_offset_txt + ".off";
    return out_name;
}

// 3D alpha wrap
Mesh _3D_alpha_wrap(const std::string output_, const double relative_alpha_, const double relative_offset_, Mesh& mesh_) {
    // std::cout << "omg omg it is a triangle mesh!;" << std::endl;
    // std::cout << "relative_alpha: " << relative_alpha_ << " and relative offset: " << relative_offset_ << std::endl;
    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh_);
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                         CGAL::square(bbox.ymax() - bbox.ymin()) +
                                         CGAL::square(bbox.zmax() - bbox.zmin()));
    // std::cout << "diagonal bbox length: " << diag_length << std::endl;
    const double alpha = diag_length / relative_alpha_;
    const double offset = diag_length / relative_offset_;
    std::cout << "Alpha wrapping with a_rel = " << relative_alpha_ << " and offset_rel = " << relative_offset_ << ", which results in alpha = " << alpha << " and offset = " << offset << std::endl;
    // Construct the wrap
    CGAL::Real_timer t;
    t.start();

    Mesh wrap;
    CGAL::alpha_wrap_3(mesh_, alpha, offset, wrap);

    t.stop();
    std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces, it took " << t.time() << " s." << std::endl;
    // std::cout << "Took " << t.time() << " s." << std::endl;

    // ------------------------------ write output mesh -------------------------------
    const std::string output_name = output_;
    std::cout << "Writing to " << output_name << std::endl;
    CGAL::IO::write_polygon_mesh(output_name, wrap, CGAL::parameters::stream_precision(17));

    return wrap;

}

// compute offset
double rel_offset_to_offset(Mesh& mesh, const double relative_offset)
{
    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                         CGAL::square(bbox.ymax() - bbox.ymin()) +
                                         CGAL::square(bbox.zmax() - bbox.zmin()));
    // std::cout << "diagonal bbox length: " << diag_length << std::endl;
    const double offset = diag_length / relative_offset;
    return offset;
}

// -------------------------------------WRAP FROM INSIDE------------------------------------
// TODO not sure if this is really reliable point placer if input mesh is not valid
Point_3 random_point_inside_mesh(const Mesh& mesh)
{
    // containment oracle
    CGAL::Side_of_triangle_mesh<Mesh, K> inside(mesh);

    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
    CGAL::Random rng;

    while (true)
    {
        double x = rng.get_double(bbox.xmin(), bbox.xmax());
        double y = rng.get_double(bbox.ymin(), bbox.ymax());
        double z = rng.get_double(bbox.zmin(), bbox.zmax());
        Point_3 p(x, y, z);

        if (inside(p) == CGAL::ON_BOUNDED_SIDE)
            return p;
    }
}

// 3D alpha wrap from inside TODO point inside generation can be better
Mesh _3D_alpha_inside_wrap(const std::string output_, const double relative_alpha_, const double relative_offset_, Mesh& mesh_) {
  // std::cout << "omg omg it is a triangle mesh!;" << std::endl;
  // std::cout << "relative_alpha: " << relative_alpha_ << " and relative offset: " << relative_offset_ << std::endl;
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh_);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));
  // std::cout << "diagonal bbox length: " << diag_length << std::endl;
  const double alpha = diag_length / relative_alpha_;
  const double offset = diag_length / relative_offset_;
    std::cout << "Alpha wrapping from inside with a_rel = " << relative_alpha_ << " and offset_rel = " << relative_offset_ << ", which results in alpha = " << alpha << " and offset = " << offset << std::endl;

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

    // // point now placed in middle of bounding box
    // std::vector<Point_3> seeds = {
    //     Point_3((bbox.xmin() + bbox.xmax()) / 2,(bbox.ymin() + bbox.ymax()) / 2,(bbox.zmin() + bbox.zmax()) / 2)
    // };

    // with random point generated inside input mesh -> not reliable method still
    std::vector<Point_3> seeds = {
        random_point_inside_mesh(mesh_)
    };

  Mesh wrap;
  CGAL::alpha_wrap_3(mesh_, alpha, offset, wrap, CGAL::parameters::seed_points(std::ref(seeds)));

  t.stop();
  std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces, it took " << t.time() << " s." << std::endl;
  // std::cout << "Took " << t.time() << " s." << std::endl;

    // ------------------------------ write out output file -------------------------------
    std::string out_path = output_;
    // Remove last 4 characters (".off")
    if (out_path.size() > 4 && out_path.substr(out_path.size() - 4) == ".off") {
        out_path = out_path.substr(0, out_path.size() - 4);
    }

    // Add cuurent alpha and offset to the output name
    std::string out_name = out_path + "_inside_wrap.off";

    std::cout << "Writing to " << out_name << std::endl;
    CGAL::IO::write_polygon_mesh(out_name, wrap, CGAL::parameters::stream_precision(17));

    // std::string offset_output = "../data/Test_3D_Alphawrap/Output/3DBAG_Buildings/offset_testtest.off";
    // Mesh offset_mesh = offset_mesh_by_5cm(wrap, offset);

    // std::cout << "Writing to " << offset_output << std::endl;
    // CGAL::IO::write_polygon_mesh(offset_output, offset_mesh, CGAL::parameters::stream_precision(17));

    return wrap;
}

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

// --------------------------------OCTREE REFINEMENT-----------------------------------

struct OctreeCell {
    CGAL::Bbox_3 bbox;
    int depth = 0;
    std::vector<Mesh::Face_index> faces;
    std::array<std::unique_ptr<OctreeCell>, 8> children;

    OctreeCell() {
        for (auto& c : children)
            c = nullptr;
    }
};

void collect_faces( OctreeCell& cell, const Tree& tree) {
    cell.faces.clear();
    std::vector<Primitive::Id> hits;
    tree.all_intersected_primitives(cell.bbox, std::back_inserter(hits));

    for (auto id : hits)
        cell.faces.push_back(id);
}

double normal_variance( const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals) {
    K::Vector_3 mean(0,0,0);

    for (auto f : faces) {
        auto n = normals.at(f);
        n = n / std::sqrt(n.squared_length());
        mean += n;
    }

    mean /= faces.size();

    double variance = 0.0;
    for (auto f : faces) {
        auto n = normals.at(f);
        n = n / std::sqrt(n.squared_length());
        variance += (n - mean).squared_length();
    }

    return variance / faces.size();
}

bool has_sharp_feature( const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals, double angle_threshold_rad ) {
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& n1 = normals.at(faces[i]);

        for (size_t j = i + 1; j < faces.size(); ++j) {
            const auto& n2 = normals.at(faces[j]);

            double cos_angle =
                (n1 * n2) / std::sqrt(n1.squared_length() * n2.squared_length());

            cos_angle = std::clamp(cos_angle, -1.0, 1.0);

            double angle = std::acos(cos_angle);

            if (angle > angle_threshold_rad)
                return true;
        }
    }
    return false;
}

double distance_to_surface( const OctreeCell& cell, const Tree& tree ) {
    const auto& b = cell.bbox;

    K::Point_3 center(
        0.5 * (b.xmin() + b.xmax()),
        0.5 * (b.ymin() + b.ymax()),
        0.5 * (b.zmin() + b.zmax())
    );

    return std::sqrt(tree.squared_distance(center));
}

bool surface_is_simple( const std::vector<Mesh::Face_index>& faces, const std::map<Mesh::Face_index, K::Vector_3>& normals, double max_angle_rad ) {
    if (faces.size() < 2)
        return true;

    const auto& n0 = normals.at(faces[0]);

    for (auto f : faces) {
        const auto& n = normals.at(f);
        double cosang = (n * n0) /
            std::sqrt(n.squared_length() * n0.squared_length());

        cosang = std::clamp(cosang, -1.0, 1.0);

        if (std::acos(cosang) > max_angle_rad)
            return false;
    }
    return true;
}

bool needs_refinement( const OctreeCell& cell, const std::map<Mesh::Face_index, K::Vector_3>& face_normals) {
    constexpr int max_depth = 9;
    constexpr double cone_angle =
        10.0 * M_PI / 180.0; // strict!

    if (cell.depth >= max_depth)
        return false;

    // 1. Empty cell → stop
    if (cell.faces.empty())
        return false;

    // 2. Surface is locally simple → stop
    if (surface_is_simple(cell.faces, face_normals, cone_angle))
        return false;

    // 3. Otherwise → refine
    return true;
}

void subdivide(OctreeCell& cell)
{
    // ------------------------------------------------------------------
    // Preconditions (important for sanity)
    // ------------------------------------------------------------------
    const CGAL::Bbox_3& b = cell.bbox;

    const double dx = b.xmax() - b.xmin();
    const double dy = b.ymax() - b.ymin();
    const double dz = b.zmax() - b.zmin();

    // This MUST hold for a true octree
    // (you can comment this out later)
    if (std::abs(dx - dy) > 1e-9 || std::abs(dy - dz) > 1e-9) {
        std::cerr << "[subdivide] Warning: non-cubic cell at depth "
                  << cell.depth << " ("
                  << dx << ", " << dy << ", " << dz << ")\n";
    }

    // ------------------------------------------------------------------
    // Compute midpoints
    // ------------------------------------------------------------------
    const double mx = 0.5 * (b.xmin() + b.xmax());
    const double my = 0.5 * (b.ymin() + b.ymax());
    const double mz = 0.5 * (b.zmin() + b.zmax());

    // ------------------------------------------------------------------
    // Helper to create a child cell
    // ------------------------------------------------------------------
    auto make_child =
        [&](double xmin, double ymin, double zmin,
            double xmax, double ymax, double zmax)
    {
        auto child = std::make_unique<OctreeCell>();
        child->bbox  = CGAL::Bbox_3(xmin, ymin, zmin,
                                    xmax, ymax, zmax);
        child->depth = cell.depth + 1;
        return child;
    };

    // ------------------------------------------------------------------
    // Create the 8 octants (standard octree layout)
    // ------------------------------------------------------------------
    cell.children[0] = make_child(b.xmin(), b.ymin(), b.zmin(), mx, my, mz);
    cell.children[1] = make_child(mx,       b.ymin(), b.zmin(), b.xmax(), my, mz);
    cell.children[2] = make_child(mx,       my,       b.zmin(), b.xmax(), b.ymax(), mz);
    cell.children[3] = make_child(b.xmin(), my,       b.zmin(), mx, b.ymax(), mz);

    cell.children[4] = make_child(b.xmin(), b.ymin(), mz, mx, my, b.zmax());
    cell.children[5] = make_child(mx,       b.ymin(), mz, b.xmax(), my, b.zmax());
    cell.children[6] = make_child(mx,       my,       mz, b.xmax(), b.ymax(), b.zmax());
    cell.children[7] = make_child(b.xmin(), my,       mz, mx, b.ymax(), b.zmax());
}

void refine_cell(OctreeCell& cell, const Tree& tree, const std::map<Mesh::Face_index, K::Vector_3>& face_normals) {
    collect_faces(cell, tree);

    if (!needs_refinement(cell, face_normals))
        return;

    subdivide(cell);

    for (auto& child : cell.children)
        refine_cell(*child, tree, face_normals);
}

struct LeafCell {
    CGAL::Bbox_3 bbox;
};

void collect_leaves( const OctreeCell& cell, std::vector<LeafCell>& leaves ) {
    bool is_leaf = true;
    for (const auto& c : cell.children)
        if (c) is_leaf = false;

    if (is_leaf) {
        leaves.push_back({cell.bbox});
        return;
    }

    for (const auto& c : cell.children)
        if (c) collect_leaves(*c, leaves);
}

void write_off_wireframe( const std::vector<LeafCell>& cells, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: cannot open output file " << filename << "\n";
        return;
    }

    constexpr int VERTS_PER_CELL = 8;
    constexpr int FACES_PER_CELL = 6; // 6 quad faces per box

    const std::size_t numVertices = cells.size() * VERTS_PER_CELL;
    const std::size_t numFaces    = cells.size() * FACES_PER_CELL;

    // ------------------------------------------------------------------
    // OFF header
    // ------------------------------------------------------------------
    out << "OFF\n";
    out << numVertices << " " << numFaces << " 0\n";

    // ------------------------------------------------------------------
    // Write vertices
    // ------------------------------------------------------------------
    for (const auto& c : cells) {
        const auto& b = c.bbox;

        const double x0 = b.xmin(), x1 = b.xmax();
        const double y0 = b.ymin(), y1 = b.ymax();
        const double z0 = b.zmin(), z1 = b.zmax();

        out << x0 << " " << y0 << " " << z0 << "\n"; // 0
        out << x1 << " " << y0 << " " << z0 << "\n"; // 1
        out << x1 << " " << y1 << " " << z0 << "\n"; // 2
        out << x0 << " " << y1 << " " << z0 << "\n"; // 3
        out << x0 << " " << y0 << " " << z1 << "\n"; // 4
        out << x1 << " " << y0 << " " << z1 << "\n"; // 5
        out << x1 << " " << y1 << " " << z1 << "\n"; // 6
        out << x0 << " " << y1 << " " << z1 << "\n"; // 7
    }

    // ------------------------------------------------------------------
    // Write quad faces (boxes)
    // ------------------------------------------------------------------
    std::size_t v = 0;
    for (std::size_t i = 0; i < cells.size(); ++i) {

        // bottom
        out << "4 " << v+0 << " " << v+1 << " " << v+2 << " " << v+3 << "\n";
        // top
        out << "4 " << v+4 << " " << v+5 << " " << v+6 << " " << v+7 << "\n";

        // sides
        out << "4 " << v+0 << " " << v+1 << " " << v+5 << " " << v+4 << "\n";
        out << "4 " << v+1 << " " << v+2 << " " << v+6 << " " << v+5 << "\n";
        out << "4 " << v+2 << " " << v+3 << " " << v+7 << " " << v+6 << "\n";
        out << "4 " << v+3 << " " << v+0 << " " << v+4 << " " << v+7 << "\n";

        v += VERTS_PER_CELL;
    }

    out.close();
}

bool is_leaf(const OctreeCell& cell) {
    for (const auto& c : cell.children)
        if (c)
            return false;
    return true;
}

void collect_intersecting_leaves( const OctreeCell& cell, std::vector<LeafCell>& out ) {
    // Case 1: leaf cell
    if (is_leaf(cell)) {
        // Only keep it if it intersects the mesh
        if (!cell.faces.empty()) {
            out.push_back({cell.bbox});
        }
        return;
    }

    // Case 2: internal node → recurse
    for (const auto& c : cell.children) {
        if (c)
            collect_intersecting_leaves(*c, out);
    }
}

bool intersects_mesh_bbox(const CGAL::Bbox_3& bbox, const Tree& tree) {
    std::vector<Primitive::Id> hits;
    tree.all_intersected_primitives(bbox, std::back_inserter(hits));
    return !hits.empty();
}

void collect_cubical_leaves_at_depth( const OctreeCell& cell, int target_depth, std::vector<LeafCell>& out ) {
    // Leaf at the right depth
    if (cell.depth == target_depth) {
        bool is_leaf = true;
        for (const auto& c : cell.children)
            if (c) { is_leaf = false; break; }

        if (is_leaf && !cell.faces.empty()) {
            out.push_back({cell.bbox});
        }
        return;
    }

    // Recurse
    for (const auto& c : cell.children)
        if (c)
            collect_cubical_leaves_at_depth(*c, target_depth, out);
}

void find_max_depth(const OctreeCell& cell, int& max_depth) {
    max_depth = std::max(max_depth, cell.depth);
    for (const auto& c : cell.children)
        if (c)
            find_max_depth(*c, max_depth);
}

void find_max_intersecting_leaf_depth( const OctreeCell& cell, int& max_depth ) {
    bool is_leaf = true;
    for (const auto& c : cell.children)
        if (c) { is_leaf = false; break; }

    if (is_leaf) {
        if (!cell.faces.empty()) {
            max_depth = std::max(max_depth, cell.depth);
        }
        return;
    }

    for (const auto& c : cell.children)
        if (c)
            find_max_intersecting_leaf_depth(*c, max_depth);
}

void collect_finest_intersecting_leaves( const OctreeCell& cell, int target_depth, std::vector<LeafCell>& out ) {
    bool is_leaf = true;
    for (const auto& c : cell.children)
        if (c) { is_leaf = false; break; }

    if (is_leaf) {
        if (cell.depth == target_depth && !cell.faces.empty()) {
            out.push_back({cell.bbox});
        }
        return;
    }

    for (const auto& c : cell.children)
        if (c)
            collect_finest_intersecting_leaves(*c, target_depth, out);
}

void collect_intersecting_leaves_verified( const OctreeCell& cell, const Tree& tree, std::vector<LeafCell>& out ){
    bool is_leaf = true;
    for (const auto& c : cell.children)
        if (c) { is_leaf = false; break; }

    if (is_leaf) {
        if (intersects_mesh_bbox(cell.bbox, tree)) {
            out.push_back({cell.bbox});
        }
        return;
    }

    for (const auto& c : cell.children)
        if (c)
            collect_intersecting_leaves_verified(*c, tree, out);
}

// ----------------------------------------------------------------------------------
// -------------------------------------MAIN-----------------------------------------
// ----------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Read the input
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("../data/Test_3D_Alphawrap/Input/demo/joep_huis.off");
    std::cout << "------------------------------------------------------------" << std::endl;
  std::cout << "Reading input: " << filename << std::endl;
    // demo 1
  const double relative_alpha = 70.; //2000. //20. //1000.
  const double relative_offset = 1500.; // 7000. //600. //12000.


    // create output name (+ folder if folder does not exist)
    std::string output_name = generate_output_name(filename, relative_alpha, relative_offset);

    std::filesystem::path p(output_name);
    std::filesystem::create_directory(p.parent_path());

    // mesh the input
    Mesh mesh;
    if(!PMP::IO::read_polygon_mesh(filename, mesh) || is_empty(mesh) || !is_triangle_mesh(mesh))
    {
        //
        return EXIT_FAILURE;
    }

    std::map<Mesh::Face_index, K::Vector_3> face_normals;

    CGAL::Polygon_mesh_processing::compute_face_normals(
        mesh,
        boost::make_assoc_property_map(face_normals)
    );


    Tree tree(faces(mesh).begin(), faces(mesh).end(), mesh);
    tree.accelerate_distance_queries();

    // ----------------------------IS MESH VALID?------------------------------------
    // std::cout << "Testing if input mesh is valid:" << std::endl;
    valid_mesh_boolean(mesh);
    std::cout << "------------------------------------------------------------" << std::endl;

    // demo 2
    // ------------------ ALPHA WRAP FROM THE OUTSIDE-----------------------------
    std::cout << "-----------------3D ALPHA WRAPPING THE INPUT:-------------" << std::endl;
    Mesh alpha_wrap = _3D_alpha_wrap( output_name,relative_alpha,relative_offset, mesh);

    // is mesh valid?
    std::cout << "Testing if 3D alpha wrapped mesh is valid:" << std::endl;
    valid_mesh_boolean(alpha_wrap);
    std::cout << "------------------------------------------------------------" << std::endl;

    // demo 3
        // ----------------- SHARPENING EDGES --------------------
    std::string weight_output_wrap =
        "../data/Test_3D_Alphawrap/Output/demo/refined.ply";
    // std::string weight_output_inner_wrap =
    //     "../data/Test_3D_Alphawrap/Output/3DBAG_Buildings/d_innerwrap_input_bk.ply";

    // // offset wrapped from inside mesh
    double _2_offset = rel_offset_to_offset(mesh, relative_offset);
    double offset = (rel_offset_to_offset(mesh, relative_offset))/2;
    // Mesh offset_inner_wrap = offset_mesh(alpha_inside_wrap, _2_offset);
    // std::cout << "Succesfully offsetted wrapped from inside mesh" << std::endl;

    Mesh offset_input = offset_mesh(mesh, offset);

    // tag + colorise outer wrap
    CGAL::Real_timer t;
    t.start();
    Mesh distance_alpha = add_midpoint_distance_tag(alpha_wrap, offset_input, 220);
    // std::cout << "Succesfully adding midpoint distance to alpha wrap mesh" << std::endl;

    // std::cout << "Writing to " << weight_output_wrap << std::endl;
    {
      std::ofstream out(weight_output_wrap, std::ios::binary);
      if (!out) {
          std::cerr << "Cannot open " << weight_output_wrap << " for writing\n";
      } else {
          CGAL::IO::write_PLY(
              out,
              distance_alpha,
              CGAL::parameters::stream_precision(17)
          );
      }
    }
    std::string refining =
    "../data/Test_3D_Alphawrap/Output/demo/sharpe_concave_corners.ply";
    Mesh refined_edges = refine_round_edges(distance_alpha, offset_input);
    t.stop();
    std::cout << "Took " << t.time() << " s." << std::endl;
    std::cout << "Succesfully refined edges" << std::endl;

    std::cout << "Writing to " << refining << std::endl;
    {
      std::ofstream out(refining, std::ios::binary);
      if (!out) {
          std::cerr << "Cannot open " << refining << " for writing\n";
      } else {
          CGAL::IO::write_PLY(
              out,
              refined_edges,
              CGAL::parameters::stream_precision(17)
          );
      }
    }
    // ----------------------------IS MESH VALID?------------------------------------
    // std::cout << "Testing if input mesh is valid:" << std::endl;
    valid_mesh_boolean(refined_edges);
    std::cout << "------------------------------------------------------------" << std::endl;

    // // -------------------------------OCTREE REFINEMENT---------------------------
    // // Create root octree cell
    // OctreeCell root;
    // // compute bbox
    // CGAL::Bbox_3 tight = CGAL::Polygon_mesh_processing::bbox(mesh);
    //
    // // Center
    // double cx = 0.5 * (tight.xmin() + tight.xmax());
    // double cy = 0.5 * (tight.ymin() + tight.ymax());
    // double cz = 0.5 * (tight.zmin() + tight.zmax());
    //
    // // Largest half-extent
    // double half = std::max({
    //     0.5 * (tight.xmax() - tight.xmin()),
    //     0.5 * (tight.ymax() - tight.ymin()),
    //     0.5 * (tight.zmax() - tight.zmin())
    // });
    //
    // // Optional looseness
    // half *= 1.1;
    //
    // // TRUE cube
    // root.bbox = CGAL::Bbox_3(
    //     cx - half, cy - half, cz - half,
    //     cx + half, cy + half, cz + half
    // );
    // root.depth = 0;
    //
    // // Run Octree refinement
    // refine_cell(root, tree, face_normals);
    //
    // // Collect all leaf cells
    // std::vector<LeafCell> leaves_all; // use these if you want all cubes
    // collect_leaves(root, leaves_all);
    //
    // // Collect leaf cells that intersect the input
    // std::vector<LeafCell> leaves;
    // collect_intersecting_leaves(root, leaves);
    //
    // // Collect finest level of detail leaves
    // int finest_depth = 0;
    // find_max_intersecting_leaf_depth(root, finest_depth);
    //
    // std::vector<LeafCell> cubes;
    // collect_finest_intersecting_leaves(root, finest_depth, cubes);
    //
    // // collect intersecting cubes 2.0
    // std::vector<LeafCell> cubes2;
    // collect_intersecting_leaves_verified(root, tree, cubes2);
    //
    // std::cout << "Leaf cells: " << cubes2.size() << std::endl; // can also use 'leaves_all' or 'leaves' or 'cubes'
    //
    // // Write .off wireframe file
    // std::string octree_refinement = "../data/Test_3D_Alphawrap/Output/3DBAG_Buildings/Octree_refinement.off";
    //
    // write_off_wireframe(cubes2, octree_refinement); // can also use 'leaves_all' or 'leaves' or 'cubes'
    //
    // std::cout << "Octree wireframe written to:\n"
    //           << octree_refinement << std::endl;



  // // ------------------ ALPHA WRAP FROM THE INSIDE -----------------------
  //   std::cout << "--------3D ALPHA WRAPPING THE INPUT FROM INSIDE:----------" << std::endl;
  // Mesh alpha_inside_wrap = _3D_alpha_inside_wrap( output_name,relative_alpha,relative_offset, mesh);

    // // is mesh valid?
    // std::cout << "Testing if from inside 3D alpha wrapped mesh is valid:" << std::endl;
    // valid_mesh_boolean(alpha_inside_wrap);
  // std::cout << "------------------------------------------------------------" << std::endl;


    //
    // // tag + colorise inner wrap
    // Mesh distance_alpha_inner = add_midpoint_distance_tag(offset_inner_wrap, offset_input, 220);
    // std::cout << "Succesfully adding midpoint distance to alpha wrap from inside mesh" << std::endl;
    //
    // std::cout << "Writing to " << weight_output_inner_wrap << std::endl;
    // {
    //   std::ofstream out(weight_output_inner_wrap, std::ios::binary);
    //   if (!out) {
    //       std::cerr << "Cannot open " << weight_output_inner_wrap << " for writing\n";
    //   } else {
    //       CGAL::IO::write_PLY(
    //           out,
    //           distance_alpha_inner,
    //           CGAL::parameters::stream_precision(17)
    //       );
    //   }
    // }

  return EXIT_SUCCESS;
}