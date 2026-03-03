//
// Created by Michel Beeren on 20/02/2026.
//

#include "hausdorff.h"
#include <CGAL/IO/Color.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Real_timer.h>
#include <CGAL/alpha_wrap_3.h>

// using namespace CGAL;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
using Vector_3 = K::Vector_3;

void hausdorff_distance(Mesh& wrapped_mesh, const Mesh& original_mesh, const double filter)
{
    using fd = Mesh::Face_index;
    using vd = Mesh::Vertex_index;
    using Color = CGAL::IO::Color;
    using hd = Mesh::Halfedge_index;

    Mesh mesh_to_tag = wrapped_mesh;

    // AABB tree on original mesh
    typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
    typedef CGAL::AABB_traits<K, Primitive>                AABB_traits;
    typedef CGAL::AABB_tree<AABB_traits>                   Tree;

    Tree tree(faces(original_mesh).begin(),
              faces(original_mesh).end(),
              original_mesh);
    tree.accelerate_distance_queries();

    // Distance property
    typedef Mesh::Property_map<fd, double> FaceWeightMap;
    FaceWeightMap fweight;
    bool created;
    boost::tie(fweight, created) =
        mesh_to_tag.add_property_map<fd, double>("f:dist_to_original", 0.0);

    double dmin = std::numeric_limits<double>::max();
    double dmax = 0.0;

    // 1) Compute distances and min/max
    for (fd f : mesh_to_tag.faces())
    {
        auto h  = mesh_to_tag.halfedge(f);
        vd v0   = mesh_to_tag.source(h);
        vd v1   = mesh_to_tag.target(h);
        vd v2   = mesh_to_tag.target(mesh_to_tag.next(h));

        const Point_3& p0 = mesh_to_tag.point(v0);
        const Point_3& p1 = mesh_to_tag.point(v1);
        const Point_3& p2 = mesh_to_tag.point(v2);

        Point_3 c( (p0.x() + p1.x() + p2.x()) / 3.0,
                   (p0.y() + p1.y() + p2.y()) / 3.0,
                   (p0.z() + p1.z() + p2.z()) / 3.0 );

        double dist = std::sqrt(tree.squared_distance(c));
        fweight[f] = dist;

        dmin = std::min(dmin, dist);
        dmax = std::max(dmax, dist);
    }

    // **Ensure the color property map exists on faces**
    Mesh::Property_map<fd, Color> fcolor;
    bool created_color;
    boost::tie(fcolor, created_color) =
        wrapped_mesh.add_property_map<fd, Color>("f:color", Color(255, 255, 255));  // Add color property map to faces

    // If the property map creation fails
    if (!created_color) {
        std::cerr << "Error: Failed to add color property map to faces.\n";
        return;
    }

    double range = (dmax > dmin) ? (dmax - dmin) : 1.0;

    // Map distance -> blue to red (0 -> blue, 0.05 -> red)
    for (fd f : mesh_to_tag.faces())
    {
        double d = fweight[f];
        double t = (d - dmin) / range;  // Normalize distance between 0 and 1

        // Apply gradient from blue to red using RGB
        unsigned char r, g, b;
        double min = 0.03;
        double max = 0.05;

        // Blue for distances smaller than 0.03
        if (d < min) {
            r = 0;
            g = 0;
            b = 255;
        }
        // For distances between 0.03 and 0.05, create a smooth gradient
        else if (d >= min && d <= max) {
            t = (d - min) / (max - min);  // Interpolate in this range
            r = static_cast<unsigned char>(255 * t);  // Red increases as distance increases
            g = 0;                                      // Green stays at 0
            b = static_cast<unsigned char>(255 * (1.0 - t)); // Blue decreases as distance increases
        }
        // Red for distances greater than 0.05
        else {
            r = 255;
            g = 0;
            b = 0;
        }

        fcolor[f] = Color(r, g, b);
    }
}

// Function to compute Hausdorff distance and apply color tagging
void hausdorff_distance2(Mesh& wrapped_mesh, const Mesh& original_mesh, const double filter)
{
    using fd = Mesh::Face_index;
    using vd = Mesh::Vertex_index;
    using Color = CGAL::IO::Color;
    using hd = Mesh::Halfedge_index;

    // Create the f:color property map on the mesh if it doesn't exist
    Mesh::Property_map<fd, Color> fcolor;
    bool created_color;
    boost::tie(fcolor, created_color) =
        wrapped_mesh.add_property_map<fd, Color>("f:color", Color(255, 255, 255));  // Add color property map to faces

    // If the property map creation fails
    if (!created_color) {
        std::cerr << "Error: Failed to add color property map to faces.\n";
        return;
    }

    Mesh mesh_to_tag = wrapped_mesh;

    // AABB tree on original mesh
    typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
    typedef CGAL::AABB_traits<K, Primitive>                AABB_traits;
    typedef CGAL::AABB_tree<AABB_traits>                   Tree;

    Tree tree(faces(original_mesh).begin(),
              faces(original_mesh).end(),
              original_mesh);
    tree.accelerate_distance_queries();

    // Distance property
    typedef Mesh::Property_map<fd, double> FaceWeightMap;
    FaceWeightMap fweight;
    bool created;
    boost::tie(fweight, created) =
        mesh_to_tag.add_property_map<fd, double>("f:dist_to_original", 0.0);

    double dmin = std::numeric_limits<double>::max();
    double dmax = 0.0;

    // Compute distances and min/max
    for (fd f : mesh_to_tag.faces())
    {
        auto h  = mesh_to_tag.halfedge(f);
        vd v0   = mesh_to_tag.source(h);
        vd v1   = mesh_to_tag.target(h);
        vd v2   = mesh_to_tag.target(mesh_to_tag.next(h));

        const Point_3& p0 = mesh_to_tag.point(v0);
        const Point_3& p1 = mesh_to_tag.point(v1);
        const Point_3& p2 = mesh_to_tag.point(v2);

        Point_3 c( (p0.x() + p1.x() + p2.x()) / 3.0,
                   (p0.y() + p1.y() + p2.y()) / 3.0,
                   (p0.z() + p1.z() + p2.z()) / 3.0 );

        double dist = std::sqrt(tree.squared_distance(c));
        fweight[f] = dist;

        dmin = std::min(dmin, dist);
        dmax = std::max(dmax, dist);
    }

    double range = (dmax > dmin) ? (dmax - dmin) : 1.0;

    // Map distance -> dark blue to red (0 -> dark blue, 1 -> red)
    for (fd f : mesh_to_tag.faces())
    {
        double d = fweight[f];
        double t = (d - dmin) / range;  // 0 = closest, 1 = farthest

        // Dark Blue (0,0,255) -> Red (255,0,0)
        unsigned char r = static_cast<unsigned char>(255 * t);  // Red increases with distance
        unsigned char g = 0;                                      // Green stays at 0
        unsigned char b = static_cast<unsigned char>(255 * (1.0 - t)); // Blue decreases with distance
        fcolor[f] = Color(r, g, b);
    }
}

void write_hausdorff_distance(const Mesh& wrapped_mesh, const std::string& filename)
{
    // Create output name
    std::string out_path = filename;
    // Remove last 4 characters (".off") if the file name ends with .off
    if (out_path.size() > 4 && out_path.substr(out_path.size() - 4) == ".off") {
        out_path = out_path.substr(0, out_path.size() - 4);
    }

    // Add current alpha and offset to the output name
    std::string out_name = out_path + "_HausDorff.ply";

    // Create a property map to store color values for faces
    typedef Mesh::Property_map<Mesh::Face_index, CGAL::IO::Color> FaceColorMap;

    // Use std::optional to capture the property map
    std::optional<FaceColorMap> face_colors_opt = wrapped_mesh.property_map<Mesh::Face_index, CGAL::IO::Color>("f:color");

    // Check if the property map exists
    if (!face_colors_opt)
    {
        std::cerr << "Error: No color property on faces.\n";
        return;
    }

    // If the map exists, use it
    FaceColorMap face_colors = *face_colors_opt;

    // Open the file for writing (in ASCII mode)
    std::ofstream out(out_name);
    if (!out) {
        std::cerr << "Cannot open " << out_name << " for writing\n";
        return;
    }

    // Write the mesh to an ASCII PLY file
    CGAL::IO::write_PLY(out, wrapped_mesh, CGAL::parameters::stream_precision(17));

    std::cout << "PLY file saved: " << out_name << std::endl;
}