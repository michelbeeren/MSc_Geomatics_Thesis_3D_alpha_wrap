//
// Created by Michel Beeren on 04/05/2026.
//

#ifndef THESIS_ALPHA_WRAPPING_H
#define THESIS_ALPHA_WRAPPING_H

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <limits>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <vector>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Real_timer.h>
#include <CGAL/alpha_wrap_3.h>
#include <map>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/helpers.h>
#include <filesystem>
#include <algorithm>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <stdexcept>

namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;

using Point_container = std::vector<Point_3>;
using Vector_3 = K::Vector_3;
using Segment_3 = K::Segment_3;
using Ray_3 = K::Ray_3;

using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using AABB_traits = CGAL::AABB_traits<K, Primitive>;
using Tree = CGAL::AABB_tree<AABB_traits>;

using face_descriptor = Mesh::Face_index;
struct MeshData {
    Mesh mesh;  // CGAL Surface_mesh type
    std::map<Mesh::Face_index, K::Vector_3> face_normals;  // Normals of the faces indexed by Face_index
    std::unique_ptr<Tree> tree;  // Optional octree structure
    std::map<Mesh::Face_index, std::set<Mesh::Face_index>> adjacency_map;  // Adjacency map of faces indexed by Face_index
};
enum class Input_kind
{
    Triangle_mesh,
    Point_cloud,
    Unknown
};

std::map<Mesh::Face_index, std::set<Mesh::Face_index>> create_adjacency_map(const Mesh& mesh) {
    std::map<Mesh::Face_index, std::set<Mesh::Face_index>> adjacency_map;

    for (auto f : mesh.faces()) {
        for (auto h : mesh.halfedges_around_face(mesh.halfedge(f))) {
            auto opp_h = mesh.opposite(h);
            if (!mesh.is_border(opp_h)) {
                Mesh::Face_index neighbor = mesh.face(opp_h);
                adjacency_map[f].insert(neighbor);
            }
        }
    }
    return adjacency_map;
}

// --------------------------------------MESH INPUT-------------------------------------
MeshData mesh_input(const std::string& filename, bool compute_normals, bool build_tree)
{
    // meshes input file, optionally also computes normals and builds a tree
    MeshData out;

    if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, out.mesh) ||
        CGAL::is_empty(out.mesh) ||
        !CGAL::is_triangle_mesh(out.mesh))
    {
        throw std::runtime_error("Failed to read mesh: " + filename);
    }

    if (compute_normals) {
        CGAL::Polygon_mesh_processing::compute_face_normals(
            out.mesh,
            boost::make_assoc_property_map(out.face_normals)
        );
    }

    if (build_tree) {
        out.tree = std::make_unique<Tree>(
            faces(out.mesh).begin(),
            faces(out.mesh).end(),
            out.mesh
        );
        out.tree->accelerate_distance_queries();

        // defensive check (should never trigger, but matches your wish)
        if (!out.tree) {
            throw std::runtime_error("Tree not built for: " + filename);
        }
    }

    // Create adjacency map
    out.adjacency_map = create_adjacency_map(out.mesh);

    std::cout << "Input file successfully meshed" << std::endl;
    return out;
}

Point_container point_input(const std::string& filename)
{
    Point_container points;

    if(!CGAL::IO::read_points(filename, std::back_inserter(points)) || points.empty())
    {
        throw std::runtime_error("Failed to read point cloud: " + filename);
    }

    std::cout << "Input: " << points.size() << " points" << std::endl;
    return points;
}

Input_kind detect_input_kind(const std::string& filename)
{
    {
        Mesh mesh;
        if(PMP::IO::read_polygon_mesh(filename, mesh) &&
           !CGAL::is_empty(mesh) &&
           CGAL::is_triangle_mesh(mesh))
        {
            return Input_kind::Triangle_mesh;
        }
    }

    {
        Point_container points;
        if(CGAL::IO::read_points(filename, std::back_inserter(points)) &&
           !points.empty())
        {
            return Input_kind::Point_cloud;
        }
    }

    return Input_kind::Unknown;
}

double smallest_point_to_point_distance_from_off(const std::string& off_filename)
{
    auto is_comment_or_empty = [](const std::string& s) -> bool
    {
        const std::size_t first = s.find_first_not_of(" \t\r\n");
        return (first == std::string::npos || s[first] == '#');
    };

    std::ifstream in(off_filename);
    if(!in)
    {
        throw std::runtime_error("Failed to open OFF file: " + off_filename);
    }

    std::string line;
    std::string header;

    std::size_t num_vertices = 0;
    std::size_t num_faces = 0;
    std::size_t num_edges = 0;
    bool got_counts = false;

    while(std::getline(in, line))
    {
        if(is_comment_or_empty(line))
            continue;

        std::istringstream ls(line);
        ls >> header;
        if(!header.empty()) {
            if(header == "OFF" && (ls >> num_vertices >> num_faces >> num_edges))
                got_counts = true;
            break;
        }
    }

    if(header != "OFF")
    {
        throw std::runtime_error("Invalid OFF header in file: " + off_filename);
    }

    while(!got_counts && std::getline(in, line))
    {
        if(is_comment_or_empty(line))
            continue;

        std::istringstream ls(line);
        if(ls >> num_vertices >> num_faces >> num_edges)
        {
            got_counts = true;
            break;
        }
    }

    if(!got_counts)
    {
        throw std::runtime_error("Missing OFF counts line in file: " + off_filename);
    }
    (void)num_faces;
    (void)num_edges;

    std::vector<Point_3> points;
    points.reserve(num_vertices);

    while(points.size() < num_vertices && std::getline(in, line))
    {
        if(is_comment_or_empty(line))
            continue;

        std::istringstream ls(line);
        double x = 0.0, y = 0.0, z = 0.0;
        if(!(ls >> x >> y >> z))
            continue;

        points.emplace_back(x, y, z);
    }

    if(points.size() != num_vertices)
    {
        throw std::runtime_error("OFF vertex section is incomplete in file: " + off_filename);
    }

    if(points.size() < 2)
    {
        throw std::runtime_error("Need at least 2 points to compute a distance: " + off_filename);
    }

    double min_sq_dist = std::numeric_limits<double>::infinity();
    for(std::size_t i = 0; i + 1 < points.size(); ++i)
    {
        for(std::size_t j = i + 1; j < points.size(); ++j)
        {
            const double d2 = CGAL::to_double(CGAL::squared_distance(points[i], points[j]));
            if(d2 < min_sq_dist)
                min_sq_dist = d2;
        }
    }

    return std::sqrt(min_sq_dist);
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

    // Remove last 4 characters (".obj")
    if (out_path.size() > 4 && out_path.substr(out_path.size() - 4) == ".obj") {
        out_path = out_path.substr(0, out_path.size() - 4);
    }

    // Add cuurent alpha and offset to the output name
    std::string out_name = out_path + "_a=" + rel_alpha_txt + "_offset=" + rel_offset_txt + ".off";
    return out_name;
}

double upper_bound_max_d_to_input(const double alpha, const double offset, Mesh& input, Mesh& output, double tau) {
    double upper_bound_max_d_to_input = 0.0;
    double t = (tau-1)*offset;
    if (t >= 0. && t <= (4./15.)*alpha) {
        upper_bound_max_d_to_input = offset + (2./3.)*alpha + t/2.;
    }
    else if (t > (4./15.)*alpha && t < (2./3.)*alpha) {
        upper_bound_max_d_to_input = offset + (3*alpha + 9*t + std::sqrt(10*alpha*alpha - (alpha + 3*t)*(alpha + 3*t)))/10.;
    }
    else if (t >= (2./3.)*alpha) {
        upper_bound_max_d_to_input = offset + alpha;
    }

    double d_output_to_input = PMP::approximate_Hausdorff_distance<CGAL::Sequential_tag>( output, input, CGAL::parameters::number_of_points_per_area_unit(100));
    std::cout << "Hausdorff distance = " << d_output_to_input << std::endl;
    return upper_bound_max_d_to_input;
}

Mesh mod_alpha_wrap_tr_mesh(const std::string filename, const double relative_alpha_, const double relative_offset_, MeshData& data_, const double max_d, bool write_out_, bool validate) {
    Mesh& mesh_ = data_.mesh;
    double max_d_to_input_in_offsets_ = max_d;

    // compute alpha and offset from a_rel and d_rel and bbox
    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh_);
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                         CGAL::square(bbox.ymax() - bbox.ymin()) +
                                         CGAL::square(bbox.zmax() - bbox.zmin()));
    const double alpha = diag_length / relative_alpha_;
    const double offset = diag_length / relative_offset_;
    std::cout << "--------------------3D ALPHA WRAPPING THE INPUT:----------------" << std::endl;
    std::cout << "alpha = " << alpha << " (a_rel = " << relative_alpha_ << ") and offset = " << offset << " (offset_rel = " << relative_offset_ << ")" << std::endl;

    // Construct the wrap
    CGAL::Real_timer t;
    t.start();

    Mesh wrap;

    std::cout << "..........Running alpha wrap algorithm with other refinement rule and steiner point placement.........." << std::endl;
    CGAL::alpha_wrap_3(mesh_, alpha, offset, wrap, CGAL::parameters::max_distance_to_input_in_offsets(max_d_to_input_in_offsets_));

    t.stop();
    std::cout << "🎁 Successfully alpha wrapped! 🎁: " << num_vertices(wrap) << " 🔘vertices🔘, " << num_faces(wrap) << " 📐faces📐, it took " << t.time() << " s.⏰" << std::endl;

    std::string output_ = generate_output_name(filename, relative_alpha_, relative_offset_);
    // Write the output mesh
    double upper_bound = alpha + offset;
    if (write_out_) {
        std::filesystem::path p(output_);
        std::filesystem::create_directories(p.parent_path());
        if (max_d_in_offets) {
            if (output_.size() > 4 && output_.substr(output_.size() - 4) == ".off") {
                output_ = output_.substr(0, output_.size() - 4);
            }
            std::ostringstream oss;
            if (max_d_to_input_in_offsets_ <= 1) {
                max_d_to_input_in_offsets_ = 1.5;
            }
            oss << std::fixed << std::setprecision(1) << max_d_to_input_in_offsets_;
            const std::string refined = oss.str();
            output_ += "_refined=" + refined + ".off";
            upper_bound = upper_bound_max_d_to_input(alpha, offset, mesh_, wrap, max_d_to_input_in_offsets_);
            std::cout << "upper_bound = " << upper_bound << std::endl;
            double g = max_d_to_input_in_offsets_*offset;
            std::cout << "g = " << g << std::endl;
        }
        if (!max_d_in_offets) {
            upper_bound_max_d_to_input(alpha, offset, mesh_, wrap, max_d_in_offets);
            std::cout << "upper_bound = " << upper_bound << std::endl;
        }
        std::cout << "📝 Writing 📝 to: " << output_ << std::endl;
        CGAL::IO::write_polygon_mesh(output_, wrap, CGAL::parameters::stream_precision(25));
        double dmin = smallest_point_to_point_distance_from_off(output_);
        double lower_bound_tr_sphere = (max_d_to_input_in_offsets_-1)*offset;
        // std::cout << "alpha = " << alpha << " , offset = " << offset << " , Beeren method lowest = " << lower_bound_tr_sphere << std::endl;
        double lower_bound = std::min({lower_bound_tr_sphere, offset, alpha});
        // std::cout << "prooven lower bound = " << lower_bound << std::endl;
        std::cout << "min distance between output points = " << dmin << " => " << lower_bound << " (proven lower bound)" << std::endl;
    }

    // ----------------------------- validate the output ------------------------------
    if (validate) {
        valid_mesh_boolean(wrap);
    }
    std::cout << "---------------------------------------------------------------" << std::endl;
    return wrap;
}

Mesh mod_alpha_wrap_pc(const std::string filename, const double relative_alpha_, const double relative_offset_, bool max_d_in_offets, const double max_d, bool write_out_, bool validate) {
    double max_d_to_input_in_offsets_ = max_d;

    Point_container points_;
    if(!CGAL::IO::read_points(filename, std::back_inserter(points_)) || points_.empty())
    {
        std::cerr << "Invalid input: " << filename << std::endl;
    }

    if (points_.empty())
    {
        throw std::runtime_error("Point cloud is empty.");
    }

    // compute alpha and offset from a_rel and d_rel and bbox
    CGAL::Bbox_3 bbox = CGAL::bbox_3(points_.begin(), points_.end());
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                         CGAL::square(bbox.ymax() - bbox.ymin()) +
                                         CGAL::square(bbox.zmax() - bbox.zmin()));
    const double alpha = diag_length / relative_alpha_;
    const double offset = diag_length / relative_offset_;
    std::cout << "--------------------3D ALPHA WRAPPING THE INPUT:----------------" << std::endl;
    std::cout << "alpha = " << alpha << " (a_rel = " << relative_alpha_ << ") and offset = " << offset << " (offset_rel = " << relative_offset_ << ")" << std::endl;

    // Construct the wrap
    CGAL::Real_timer t;
    t.start();

    Mesh wrap;
    if (max_d_in_offets) {
        std::cout << "..........Running alpha wrap algorithm with other refinement rule and steiner point placement.........." << std::endl;
        CGAL::alpha_wrap_3(points_, alpha, offset, wrap, CGAL::parameters::max_distance_to_input_in_offsets(max_d_to_input_in_offsets_));
    }
    if (!max_d_in_offets) {
        std::cout << "..........Running normal alpha wrap algorithm.........." << std::endl;
        CGAL::alpha_wrap_3(points_, alpha, offset, wrap);
    }
    t.stop();
    std::cout << "🎁 Successfully alpha wrapped! 🎁: " << num_vertices(wrap) << " 🔘vertices🔘, " << num_faces(wrap) << " 📐faces📐, it took " << t.time() << " s.⏰" << std::endl;

    std::string output_ = generate_output_name(filename, relative_alpha_, relative_offset_);
    // Write the output mesh
    if (write_out_) {
        std::filesystem::path p(output_);
        std::filesystem::create_directories(p.parent_path());
        if (max_d_in_offets) {
            if (output_.size() > 4 && output_.substr(output_.size() - 4) == ".off") {
                output_ = output_.substr(0, output_.size() - 4);
            }
            std::ostringstream oss;
            if (max_d_to_input_in_offsets_ <= 1) {
                max_d_to_input_in_offsets_ = 1.5;
            }
            oss << std::fixed << std::setprecision(1) << max_d_to_input_in_offsets_;
            const std::string refined = oss.str();
            output_ += "_refined=" + refined + ".off";
        }
        std::cout << "📝 Writing 📝 to: " << output_ << std::endl;
        CGAL::IO::write_polygon_mesh(output_, wrap, CGAL::parameters::stream_precision(25));
    }

    // ----------------------------- validate the output ------------------------------
    if (validate) {
        valid_mesh_boolean(wrap);
    }
    std::cout << "---------------------------------------------------------------" << std::endl;
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

#endif //THESIS_ALPHA_WRAPPING_H
