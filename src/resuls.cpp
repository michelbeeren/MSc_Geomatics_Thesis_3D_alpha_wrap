//
// Created by Michel Beeren on 13/05/2026.
//

#include "resuls.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include "val3dity.h"

namespace PMP = CGAL::Polygon_mesh_processing;

namespace
{
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using AABB_traits = CGAL::AABB_traits_3<K, Primitive>;
using Tree = CGAL::AABB_tree<AABB_traits>;

void ensure_triangle_mesh_for_distance(const Mesh& mesh, const std::string& mesh_name)
{
    if (CGAL::is_empty(mesh)) {
        throw std::invalid_argument(mesh_name + " is empty");
    }
    if (!CGAL::is_triangle_mesh(mesh)) {
        throw std::invalid_argument(mesh_name + " is not a triangle mesh");
    }
}

unsigned int to_unsigned_sample_count(std::size_t number_of_samples)
{
    if (number_of_samples == 0) {
        return 0;
    }

    return static_cast<unsigned int>(
        std::min<std::size_t>(number_of_samples, std::numeric_limits<unsigned int>::max()));
}

void ensure_positive_finite(const double value, const std::string& name)
{
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(name + " must be positive and finite");
    }
}

void ensure_triangle_mesh_in_place(Mesh& mesh, const std::string& mesh_name)
{
    if (CGAL::is_empty(mesh)) {
        throw std::invalid_argument(mesh_name + " is empty");
    }

    if (!CGAL::is_triangle_mesh(mesh)) {
        PMP::triangulate_faces(mesh);
    }

    if (!CGAL::is_triangle_mesh(mesh)) {
        throw std::invalid_argument(mesh_name + " could not be triangulated");
    }
}

double bbox_diagonal(const Mesh& mesh)
{
    const CGAL::Bbox_3 bbox = PMP::bbox(mesh);
    const double dx = bbox.xmax() - bbox.xmin();
    const double dy = bbox.ymax() - bbox.ymin();
    const double dz = bbox.zmax() - bbox.zmin();
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::string statistics_output_name(
    const double relative_alpha,
    const double relative_offset,
    const double tau,
    const bool use_beeren_method)
{
    const int rel_alpha_as_int = static_cast<int>(relative_alpha);
    const int rel_offset_as_int = static_cast<int>(relative_offset);

    std::string output = "../data/Output/statistics/wrap_a=" +
                         std::to_string(rel_alpha_as_int) +
                         "_offset=" + std::to_string(rel_offset_as_int);

    if (use_beeren_method) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << tau;
        output += "_refined=" + oss.str();
    }

    output += ".off";
    return output;
}

std::vector<double> distances_to_mesh(
    const std::vector<Point_3>& sample_points,
    const Mesh& target_mesh)
{
    ensure_triangle_mesh_for_distance(target_mesh, "target_mesh");

    Tree tree(faces(target_mesh).begin(), faces(target_mesh).end(), target_mesh);
    tree.accelerate_distance_queries();

    std::vector<double> distances;
    distances.reserve(sample_points.size());

    for (const Point_3& p : sample_points) {
        distances.push_back(std::sqrt(tree.squared_distance(p)));
    }

    return distances;
}

double mean_of_distances(const std::vector<double>& distances)
{
    if (distances.empty()) {
        return 0.0;
    }

    const double total = std::accumulate(distances.begin(), distances.end(), 0.0);
    return total / static_cast<double>(distances.size());
}

double max_of_distances(const std::vector<double>& distances)
{
    if (distances.empty()) {
        return 0.0;
    }

    return *std::max_element(distances.begin(), distances.end());
}

int validate_per_connected_component_binary(const Mesh& mesh)
{
    std::vector<Mesh> components;
    PMP::split_connected_components(mesh, components);

    if (components.empty()) {
        return 0;
    }

    bool all_valid = true;
    for (const Mesh& component : components) {
        all_valid = all_valid && valid_mesh_boolean(component);
    }

    return all_valid ? 1 : 0;
}

constexpr std::size_t k_default_statistics_sample_count = 10000;

std::string optional_int_to_csv_cell(const std::optional<int>& value)
{
    return value ? std::to_string(*value) : "";
}

std::string optional_double_to_csv_cell(const std::optional<double>& value)
{
    if (!value) {
        return "";
    }

    std::ostringstream oss;
    oss << std::setprecision(17) << *value;
    return oss.str();
}

std::string distances_to_csv_cell(const std::vector<double>& distances)
{
    std::ostringstream oss;
    oss << std::setprecision(17);
    for (std::size_t i = 0; i < distances.size(); ++i) {
        if (i > 0) {
            oss << ';';
        }
        oss << distances[i];
    }
    return oss.str();
}

void write_statistics_csv_header(std::ofstream& csv)
{
    csv << "relative_alpha,alpha,relative_offset,offset,tau,absolute_tau,runtime_s,"
           "total_output_face_count,total_output_vertex_count,valid_binary,"
           "directed_chamfer_distance,directed_hausdorff_distance,"
           "sampled_output_to_input_distances\n";
}

void append_statistics_csv_row(std::ofstream& csv, const Statistics_result& result)
{
    csv << result.relative_alpha << ','
        << result.alpha << ','
        << result.relative_offset << ','
        << result.offset << ','
        << result.tau << ','
        << result.absolute_tau << ','
        << result.runtime << ','
        << result.total_output_face_count << ','
        << result.total_output_vertex_count << ','
        << optional_int_to_csv_cell(result.valid_binary) << ','
        << optional_double_to_csv_cell(result.directed_chamfer_distance) << ','
        << optional_double_to_csv_cell(result.directed_hausdorff_distance) << ",\""
        << distances_to_csv_cell(result.sampled_output_to_input_distances) << "\"\n";
}

std::ofstream open_statistics_csv_or_throw(const std::string& csv_output_path)
{
    std::filesystem::path out_path(csv_output_path);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    std::ofstream csv(csv_output_path);
    if (!csv) {
        throw std::runtime_error("Could not open CSV output path: " + csv_output_path);
    }

    return csv;
}
} // namespace

std::vector<Point_3> random_surface_samples_on_mesh(
    const Mesh& mesh,
    std::size_t number_of_points,
    bool sample_vertices,
    bool sample_edges)
{
    ensure_triangle_mesh_for_distance(mesh, "mesh");

    const unsigned int sample_count = to_unsigned_sample_count(number_of_points);
    if (sample_count == 0) {
        return {};
    }

    std::vector<Point_3> samples;
    samples.reserve(static_cast<std::size_t>(sample_count));

    PMP::sample_triangle_mesh(
        mesh,
        std::back_inserter(samples),
        CGAL::parameters::use_random_uniform_sampling(true)
            .number_of_points_on_faces(sample_count)
            .do_sample_faces(true)
            .do_sample_vertices(sample_vertices)
            .do_sample_edges(sample_edges));

    return samples;
}

double directed_chamfer_distance(
    const Mesh& source_mesh,
    const Mesh& target_mesh,
    std::size_t number_of_samples)
{
    ensure_triangle_mesh_for_distance(source_mesh, "source_mesh");
    ensure_triangle_mesh_for_distance(target_mesh, "target_mesh");

    const std::vector<Point_3> samples =
        random_surface_samples_on_mesh(source_mesh, number_of_samples, false, false);
    if (samples.empty()) {
        return 0.0;
    }

    const std::vector<double> distances = distances_to_mesh(samples, target_mesh);
    return mean_of_distances(distances);
}

double directed_hausdorff_distance(
    const Mesh& source_mesh,
    const Mesh& target_mesh,
    std::size_t number_of_samples)
{
    ensure_triangle_mesh_for_distance(source_mesh, "source_mesh");
    ensure_triangle_mesh_for_distance(target_mesh, "target_mesh");

    const unsigned int sample_count = to_unsigned_sample_count(number_of_samples);
    if (sample_count == 0) {
        return 0.0;
    }

    const std::vector<Point_3> samples =
        random_surface_samples_on_mesh(source_mesh, sample_count, false, false);
    const std::vector<double> distances = distances_to_mesh(samples, target_mesh);
    return max_of_distances(distances);
}

bool run_val3dity_test(const Mesh& mesh)
{
    std::filesystem::create_directories("../data/Output/val3dity_check");
    return valid_mesh_boolean(mesh);
}

bool run_val3dity_test(const std::string& input_path)
{
    std::filesystem::create_directories("../data/Output/val3dity_check");
    return valid_file_boolean(input_path);
}

std::size_t output_mesh_vertex_count(const Mesh& mesh)
{
    return num_vertices(mesh);
}

std::size_t output_mesh_face_count(const Mesh& mesh)
{
    return num_faces(mesh);
}

Statistics_result statisctics(
    const double relative_alpha,
    const double relative_offset,
    const double tau,
    Mesh input,
    const bool validate,
    const bool use_beeren_method,
    const bool statistics,
    const bool write_output)
{
    ensure_positive_finite(relative_alpha, "relative_alpha");
    ensure_positive_finite(relative_offset, "relative_offset");
    ensure_positive_finite(tau, "tau");
    if (use_beeren_method && (!std::isfinite(tau) || tau <= 1.0)) {
        throw std::invalid_argument("tau must be finite and > 1.0 when use_beeren_method=true");
    }

    ensure_triangle_mesh_in_place(input, "input");

    Statistics_result result;
    result.relative_alpha = relative_alpha;
    result.relative_offset = relative_offset;
    result.tau = tau;

    const double diagonal = bbox_diagonal(input);
    result.alpha = diagonal / relative_alpha;
    result.offset = diagonal / relative_offset;
    result.absolute_tau = tau * result.offset;

    Mesh wrap;
    CGAL::Real_timer timer;
    timer.start();
    if (use_beeren_method) {
        CGAL::alpha_wrap_3(
            input,
            result.alpha,
            result.offset,
            wrap,
            CGAL::parameters::max_distance_to_input_in_offsets(tau));
    } else {
        CGAL::alpha_wrap_3(input, result.alpha, result.offset, wrap);
    }
    timer.stop();
    result.runtime = timer.time();

    result.total_output_face_count = num_faces(wrap);
    result.total_output_vertex_count = num_vertices(wrap);

    if (validate) {
        result.valid_binary = validate_per_connected_component_binary(wrap);
    }

    if (statistics) {
        const std::vector<Point_3> sampled_output_points =
            random_surface_samples_on_mesh(wrap, k_default_statistics_sample_count, false, false);

        result.sampled_output_to_input_distances = distances_to_mesh(sampled_output_points, input);
        result.directed_chamfer_distance = mean_of_distances(result.sampled_output_to_input_distances);
        result.directed_hausdorff_distance = max_of_distances(result.sampled_output_to_input_distances);
    }

    if (write_output) {
        std::string output_path =
            statistics_output_name(relative_alpha, relative_offset, tau, use_beeren_method);
        std::filesystem::path output_fs_path(output_path);
        std::filesystem::create_directories(output_fs_path.parent_path());
        CGAL::IO::write_polygon_mesh(output_path, wrap, CGAL::parameters::stream_precision(25));
    }

    return result;
}

Statistics_result statistics(
    const double relative_alpha,
    const double relative_offset,
    const double tau,
    Mesh input,
    const bool validate,
    const bool use_beeren_method,
    const bool compute_statistics,
    const bool write_output)
{
    return statisctics(
        relative_alpha,
        relative_offset,
        tau,
        std::move(input),
        validate,
        use_beeren_method,
        compute_statistics,
        write_output);
}

void statistics_over_relative_alpha_to_csv(
    const std::vector<double>& relative_alpha_values,
    const double relative_offset,
    const double tau,
    const Mesh& input,
    const bool validate,
    const bool use_beeren_method,
    const bool compute_statistics,
    const bool write_output,
    const std::string& csv_output_path)
{
    if (relative_alpha_values.empty()) {
        throw std::invalid_argument("relative_alpha_values is empty");
    }

    std::ofstream csv = open_statistics_csv_or_throw(csv_output_path);
    write_statistics_csv_header(csv);

    for (const double relative_alpha : relative_alpha_values) {
        const Statistics_result result = statisctics(
            relative_alpha,
            relative_offset,
            tau,
            input,
            validate,
            use_beeren_method,
            compute_statistics,
            write_output);
        append_statistics_csv_row(csv, result);
    }
}

void statistics_over_relative_offset_to_csv(
    const double relative_alpha,
    const std::vector<double>& relative_offset_values,
    const double tau,
    const Mesh& input,
    const bool validate,
    const bool use_beeren_method,
    const bool compute_statistics,
    const bool write_output,
    const std::string& csv_output_path)
{
    if (relative_offset_values.empty()) {
        throw std::invalid_argument("relative_offset_values is empty");
    }

    std::ofstream csv = open_statistics_csv_or_throw(csv_output_path);
    write_statistics_csv_header(csv);

    for (const double relative_offset : relative_offset_values) {
        const Statistics_result result = statisctics(
            relative_alpha,
            relative_offset,
            tau,
            input,
            validate,
            use_beeren_method,
            compute_statistics,
            write_output);
        append_statistics_csv_row(csv, result);
    }
}

void statistics_over_tau_to_csv(
    const double relative_alpha,
    const double relative_offset,
    const std::vector<double>& tau_values,
    const Mesh& input,
    const bool validate,
    const bool use_beeren_method,
    const bool compute_statistics,
    const bool write_output,
    const std::string& csv_output_path)
{
    if (tau_values.empty()) {
        throw std::invalid_argument("tau_values is empty");
    }

    std::ofstream csv = open_statistics_csv_or_throw(csv_output_path);
    write_statistics_csv_header(csv);

    for (const double tau : tau_values) {
        const Statistics_result result = statisctics(
            relative_alpha,
            relative_offset,
            tau,
            input,
            validate,
            use_beeren_method,
            compute_statistics,
            write_output);
        append_statistics_csv_row(csv, result);
    }
}
