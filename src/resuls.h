//
// Created by Michel Beeren on 13/05/2026.
//

#ifndef THESIS_RESULS_H
#define THESIS_RESULS_H

#include <cstddef>
#include <optional>
#include <string>
#include <vector>
#include <utility>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Surface_mesh.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;

struct Statistics_result
{
    double relative_alpha = 0.0;
    double alpha = 0.0;
    double relative_offset = 0.0;
    double offset = 0.0;
    double tau = 0.0;
    double absolute_tau = 0.0;
    double upper_bound = 0.0;
    double runtime = 0.0;
    std::size_t total_output_face_count = 0;
    std::size_t total_output_vertex_count = 0;
    std::optional<int> valid_binary; // 1 = valid, 0 = invalid
    std::optional<double> directed_chamfer_distance;
    std::optional<double> directed_hausdorff_distance;
    std::vector<double> sampled_output_to_input_distances;
};

std::vector<Point_3> random_surface_samples_on_mesh(
    const Mesh& mesh,
    std::size_t number_of_points,
    bool sample_vertices = false,
    bool sample_edges = false);

double directed_chamfer_distance(
    const Mesh& source_mesh,
    const Mesh& target_mesh,
    std::size_t number_of_samples = 10000);

double directed_hausdorff_distance(
    const Mesh& source_mesh,
    const Mesh& target_mesh,
    std::size_t number_of_samples = 10000);

bool run_val3dity_test(const Mesh& mesh);
bool run_val3dity_test(const std::string& input_path);

std::size_t output_mesh_vertex_count(const Mesh& mesh);
std::size_t output_mesh_face_count(const Mesh& mesh);

Mesh exploder(
    const std::string& input_off_path,
    double mean_shift_distance,
    double stddev_shift_distance);

bool exploder_to_off(
    const std::string& input_off_path,
    const std::string& output_off_path,
    double mean_shift_distance,
    double stddev_shift_distance);

Statistics_result statisctics(
    double relative_alpha,
    double relative_offset,
    double tau,
    Mesh input,
    bool validate,
    bool use_beeren_method,
    bool statistics,
    bool write_output,
    std::size_t wrapping_time_repetitions = 10);

Statistics_result statistics(
    double relative_alpha,
    double relative_offset,
    double tau,
    Mesh input,
    bool validate,
    bool use_beeren_method,
    bool compute_statistics,
    bool write_output,
    std::size_t wrapping_time_repetitions = 10);

void statistics_over_relative_alpha_to_csv(
    const std::vector<double>& relative_alpha_values,
    double relative_offset,
    double tau,
    const Mesh& input,
    bool validate,
    bool use_beeren_method,
    bool compute_statistics,
    bool write_output,
    const std::string& csv_output_path,
    std::size_t wrapping_time_repetitions = 10);

void statistics_over_relative_offset_to_csv(
    double relative_alpha,
    const std::vector<double>& relative_offset_values,
    double tau,
    const Mesh& input,
    bool validate,
    bool use_beeren_method,
    bool compute_statistics,
    bool write_output,
    const std::string& csv_output_path,
    std::size_t wrapping_time_repetitions = 10);

void statistics_over_tau_to_csv(
    double relative_alpha,
    double relative_offset,
    const std::vector<double>& tau_values,
    const Mesh& input,
    bool validate,
    bool use_beeren_method,
    bool compute_statistics,
    bool write_output,
    const std::string& csv_output_path,
    std::size_t wrapping_time_repetitions = 10);

template <typename Callable>
double runtime_seconds(Callable&& callable)
{
    CGAL::Real_timer timer;
    timer.start();
    std::forward<Callable>(callable)();
    timer.stop();
    return timer.time();
}

#endif //THESIS_RESULS_H
