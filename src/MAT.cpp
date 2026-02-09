//
// Created by Michel Beeren on 02/02/2026.
//

#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/measure.h>   // area()
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <array>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include "cnpy.h"
#include <climits>
#include <cstdint>

namespace PMP = CGAL::Polygon_mesh_processing;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
using Primitive   = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using AABB_traits = CGAL::AABB_traits<K, Primitive>;
using Tree        = CGAL::AABB_tree<AABB_traits>;

#include "MAT.h"

// Blue noise sample input mesh
std::vector<Point_3> _surface_sampling(const Mesh& mesh, const double relative_alpha_) {
    std::vector<Point_3> surface_points;

    // compute alpha and offset from a_rel and d_rel and bbox
    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                         CGAL::square(bbox.ymax() - bbox.ymin()) +
                                         CGAL::square(bbox.zmax() - bbox.zmin()));
    std::cout << "diagonal bbox length: " << diag_length << std::endl;
    const double alpha = diag_length / relative_alpha_;
    const double density = 300. * std::pow((1-alpha),10);

    // ✅ ADD THIS BLOCK RIGHT HERE
    std::cout
        << "tri=" << CGAL::is_triangle_mesh(mesh)
        << " faces=" << num_faces(mesh)
        << " area=" << PMP::area(mesh)
        << " alpha=" << alpha
        << " density=" << density
        << " expected_points~" << PMP::area(mesh) * density
        << std::endl;

    // (optional but very useful safety checks)
    if (!std::isfinite(alpha) || alpha <= 0.0)
        throw std::runtime_error("alpha invalid");
    if (!std::isfinite(density) || density <= 0.0)
        throw std::runtime_error("density invalid");

    PMP::sample_triangle_mesh(
        mesh,
        std::back_inserter(surface_points),
        CGAL::parameters::use_random_uniform_sampling(true)
                        .number_of_points_per_area_unit(density)
                        .do_sample_vertices(false)
                        .do_sample_edges(true)
                        .do_sample_faces(true)
    );

    std::cout << "Sampled " << surface_points.size() << std::endl;

    return surface_points;
}

// For visualization of the blue noise samling
void write_points_as_off(const std::string& filename, const std::vector<Point_3>& points)
{
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Cannot open file " << filename << "\n";
        return;
    }

    out << std::setprecision(17) << std::fixed; // <<<<<< IMPORTANT

    out << "OFF\n";
    out << points.size() << " 0 0\n";

    for (const auto& p : points) {
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
}

std::vector<Vector_3> normals_for_points_from_closest_face(const Mesh& mesh, const Tree& tree, const std::map<Mesh::Face_index, Vector_3>& face_normals, const std::vector<Point_3>& points)
{
    std::vector<Vector_3> normals;
    normals.reserve(points.size());

    for (const auto& p : points) {
        Mesh::Face_index f = tree.closest_point_and_primitive(p).second;

        auto it = face_normals.find(f);
        if (it == face_normals.end()) {
            throw std::runtime_error("Face normal not found for closest face");
        }

        Vector_3 n = it->second;

        // normalize for safety
        const double len2 = CGAL::to_double(n.squared_length());
        if (len2 > 0.0) n = n / std::sqrt(len2);

        normals.push_back(n);
    }

    return normals;
}

void write_coords_npy(const std::string& coords_path, const std::vector<Point_3>& points)
{
    const size_t N = points.size();
    if (N == 0) throw std::runtime_error("No points to write");

    std::filesystem::path p(coords_path);
    std::filesystem::create_directories(p.parent_path());

    std::vector<float> data;
    data.reserve(N * 3);

    for (const auto& pt : points) {
        data.push_back((float)CGAL::to_double(pt.x()));
        data.push_back((float)CGAL::to_double(pt.y()));
        data.push_back((float)CGAL::to_double(pt.z()));
    }

    if (N > std::numeric_limits<unsigned int>::max())
        throw std::runtime_error("Too many points for cnpy unsigned int shape");

    unsigned int shape[2] = { (unsigned int)N, 3u };
    cnpy::npy_save<float>(coords_path, data.data(), shape, 2, "w");
}

void write_normals_npy(const std::string& normals_path, const std::vector<Vector_3>& normals)
{
    const size_t N = normals.size();
    if (N == 0) throw std::runtime_error("No normals to write");

    std::filesystem::path p(normals_path);
    std::filesystem::create_directories(p.parent_path());

    std::vector<float> data;
    data.reserve(N * 3);

    for (const auto& n : normals) {
        data.push_back((float)CGAL::to_double(n.x()));
        data.push_back((float)CGAL::to_double(n.y()));
        data.push_back((float)CGAL::to_double(n.z()));
    }

    if (N > std::numeric_limits<unsigned int>::max())
        throw std::runtime_error("Too many normals for cnpy unsigned int shape");

    unsigned int shape[2] = { (unsigned int)N, 3u };
    cnpy::npy_save<float>(normals_path, data.data(), shape, 2, "w");
}

void run_masbcpp_compute_ma(const std::string& in_dir, const std::string& out_dir)
{
    std::filesystem::create_directories(out_dir);

    std::string cmd =
        std::string("\"") + MASBCPP_COMPUTE_MA_PATH + "\" \"" + in_dir + "\" \"" + out_dir + "\"";

    // Redirect stderr->stdout so you see the real error in CLion
    cmd += " 2>&1";

    std::cout << "Running: " << cmd << "\n";

    int rc = std::system(cmd.c_str());
    std::cout << "compute_ma rc=" << rc << "\n";

    if (rc != 0) throw std::runtime_error("masbcpp compute_ma failed");
}

static void require_shape_N3(const cnpy::NpyArray& a, const std::string& name)
{
    if (a.shape.size() != 2 || a.shape[1] != 3)
        throw std::runtime_error(name + " must have shape (N, 3)");
}

static void require_shape_N(const cnpy::NpyArray& a, size_t N, const std::string& name)
{
    if (a.shape.size() != 1 || a.shape[0] != N)
        throw std::runtime_error(name + " must have shape (N,)");
}

// Convert qidx array (int32 or int64) into a uniform int64_t vector
static std::vector<int64_t> load_qidx_as_i64(const cnpy::NpyArray& arr, const std::string& name)
{
    std::vector<int64_t> out;
    out.resize(arr.shape[0]);

    if (arr.word_size == sizeof(int64_t))
    {
        const int64_t* p = reinterpret_cast<const int64_t*>(arr.data);
        std::copy(p, p + arr.shape[0], out.begin());
    }
    else if (arr.word_size == sizeof(int32_t))
    {
        const int32_t* p = reinterpret_cast<const int32_t*>(arr.data);
        for (size_t i = 0; i < arr.shape[0]; ++i) out[i] = static_cast<int64_t>(p[i]);
    }
    else
    {
        throw std::runtime_error(name + " has unsupported dtype: word_size=" + std::to_string(arr.word_size) +
                                 " (expected int32 or int64)");
    }
    return out;
}

void generate_lfs_ply(const std::string& in_dir,
                      const std::string& out_dir,
                      const std::string& output_path)
{
    // 1) Load NPY files
    cnpy::NpyArray arr_coords = cnpy::npy_load(in_dir  + "/coords.npy");
    cnpy::NpyArray arr_ma_in  = cnpy::npy_load(out_dir + "/ma_coords_in.npy");
    cnpy::NpyArray arr_ma_out = cnpy::npy_load(out_dir + "/ma_coords_out.npy");
    cnpy::NpyArray arr_qi     = cnpy::npy_load(out_dir + "/ma_qidx_in.npy");
    cnpy::NpyArray arr_qo     = cnpy::npy_load(out_dir + "/ma_qidx_out.npy");

    // 2) Validate shapes
    require_shape_N3(arr_coords, "coords");
    require_shape_N3(arr_ma_in,  "ma_coords_in");
    require_shape_N3(arr_ma_out, "ma_coords_out");

    const size_t N = arr_coords.shape[0];
    require_shape_N(arr_qi, N, "ma_qidx_in");
    require_shape_N(arr_qo, N, "ma_qidx_out");

    const size_t Nin_ma  = arr_ma_in.shape[0];
    const size_t Nout_ma = arr_ma_out.shape[0];

    // 3) Expect float32 for coords + MA
    if (arr_coords.word_size != sizeof(float))
        throw std::runtime_error("coords.npy is not float32 (word_size=" + std::to_string(arr_coords.word_size) + ")");
    if (arr_ma_in.word_size != sizeof(float) || arr_ma_out.word_size != sizeof(float))
        throw std::runtime_error("ma_coords_*.npy are not float32");

    // 4) Access float pointers
    const float* coords = reinterpret_cast<const float*>(arr_coords.data);
    const float* ma_in  = reinterpret_cast<const float*>(arr_ma_in.data);
    const float* ma_out = reinterpret_cast<const float*>(arr_ma_out.data);

    // 5) Load qidx as int64 (supports int32 or int64 input)
    std::vector<int64_t> qidx_in  = load_qidx_as_i64(arr_qi, "ma_qidx_in");
    std::vector<int64_t> qidx_out = load_qidx_as_i64(arr_qo, "ma_qidx_out");

    // 6) Write PLY header
    std::ofstream out(output_path);
    if (!out.is_open())
        throw std::runtime_error("Could not open output PLY file: " + output_path);

    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << N << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "property float feature_size\n";
    out << "end_header\n";
    out << std::fixed << std::setprecision(7);

    // 7) Compute LFS with -1 + bounds checks
    for (size_t i = 0; i < N; ++i)
    {
        const float px = coords[i*3 + 0];
        const float py = coords[i*3 + 1];
        const float pz = coords[i*3 + 2];

        float dist_in  = std::numeric_limits<float>::infinity();
        float dist_out = std::numeric_limits<float>::infinity();

        const int64_t ii = qidx_in[i];
        if (ii >= 0 && static_cast<size_t>(ii) < Nin_ma)
        {
            const float dx = px - ma_in[ii*3 + 0];
            const float dy = py - ma_in[ii*3 + 1];
            const float dz = pz - ma_in[ii*3 + 2];
            dist_in = std::sqrt(dx*dx + dy*dy + dz*dz);
        }

        const int64_t io = qidx_out[i];
        if (io >= 0 && static_cast<size_t>(io) < Nout_ma)
        {
            const float dx = px - ma_out[io*3 + 0];
            const float dy = py - ma_out[io*3 + 1];
            const float dz = pz - ma_out[io*3 + 2];
            dist_out = std::sqrt(dx*dx + dy*dy + dz*dz);
        }

        float lfs = std::min(dist_in, dist_out);
        if (!std::isfinite(lfs)) lfs = -1.0f;

        out << px << " " << py << " " << pz << " " << lfs << "\n";
    }

    out.close();
}

void generate_lfs_ply_concave(const std::string& in_dir,
                      const std::string& out_dir,
                      const std::string& output_path)
{
    // 1) Load NPY files
    cnpy::NpyArray arr_coords = cnpy::npy_load(in_dir  + "/coords.npy");
    cnpy::NpyArray arr_ma_out = cnpy::npy_load(out_dir + "/ma_coords_out.npy");
    cnpy::NpyArray arr_qo     = cnpy::npy_load(out_dir + "/ma_qidx_out.npy");

    // 2) Validate shapes
    require_shape_N3(arr_coords, "coords");
    require_shape_N3(arr_ma_out, "ma_coords_out");

    const size_t N = arr_coords.shape[0];
    require_shape_N(arr_qo, N, "ma_qidx_out");

    const size_t Nout_ma = arr_ma_out.shape[0];

    // 3) Expect float32 for coords + MA
    if (arr_coords.word_size != sizeof(float))
        throw std::runtime_error("coords.npy is not float32 (word_size=" + std::to_string(arr_coords.word_size) + ")");
    if (arr_ma_out.word_size != sizeof(float))
        throw std::runtime_error("ma_coords_out.npy is not float32");

    // 4) Access float pointers
    const float* coords = reinterpret_cast<const float*>(arr_coords.data);
    const float* ma_out = reinterpret_cast<const float*>(arr_ma_out.data);

    // 5) Load qidx as int64 (supports int32 or int64 input)
    std::vector<int64_t> qidx_out = load_qidx_as_i64(arr_qo, "ma_qidx_out");

    // 6) Write PLY header
    std::ofstream out(output_path);
    if (!out.is_open())
        throw std::runtime_error("Could not open output PLY file: " + output_path);

    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << N << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "property float feature_size\n";
    out << "end_header\n";
    out << std::fixed << std::setprecision(7);

    // 7) Compute LFS using only the outer medial axis (ignore the inner)
    for (size_t i = 0; i < N; ++i)
    {
        const float px = coords[i*3 + 0];
        const float py = coords[i*3 + 1];
        const float pz = coords[i*3 + 2];

        float dist_out = std::numeric_limits<float>::infinity();

        const int64_t io = qidx_out[i];
        if (io >= 0 && static_cast<size_t>(io) < Nout_ma)
        {
            const float dx = px - ma_out[io*3 + 0];
            const float dy = py - ma_out[io*3 + 1];
            const float dz = pz - ma_out[io*3 + 2];
            dist_out = std::sqrt(dx*dx + dy*dy + dz*dz);
        }

        // Since we're using only the outer MA, dist_in is not considered (set to inf)
        float lfs = dist_out;
        if (!std::isfinite(lfs)) lfs = -1.0f;

        out << px << " " << py << " " << pz << " " << lfs << "\n";
    }

    out.close();
}