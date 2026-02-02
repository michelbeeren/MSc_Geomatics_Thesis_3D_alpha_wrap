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
    // std::cout << "diagonal bbox length: " << diag_length << std::endl;
    const double alpha = diag_length / relative_alpha_;
    const double density = 10.0 / (std::sqrt(3.0) * alpha * alpha);

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
                        .do_sample_edges(false)
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

void write_points_as_npy(const std::string& filename, const std::vector<Point_3>& points)
{
    const size_t N = points.size();
    if (N == 0) throw std::runtime_error("write_points_as_npy: no points");

    std::vector<double> data;
    data.reserve(N * 3);

    for (const auto& p : points) {
        data.push_back(CGAL::to_double(p.x()));
        data.push_back(CGAL::to_double(p.y()));
        data.push_back(CGAL::to_double(p.z()));
    }

    // cnpy in masbcpp expects unsigned int* + ndims
    if (N > std::numeric_limits<unsigned int>::max())
        throw std::runtime_error("Too many points for cnpy shape (unsigned int overflow)");

    unsigned int shape[2] = {
        static_cast<unsigned int>(N),
        3u
    };

    cnpy::npy_save<double>(filename, data.data(), shape, 2, "w");
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

void write_npy_points_as_off(const std::string& npy_file, const std::string& off_file)
{
    cnpy::NpyArray arr = cnpy::npy_load(npy_file);

    if (arr.shape.size() != 2 || arr.shape[1] != 3)
        throw std::runtime_error("Expected Nx3 array");

    if (arr.word_size != sizeof(float))
        throw std::runtime_error("Expected float32 array");

    const size_t N = arr.shape[0];
    const float* data = reinterpret_cast<const float*>(arr.data);

    std::ofstream out(off_file);
    if (!out)
        throw std::runtime_error("Cannot open " + off_file);

    out << std::setprecision(17) << std::fixed;
    out << "OFF\n";
    out << N << " 0 0\n";

    for (size_t i = 0; i < N; ++i) {
        out << data[3*i + 0] << " "
            << data[3*i + 1] << " "
            << data[3*i + 2] << "\n";
    }
}

#include <climits>
#include <cstdint>
#include <limits>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <algorithm>

void write_mat_colored_coff(const std::string& ma_coords_out_npy,const std::string& ma_qidx_out_npy,const std::string& coords_indexed_by_qidx_npy,const std::string& coff_file)
{
    cnpy::NpyArray ma_arr   = cnpy::npy_load(ma_coords_out_npy);
    cnpy::NpyArray qidx_arr = cnpy::npy_load(ma_qidx_out_npy);
    cnpy::NpyArray in_arr   = cnpy::npy_load(coords_indexed_by_qidx_npy);

    // --- checks: coords arrays must be Nx3 float32 ---
    auto check_coords = [&](const cnpy::NpyArray& a, const std::string& name) {
        if (a.shape.size() != 2 || a.shape[1] != 3)
            throw std::runtime_error(name + " must be Nx3");
        if (a.word_size != sizeof(float))
            throw std::runtime_error(name + " must be float32");
    };

    check_coords(ma_arr, "ma_coords_out");
    check_coords(in_arr, "coords_indexed_by_qidx");

    const size_t Nma = ma_arr.shape[0];
    const size_t Nin = in_arr.shape[0];

    // --- qidx shape: allow 1D (N), or 2D (N x 1/2) ---
    size_t qcols = 0;
    if (qidx_arr.shape.size() == 1) {
        if (qidx_arr.shape[0] != Nma)
            throw std::runtime_error("ma_qidx_out length must match ma_coords_out rows");
        qcols = 1;
    } else if (qidx_arr.shape.size() == 2) {
        if (qidx_arr.shape[0] != Nma)
            throw std::runtime_error("ma_qidx_out rows must match ma_coords_out rows");
        qcols = qidx_arr.shape[1];
        if (qcols != 1 && qcols != 2)
            throw std::runtime_error("ma_qidx_out must have 1 or 2 columns");
    } else {
        throw std::runtime_error("ma_qidx_out must be 1D or 2D");
    }

    const float* ma = reinterpret_cast<const float*>(ma_arr.data);
    const float* in = reinterpret_cast<const float*>(in_arr.data);

    // Read a flattened qidx element as signed long long, supporting int/uint 32/64.
    auto q_read_as_ll = [&](size_t flat_idx) -> long long {
        const char* raw = reinterpret_cast<const char*>(qidx_arr.data);

        if (qidx_arr.word_size == 4) {
            // Try to interpret as signed first, but also allow unsigned values.
            const int32_t* qs = reinterpret_cast<const int32_t*>(raw);
            int32_t v = qs[flat_idx];
            // If negative, keep it (sentinel like -1).
            if (v < 0) return (long long)v;

            // Otherwise, it could still be uint32, but reading as int32 is fine for <2^31.
            // If your values exceed 2^31-1, you likely want uint32:
            const uint32_t* qu = reinterpret_cast<const uint32_t*>(raw);
            return (long long)qu[flat_idx];
        }

        if (qidx_arr.word_size == 8) {
            const int64_t* qs = reinterpret_cast<const int64_t*>(raw);
            int64_t v = qs[flat_idx];
            if (v < 0) return (long long)v;

            const uint64_t* qu = reinterpret_cast<const uint64_t*>(raw);
            // Clamp huge unsigned values into signed range (should not happen for indices)
            if (qu[flat_idx] > (uint64_t)LLONG_MAX) return LLONG_MAX;
            return (long long)qu[flat_idx];
        }

        throw std::runtime_error("Unsupported qidx dtype (expected 4 or 8 bytes)");
    };

    // Helper: access qidx at (row i, col c)
    auto q_at = [&](size_t i, size_t c) -> long long {
        const size_t flat_idx = (qidx_arr.shape.size() == 1) ? i : (i * qcols + c);
        return q_read_as_ll(flat_idx);
    };

    // --- compute qmin/qmax for debugging + 1-based detection ---
    long long qmin = LLONG_MAX, qmax = LLONG_MIN;
    const size_t qN = (qidx_arr.shape.size() == 1) ? qidx_arr.shape[0] : qidx_arr.shape[0] * qidx_arr.shape[1];

    for (size_t i = 0; i < qN; ++i) {
        long long v = q_read_as_ll(i);
        qmin = std::min(qmin, v);
        qmax = std::max(qmax, v);
    }

    std::cout << "qidx stats: Nin=" << Nin
              << " qmin=" << qmin
              << " qmax=" << qmax
              << " qcols=" << qcols
              << " qidx_word_size=" << qidx_arr.word_size
              << " qidx_dims=" << qidx_arr.shape.size()
              << "\n";

    // Detect 1-based indexing in the *valid* index range.
    // Typical: qmin>=1 and qmax==Nin (or <=Nin).
    const bool one_based = (qmin >= 1 && qmax <= (long long)Nin && qmax > 0);

    auto to_zero_based = [&](long long qi) -> long long {
        if (one_based) return qi - 1;
        return qi;
    };

    // --- radius proxy per MA point: distance to its (one or two) indexed input points ---
    std::vector<double> r;
    r.reserve(Nma);

    double rmin = std::numeric_limits<double>::infinity();
    double rmax = 0.0;

    auto dist_to = [&](double mx, double my, double mz, long long qi_raw) -> double {
        long long qi = to_zero_based(qi_raw);

        // allow sentinels like -1
        if (qi < 0) return std::numeric_limits<double>::infinity();
        if ((size_t)qi >= Nin) return std::numeric_limits<double>::infinity();

        const size_t j = (size_t)qi;
        const double sx = in[3*j + 0];
        const double sy = in[3*j + 1];
        const double sz = in[3*j + 2];
        const double dx = mx - sx, dy = my - sy, dz = mz - sz;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    };

    for (size_t i = 0; i < Nma; ++i) {
        const double mx = ma[3*i + 0];
        const double my = ma[3*i + 1];
        const double mz = ma[3*i + 2];

        const long long q0 = q_at(i, 0);
        double d = dist_to(mx, my, mz, q0);

        if (qcols == 2) {
            const long long q1 = q_at(i, 1);
            const double d1 = dist_to(mx, my, mz, q1);
            d = std::min(d, d1);
        }

        // If both indices are invalid, choose a default
        if (!std::isfinite(d)) {
            // You can skip instead, but skipping would desync Nma vs written count.
            // So we give a neutral radius:
            d = 0.0;
        }

        r.push_back(d);
        rmin = std::min(rmin, d);
        rmax = std::max(rmax, d);
    }

    // --- write COFF ---
    std::filesystem::create_directories(std::filesystem::path(coff_file).parent_path());

    std::ofstream out(coff_file);
    if (!out) throw std::runtime_error("Cannot open " + coff_file);

    out << std::setprecision(17) << std::fixed;
    out << "COFF\n";
    out << Nma << " 0 0\n";

    const double denom = (rmax - rmin) + 1e-12;

    for (size_t i = 0; i < Nma; ++i) {
        const double t = (r[i] - rmin) / denom;  // 0..1

        // small radius => red, large => blue
        const int R = (int)std::lround(255.0 * (1.0 - t));
        const int G = 0;
        const int B = (int)std::lround(255.0 * t);
        const int A = 255;

        out << ma[3*i + 0] << " "
            << ma[3*i + 1] << " "
            << ma[3*i + 2] << " "
            << R << " " << G << " " << B << " " << A << "\n";
    }

    std::cout << "Wrote " << Nma << " colored MA points to " << coff_file
              << " (one_based=" << (one_based ? "true" : "false") << ")\n";
}



