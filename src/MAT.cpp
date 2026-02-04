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

void write_coords_npy2(const std::string& coords_path, const std::vector<Point_3>& points)
{
    const size_t N = points.size();
    if (N == 0) throw std::runtime_error("No points to write");

    std::filesystem::path p(coords_path);
    std::filesystem::create_directories(p.parent_path());

    // --- CHANGE: Use double instead of float ---
    std::vector<double> data;
    data.reserve(N * 3);

    for (const auto& pt : points) {
        // Remove the (float) cast to keep full precision
        data.push_back(CGAL::to_double(pt.x()));
        data.push_back(CGAL::to_double(pt.y()));
        data.push_back(CGAL::to_double(pt.z()));
    }

    unsigned int shape[2] = { (unsigned int)N, 3u };
    // --- CHANGE: Save as double ---
    cnpy::npy_save<double>(coords_path, data.data(), shape, 2, "w");
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

void write_normals_npy2(const std::string& normals_path, const std::vector<Vector_3>& normals)
{
    const size_t N = normals.size();
    if (N == 0) throw std::runtime_error("No normals to write");

    std::filesystem::path p(normals_path);
    std::filesystem::create_directories(p.parent_path());

    // Switch to double to preserve CGAL precision
    std::vector<double> data;
    data.reserve(N * 3);

    for (const auto& n : normals) {
        data.push_back(CGAL::to_double(n.x()));
        data.push_back(CGAL::to_double(n.y()));
        data.push_back(CGAL::to_double(n.z()));
    }

    unsigned int shape[2] = { (unsigned int)N, 3u };
    // Save as double
    cnpy::npy_save<double>(normals_path, data.data(), shape, 2, "w");
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

    size_t bad = 0, neg = 0;
    for (size_t i = 0; i < Nma; ++i) {
        for (size_t c = 0; c < qcols; ++c) {
            long long qi = to_zero_based(q_at(i,c));
            if (qi < 0) neg++;
            else if ((size_t)qi >= Nin) bad++;
        }
    }
    std::cout << "qidx: neg=" << neg << " out_of_range=" << bad << "\n";

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

// Color surface points (coords) by proxy ball radius inferred from MA point(s) that reference them.
void write_surface_colored_by_ball_proxy(const std::string& ma_coords_out_npy, const std::string& ma_qidx_out_npy, const std::string& coords_indexed_by_qidx_npy, const std::string& coff_file)
{
    cnpy::NpyArray ma_arr   = cnpy::npy_load(ma_coords_out_npy);
    cnpy::NpyArray qidx_arr = cnpy::npy_load(ma_qidx_out_npy);
    cnpy::NpyArray in_arr   = cnpy::npy_load(coords_indexed_by_qidx_npy);

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

    // qidx shape: allow 1D (N), or 2D (N x 1/2)
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
            const int32_t* qs = reinterpret_cast<const int32_t*>(raw);
            int32_t v = qs[flat_idx];
            if (v < 0) return (long long)v; // sentinel like -1
            const uint32_t* qu = reinterpret_cast<const uint32_t*>(raw);
            return (long long)qu[flat_idx];
        }

        if (qidx_arr.word_size == 8) {
            const int64_t* qs = reinterpret_cast<const int64_t*>(raw);
            int64_t v = qs[flat_idx];
            if (v < 0) return (long long)v;
            const uint64_t* qu = reinterpret_cast<const uint64_t*>(raw);
            if (qu[flat_idx] > (uint64_t)LLONG_MAX) return LLONG_MAX;
            return (long long)qu[flat_idx];
        }

        throw std::runtime_error("Unsupported qidx dtype (expected 4 or 8 bytes)");
    };

    auto q_at = [&](size_t i, size_t c) -> long long {
        const size_t flat_idx = (qidx_arr.shape.size() == 1) ? i : (i * qcols + c);
        return q_read_as_ll(flat_idx);
    };

    // --- compute qmin/qmax for debugging + detect 1-based ---
    long long qmin = LLONG_MAX, qmax = LLONG_MIN;
    const size_t qN = (qidx_arr.shape.size() == 1) ? qidx_arr.shape[0] : qidx_arr.shape[0] * qidx_arr.shape[1];
    for (size_t i = 0; i < qN; ++i) {
        long long v = q_read_as_ll(i);
        qmin = std::min(qmin, v);
        qmax = std::max(qmax, v);
    }

    // Detect 1-based indexing if values look like [1..Nin]
    const bool one_based = (qmin >= 1 && qmax <= (long long)Nin && qmax > 0);

    std::cout << "surface-color debug: Nin=" << Nin
              << " Nma=" << Nma
              << " qmin=" << qmin << " qmax=" << qmax
              << " qcols=" << qcols
              << " one_based=" << (one_based ? "true" : "false")
              << "\n";

    auto to_zero_based = [&](long long qi) -> long long {
        return one_based ? (qi - 1) : qi;
    };

    // For each surface point j, store an aggregated radius proxy (min is robust)
    std::vector<double> r_per_in(Nin, std::numeric_limits<double>::infinity());
    std::vector<int>    hits(Nin, 0);

    auto dist_ma_to_in = [&](size_t i_ma, size_t j_in) -> double {
        const double mx = ma[3*i_ma + 0], my = ma[3*i_ma + 1], mz = ma[3*i_ma + 2];
        const double sx = in[3*j_in + 0], sy = in[3*j_in + 1], sz = in[3*j_in + 2];
        const double dx = mx - sx, dy = my - sy, dz = mz - sz;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    };

    // Push MA radius proxies back onto the referenced surface indices
    for (size_t i = 0; i < Nma; ++i) {
        // compute proxy radius for this MA point:
        // distance to its referenced surface point(s)
        double proxy_r = std::numeric_limits<double>::infinity();

        for (size_t c = 0; c < qcols; ++c) {
            long long qi_raw = q_at(i, c);
            long long qi = to_zero_based(qi_raw);

            if (qi < 0) continue;
            if ((size_t)qi >= Nin) continue;

            proxy_r = std::min(proxy_r, dist_ma_to_in(i, (size_t)qi));
        }

        if (!std::isfinite(proxy_r)) continue; // this MA point had no valid refs

        // assign this proxy radius to each referenced input point
        for (size_t c = 0; c < qcols; ++c) {
            long long qi_raw = q_at(i, c);
            long long qi = to_zero_based(qi_raw);

            if (qi < 0) continue;
            if ((size_t)qi >= Nin) continue;

            // aggregation choice:
            // MIN is stable (if multiple balls reference the same surface point, keep smallest).
            r_per_in[(size_t)qi] = std::min(r_per_in[(size_t)qi], proxy_r);
            hits[(size_t)qi] += 1;
        }
    }

    // Replace untouched points (never referenced) with something safe
    // so they still get a color (here: use 0.0 => will map to red)
    size_t unhit = 0;
    for (size_t j = 0; j < Nin; ++j) {
        if (!std::isfinite(r_per_in[j])) {
            r_per_in[j] = 0.0;
            unhit++;
        }
    }
    std::cout << "surface-color: unhit surface points = " << unhit << " / " << Nin << "\n";

    // Normalize radii to 0..1 for coloring
    double rmin = std::numeric_limits<double>::infinity();
    double rmax = 0.0;
    for (size_t j = 0; j < Nin; ++j) {
        rmin = std::min(rmin, r_per_in[j]);
        rmax = std::max(rmax, r_per_in[j]);
    }
    const double denom = (rmax - rmin) + 1e-12;

    // Write COFF with SURFACE coords
    std::filesystem::create_directories(std::filesystem::path(coff_file).parent_path());
    std::ofstream out(coff_file);
    if (!out) throw std::runtime_error("Cannot open " + coff_file);

    out << std::setprecision(17) << std::fixed;
    out << "COFF\n";
    out << Nin << " 0 0\n";

    for (size_t j = 0; j < Nin; ++j) {
        const double t = (r_per_in[j] - rmin) / denom; // 0..1

        // red (small) -> blue (large)
        const int R = (int)std::lround(255.0 * (1.0 - t));
        const int G = 0;
        const int B = (int)std::lround(255.0 * t);
        const int A = 255;

        out << in[3*j + 0] << " "
            << in[3*j + 1] << " "
            << in[3*j + 2] << " "
            << R << " " << G << " " << B << " " << A << "\n";
    }

    std::cout << "Wrote colored SURFACE points to " << coff_file << "\n";
}

void write_normals_as_obj_lines(const std::string& obj_file,
                                const std::vector<Point_3>& pts,
                                const std::vector<Vector_3>& normals,
                                double scale)   // e.g. 1.0, 5.0, 10.0 depending on your units
{
    if (pts.size() != normals.size())
        throw std::runtime_error("write_normals_as_obj_lines: pts/normals size mismatch");

    std::ofstream out(obj_file);
    if (!out) throw std::runtime_error("Cannot open " + obj_file);

    out << std::setprecision(17) << std::fixed;

    // We write 2 vertices per normal: start and end
    // Then a line element connecting them
    for (size_t i = 0; i < pts.size(); ++i) {
        const auto& p = pts[i];
        Vector_3 n = normals[i];

        // normalize (safety)
        const double len2 = CGAL::to_double(n.squared_length());
        if (len2 > 0.0) n = n / std::sqrt(len2);

        Point_3 q = p + scale * n;

        out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
        out << "v " << q.x() << " " << q.y() << " " << q.z() << "\n";
    }

    // OBJ indices are 1-based
    for (size_t i = 0; i < pts.size(); ++i) {
        const size_t v0 = 2*i + 1;
        const size_t v1 = 2*i + 2;
        out << "l " << v0 << " " << v1 << "\n";
    }

    std::cout << "Wrote normals as OBJ lines: " << obj_file
              << " (lines=" << pts.size() << ", scale=" << scale << ")\n";
}

double bbox_diag(const Mesh& mesh) {
    CGAL::Bbox_3 b = PMP::bbox(mesh);
    const double dx = b.xmax() - b.xmin();
    const double dy = b.ymax() - b.ymin();
    const double dz = b.zmax() - b.zmin();
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

void write_surface_colored_by_power_distance(const std::string& ma_coords_out_npy,
                                             const std::string& ma_qidx_out_npy,
                                             const std::string& surface_coords_npy,
                                             const std::string& coff_file)
{
    cnpy::NpyArray ma_arr   = cnpy::npy_load(ma_coords_out_npy);
    cnpy::NpyArray qidx_arr = cnpy::npy_load(ma_qidx_out_npy);
    cnpy::NpyArray in_arr   = cnpy::npy_load(surface_coords_npy);

    auto check_coords = [&](const cnpy::NpyArray& a, const std::string& name) {
        if (a.shape.size() != 2 || a.shape[1] != 3) throw std::runtime_error(name + " must be Nx3");
        if (a.word_size != sizeof(float)) throw std::runtime_error(name + " must be float32");
    };
    check_coords(ma_arr, "ma_coords_out");
    check_coords(in_arr, "surface_coords");

    const size_t Nma = ma_arr.shape[0];
    const size_t Nin = in_arr.shape[0];

    if (!((qidx_arr.shape.size()==1 && qidx_arr.shape[0]==Nma) ||
          (qidx_arr.shape.size()==2 && qidx_arr.shape[0]==Nma && qidx_arr.shape[1]==1)))
        throw std::runtime_error("ma_qidx_out must be length Nma (1D) or Nma x 1");

    const float* ma = reinterpret_cast<const float*>(ma_arr.data);
    const float* in = reinterpret_cast<const float*>(in_arr.data);

    auto q_at = [&](size_t i) -> long long {
        const char* raw = reinterpret_cast<const char*>(qidx_arr.data);
        const size_t flat = (qidx_arr.shape.size()==1) ? i : i;
        if (qidx_arr.word_size == 4) return (long long)reinterpret_cast<const int32_t*>(raw)[flat];
        if (qidx_arr.word_size == 8) return (long long)reinterpret_cast<const int64_t*>(raw)[flat];
        throw std::runtime_error("Unsupported qidx dtype (expected 4 or 8 bytes)");
    };

    // Build valid balls: centers + radius proxy
    struct Ball { float cx, cy, cz; float r; };
    std::vector<Ball> balls;
    balls.reserve(Nma);

    size_t neg = 0, bad = 0;

    for (size_t i = 0; i < Nma; ++i) {
        long long qi = q_at(i);
        if (qi < 0) { neg++; continue; }
        if ((size_t)qi >= Nin) { bad++; continue; }

        const float cx = ma[3*i + 0], cy = ma[3*i + 1], cz = ma[3*i + 2];
        const size_t j = (size_t)qi;
        const float sx = in[3*j + 0], sy = in[3*j + 1], sz = in[3*j + 2];

        const float dx = cx - sx, dy = cy - sy, dz = cz - sz;
        const float r  = std::sqrt(dx*dx + dy*dy + dz*dz);

        balls.push_back({cx,cy,cz,r});
    }

    std::cout << "power-distance: Nin=" << Nin
              << " Nma=" << Nma
              << " valid_balls=" << balls.size()
              << " qidx_neg=" << neg
              << " qidx_bad=" << bad << "\n";

    if (balls.empty())
        throw std::runtime_error("No valid balls (all qidx are -1). Cannot assign ownership.");

    // Assign each surface point to a ball by min power distance: ||p-c||^2 - r^2
    std::vector<float> owned_r(Nin, 0.f);
    float rmin = std::numeric_limits<float>::infinity();
    float rmax = 0.f;

    for (size_t j = 0; j < Nin; ++j) {
        const float px = in[3*j + 0], py = in[3*j + 1], pz = in[3*j + 2];

        float best_pd = std::numeric_limits<float>::infinity();
        size_t best_i = 0;

        for (size_t i = 0; i < balls.size(); ++i) {
            const float dx = px - balls[i].cx;
            const float dy = py - balls[i].cy;
            const float dz = pz - balls[i].cz;

            const float d2 = dx*dx + dy*dy + dz*dz;
            const float pd = d2 - balls[i].r * balls[i].r;

            if (pd < best_pd) { best_pd = pd; best_i = i; }
        }

        const float r = balls[best_i].r;
        owned_r[j] = r;
        rmin = std::min(rmin, r);
        rmax = std::max(rmax, r);
    }

    // Write COFF at surface coords, colored by owned radius (red small -> blue large)
    std::filesystem::create_directories(std::filesystem::path(coff_file).parent_path());
    std::ofstream out(coff_file);
    if (!out) throw std::runtime_error("Cannot open " + coff_file);

    out << std::setprecision(17) << std::fixed;
    out << "COFF\n";
    out << Nin << " 0 0\n";

    const float denom = (rmax - rmin) + 1e-12f;

    for (size_t j = 0; j < Nin; ++j) {
        const float t = (owned_r[j] - rmin) / denom;

        const int R = (int)std::lround(255.0 * (1.0 - t));
        const int G = 0;
        const int B = (int)std::lround(255.0 * t);
        const int A = 255;

        out << in[3*j + 0] << " "
            << in[3*j + 1] << " "
            << in[3*j + 2] << " "
            << R << " " << G << " " << B << " " << A << "\n";
    }

    std::cout << "Wrote power-distance owned surface COFF to " << coff_file
              << " rmin=" << rmin << " rmax=" << rmax << "\n";
}


// ----------------------------------GEMINI-----------------------------------
// Helper to calculate squared Euclidean distance
float dist_sq(float x1, float y1, float z1, float x2, float y2, float z2) {
    float dx = x1 - x2;
    float dy = y1 - y2;
    float dz = z1 - z2;
    return dx*dx + dy*dy + dz*dz;
}

static void require_2d_n3(const cnpy::NpyArray& a, const std::string& name)
{
    if (a.shape.size() != 2 || a.shape[1] != 3)
        throw std::runtime_error(name + " must have shape (N,3)");
}

static void require_1d_n(const cnpy::NpyArray& a, size_t n, const std::string& name)
{
    if (a.shape.size() != 1 || a.shape[0] != n)
        throw std::runtime_error(name + " must have shape (N,) matching coords");
}

#include <cnpy.h>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>

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




