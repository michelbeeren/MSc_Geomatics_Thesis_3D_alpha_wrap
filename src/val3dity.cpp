//
// Created by Michel Beeren on 26/01/2026.
//
#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <iomanip>
#include <filesystem>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;

// =================== INPUT FILE TO CHECK VALIDITY ====================
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

// ================== MAKE AND CHECK VAL3DITY REPORT ====================
bool report_is_valid(const std::string& report_path)
{
    std::ifstream file(report_path);
    if (!file.is_open())
        return false;

    std::string content(
        (std::istreambuf_iterator<char>(file)),
         std::istreambuf_iterator<char>()
    );

    // val3dity writes: "validity": true / false
    return content.find("\"validity\": true") != std::string::npos;
}

bool run_val3dity_and_check(const std::string& input_path, const std::string& report_path, const std::string& primitive = "")
{
    std::string exe = "/opt/homebrew/bin/val3dity";

    std::string cmd = exe + " --snap_tol 1e-06 --report " + report_path + " ";

    if (!primitive.empty())
        cmd += "-p " + primitive + " ";

    cmd += input_path + " > /dev/null 2>&1";

    int system_code = std::system(cmd.c_str());
    if (system_code != 0) return false;

    return report_is_valid(report_path); // parses JSON "validity"
}

bool valid_mesh_boolean(const Mesh& mesh_)
{
    const std::string base =
        "../data/Test_3D_Alphawrap/Output/val3dity_check/";

    const std::string cityjson_path = base + "check_me.json";
    const std::string obj_path      = base + "check_me.obj";
    const std::string report_path   = base + "report.json";

    bool mesh_is_valid = false;

    try {
        // OPTION A — CityJSON (what you already have)
        write_cityjson_from_mesh(mesh_, cityjson_path);
        mesh_is_valid = run_val3dity_and_check(cityjson_path, report_path);

        // OPTION B — OBJ instead (uncomment if you want this path)
        /*
        write_obj_from_mesh(mesh_, obj_path);
        mesh_is_valid = run_val3dity_and_check(obj_path, report_path, "Solid");
        */
    }
    catch (const std::exception&) {
        mesh_is_valid = false;
    }

    //
    if (mesh_is_valid) {
        std::cout << "❤️❤️❤️Mesh is ✅VALID ✅ :) 😍😎👍🏼" << std::endl;
        return true;
    } else {
        std::cout << " ⚠️⚠️⚠️ Mesh is ❌INVALID ❌ :( 👎🏼👎🏼👎🏼" << std::endl;
        return false;
    }
}

bool valid_file_boolean(const std::string& input_path)
{
    const std::string report_path =
        "../data/Test_3D_Alphawrap/Output/val3dity_check/report.json";

    // detect extension
    std::string ext;
    auto pos = input_path.find_last_of('.');
    if (pos != std::string::npos)
        ext = input_path.substr(pos + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });

    bool is_valid = false;

    try {
        if (ext == "obj" || ext == "off") {
            // OBJ/OFF need an explicit primitive
            is_valid = run_val3dity_and_check(input_path, report_path, "Solid");
        } else {
            // CityJSON (and other formats val3dity understands with embedded type)
            is_valid = run_val3dity_and_check(input_path, report_path);
        }
    } catch (...) {
        is_valid = false;
    }

    if (is_valid) {
        std::cout << "Mesh is VALID :)" << std::endl;
        return true;
    } else {
        std::cout << "Mesh is INVALID :(" << std::endl;
        return false;
    }
}



