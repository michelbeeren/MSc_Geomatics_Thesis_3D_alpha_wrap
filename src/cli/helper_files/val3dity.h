#ifndef THESIS_CLI_VAL3DITY_H
#define THESIS_CLI_VAL3DITY_H

#include "alpha_wrapping.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace cli_helpers
{

inline void write_cityjson_from_mesh(const Mesh& mesh, const std::string& filename)
{
  std::ofstream out(filename);
  if(!out)
    throw std::runtime_error("Cannot open output CityJSON file: " + filename);

  out << std::fixed << std::setprecision(17);

  using VI = Mesh::Vertex_index;
  using FI = Mesh::Face_index;

  std::vector<VI> index_to_vertex;
  index_to_vertex.reserve(num_vertices(mesh));
  std::unordered_map<VI, int> vidx;

  int idx = 0;
  for(VI v : mesh.vertices())
  {
    index_to_vertex.push_back(v);
    vidx[v] = idx++;
  }

  std::vector<Point_3> verts;
  verts.reserve(index_to_vertex.size());
  for(VI v : index_to_vertex)
    verts.push_back(mesh.point(v));

  std::vector<std::vector<int>> faces;
  faces.reserve(num_faces(mesh));
  for(FI f : mesh.faces())
  {
    std::vector<int> face;
    for(VI v : CGAL::vertices_around_face(mesh.halfedge(f), mesh))
      face.push_back(vidx[v]);
    if(face.size() >= 3)
      faces.push_back(std::move(face));
  }

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
  out << "            [\n";

  for(std::size_t i = 0; i < faces.size(); ++i)
  {
    out << "              [ [";
    for(std::size_t j = 0; j < faces[i].size(); ++j)
    {
      out << faces[i][j];
      if(j + 1 < faces[i].size())
        out << ", ";
    }
    out << "] ]";
    if(i + 1 < faces.size())
      out << ",";
    out << "\n";
  }

  out << "            ]\n";
  out << "          ]\n";
  out << "        }\n";
  out << "      ]\n";
  out << "    }\n";
  out << "  },\n";
  out << "  \"vertices\": [\n";

  for(std::size_t k = 0; k < verts.size(); ++k)
  {
    out << "    [" << verts[k].x() << ", " << verts[k].y() << ", " << verts[k].z() << "]";
    if(k + 1 < verts.size())
      out << ",";
    out << "\n";
  }

  out << "  ]\n";
  out << "}\n";
}

inline bool report_is_valid(const std::string& report_path)
{
  std::ifstream file(report_path);
  if(!file.is_open())
    return false;

  const std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
  return content.find("\"validity\": true") != std::string::npos;
}

inline std::string val3dity_executable()
{
#ifdef VAL3DITY_PATH
  return std::string(VAL3DITY_PATH);
#else
  return std::string("/opt/homebrew/bin/val3dity");
#endif
}

inline bool run_val3dity_and_check(const std::string& input_path,
                                   const std::string& report_path,
                                   const std::string& primitive = "")
{
  const std::string exe = val3dity_executable();
  std::string cmd = "\"" + exe + "\" --snap_tol 1e-06 --report \"" + report_path + "\" ";

  if(!primitive.empty())
    cmd += "-p " + primitive + " ";

  cmd += "\"" + input_path + "\" > /dev/null 2>&1";
  const int system_code = std::system(cmd.c_str());
  if(system_code != 0)
    return false;

  return report_is_valid(report_path);
}

inline bool valid_mesh_boolean(const Mesh& mesh)
{
  const std::filesystem::path base = std::filesystem::temp_directory_path() / "alpha_wrap_cli_val3dity";
  std::filesystem::create_directories(base);

  const std::string cityjson_path = (base / "check_me.json").string();
  const std::string report_path = (base / "report.json").string();

  bool mesh_is_valid = false;
  try
  {
    write_cityjson_from_mesh(mesh, cityjson_path);
    mesh_is_valid = run_val3dity_and_check(cityjson_path, report_path);
  }
  catch(const std::exception&)
  {
    mesh_is_valid = false;
  }

  std::cout << (mesh_is_valid ? "Validation result: VALID\n" : "Validation result: INVALID\n");
  return mesh_is_valid;
}

inline bool valid_file_boolean(const std::string& input_path)
{
  const std::filesystem::path base = std::filesystem::temp_directory_path() / "alpha_wrap_cli_val3dity";
  std::filesystem::create_directories(base);
  const std::string report_path = (base / "report.json").string();

  std::string ext;
  const auto pos = input_path.find_last_of('.');
  if(pos != std::string::npos)
    ext = input_path.substr(pos + 1);

  std::transform(ext.begin(), ext.end(), ext.begin(),
                 [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });

  bool is_valid = false;
  try
  {
    if(ext == "obj" || ext == "off")
      is_valid = run_val3dity_and_check(input_path, report_path, "Solid");
    else
      is_valid = run_val3dity_and_check(input_path, report_path);
  }
  catch(...)
  {
    is_valid = false;
  }

  std::cout << (is_valid ? "Validation result: VALID\n" : "Validation result: INVALID\n");
  return is_valid;
}

} // namespace cli_helpers

#endif // THESIS_CLI_VAL3DITY_H
