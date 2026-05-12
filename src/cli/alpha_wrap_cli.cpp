#include "helper_files/alpha_wrapping.h"
#include "helper_files/hausdorff.h"
#include "helper_files/val3dity.h"

#include <chrono>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

struct Arguments
{
  std::string input_path;
  std::string output_path;
  double alpha = 0.0;
  double offset = 0.0;
  double tau = 0.0;
  bool validate = false;
};

void print_usage(const char* executable)
{
  std::cerr << "Usage:\n"
            << "  " << executable << " <input_path> <output_path> <alpha> <offset> <tau> [--validate]\n\n"
            << "Notes:\n"
            << "  - alpha and offset are relative to the input bbox diagonal\n"
            << "  - tau is max_distance_to_input_in_offsets and must be > 1.0\n"
            << "  - add --validate to run mesh validation\n";
}

bool parse_arguments(const int argc, char** argv, Arguments& args)
{
  if(argc != 6 && argc != 7)
    return false;

  args.input_path = argv[1];
  args.output_path = argv[2];

  try
  {
    args.alpha = std::stod(argv[3]);
    args.offset = std::stod(argv[4]);
    args.tau = std::stod(argv[5]);
  }
  catch(const std::exception&)
  {
    return false;
  }

  if(argc == 7)
  {
    if(std::string(argv[6]) != "--validate")
      return false;
    args.validate = true;
  }

  return true;
}

int main(int argc, char** argv)
{
  Arguments args;
  if(!parse_arguments(argc, argv, args))
  {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  try
  {
    const cli_helpers::Wrap_request request{
        args.input_path,
        args.alpha,
        args.offset,
        args.tau};

    const auto wrap_start = std::chrono::steady_clock::now();
    cli_helpers::Wrap_result result = cli_helpers::run_alpha_wrap(request);
    const auto wrap_end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> wrap_elapsed = wrap_end - wrap_start;

    std::cout << "Wrap result: " << num_vertices(result.wrap_mesh) << " vertices, "
              << num_faces(result.wrap_mesh) << " faces\n";
    std::cout << "Total wrapping time: " << std::fixed << std::setprecision(3)
              << wrap_elapsed.count() << " s\n";
    std::cout << "Writing output mesh: " << args.output_path << "\n";
    if(!cli_helpers::write_output_mesh(args.output_path, result.wrap_mesh))
    {
      std::cerr << "Error: could not write output mesh.\n";
      return EXIT_FAILURE;
    }

    if(args.validate)
    {
      const bool is_valid = cli_helpers::valid_mesh_boolean(result.wrap_mesh);
      if(!is_valid)
        std::cerr << "Warning: output mesh did not pass validation.\n";
    }

    std::cout << "Done.\n";
    return EXIT_SUCCESS;
  }
  catch(const std::exception& e)
  {
    std::cerr << "Error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
}
