#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cinttypes>
#include <stdexcept>
#include <numeric>
#include <algorithm>

#include "nlohmann-json/json.hpp"

#include "utils/utils.hpp"

namespace amracut_testbench
{
  template <>
  std::string VectorToString(std::vector<uint8_t> vec)
  {
    std::ostringstream output;
    output << "[ ";
    for (auto element : vec)
    {
      output << +element << ", ";
    }
    output << "]\n";
    return output.str();
  }

  void ExportMetricsToJson(
      std::string mesh_file,
      int partition_count, uint64_t global_vertex_count,
      uint64_t sfc_partition_time,
      std::vector<uint64_t> &sfc_partition_sizes, std::vector<uint64_t> &sfc_partition_boundaries, std::vector<uint64_t> &sfc_partition_cuts,
      uint64_t amracut_labeling_time,
      std::vector<uint64_t> &amracut_partition_sizes, std::vector<uint64_t> &amracut_partition_boundaries, std::vector<uint64_t> &amracut_partition_cuts,
      std::string metrics_out_file_path)
  {

    uint64_t total_weight_sum = std::accumulate(sfc_partition_sizes.begin(), sfc_partition_sizes.end(), (uint64_t)0);
    uint64_t ideal_partition_weight = total_weight_sum / partition_count;

    nlohmann::json output_json =
        {
            {"mesh_file", mesh_file},
            {"np", partition_count},
            {"n", global_vertex_count},
            {"total_weight_sum", total_weight_sum},

            {"SFC_morton_boundary_ratio", std::accumulate(sfc_partition_boundaries.begin(), sfc_partition_boundaries.end(), static_cast<uint64_t>(0)) / static_cast<float>(total_weight_sum)},
            {"SFC_morton_cut_ratio", std::accumulate(sfc_partition_cuts.begin(), sfc_partition_cuts.end(), static_cast<uint64_t>(0)) / static_cast<float>(total_weight_sum * 2)}, // dividing by 2 because egdes are counted twice
            {"SFC_morton_rho_max", static_cast<float>(*std::max_element(sfc_partition_sizes.begin(), sfc_partition_sizes.end())) / ideal_partition_weight},
            {"SFC_morton_rho_min", static_cast<float>(*std::min_element(sfc_partition_sizes.begin(), sfc_partition_sizes.end())) / ideal_partition_weight},
            {"SFC_morton_partition_sizes", sfc_partition_sizes},
            {"SFC_morton_partition_boundaries", sfc_partition_boundaries},
            {"SFC_morton_partition_cuts", sfc_partition_cuts},
            {"SFC_morton_partition_time", sfc_partition_time},

            {"AMRaCut_boundary_ratio", std::accumulate(amracut_partition_boundaries.begin(), amracut_partition_boundaries.end(), static_cast<uint32_t>(0)) / static_cast<float>(total_weight_sum)},
            {"AMRaCut_cut_ratio", std::accumulate(amracut_partition_cuts.begin(), amracut_partition_cuts.end(), static_cast<uint32_t>(0)) / static_cast<float>(total_weight_sum * 2)},
            {"AMRaCut_rho_max", static_cast<float>(*std::max_element(amracut_partition_sizes.begin(), amracut_partition_sizes.end())) / ideal_partition_weight},
            {"AMRaCut_rho_min", static_cast<float>(*std::min_element(amracut_partition_sizes.begin(), amracut_partition_sizes.end())) / ideal_partition_weight},
            {"AMRaCut_partition_sizes", amracut_partition_sizes},
            {"AMRaCut_partition_boundaries", amracut_partition_boundaries},
            {"AMRaCut_partition_cuts", amracut_partition_cuts},
            {"AMRaCut_labeling_time", amracut_labeling_time}

        };
    print_log("SFC_morton_rho_max: ", output_json["SFC_morton_rho_max"]);
    print_log("amracut_rho_max: ", output_json["AMRaCut_rho_max"]);
    print_log("");

    print_log("SFC_morton_cut_ratio: ", output_json["SFC_morton_cut_ratio"]);
    print_log("amracut_cut_ratio: ", output_json["AMRaCut_cut_ratio"]);
    print_log("");

    print_log("SFC_morton_boundary_ratio: ", output_json["SFC_morton_boundary_ratio"]);
    print_log("amracut_boundary_ratio: ", output_json["AMRaCut_boundary_ratio"]);
    print_log("");

    print_log("SFC_morton_bdry_max: ", (*std::max_element(sfc_partition_boundaries.begin(), sfc_partition_boundaries.end())));
    print_log("amracut_bdry_max: ", (*std::max_element(amracut_partition_boundaries.begin(), amracut_partition_boundaries.end())));
    print_log("");

    print_log("SFC_morton_cut_max: ", (*std::max_element(sfc_partition_cuts.begin(), sfc_partition_cuts.end())));
    print_log("amracut_cut_max: ", (*std::max_element(amracut_partition_cuts.begin(), amracut_partition_cuts.end())));
    print_log("");

    std::ofstream file;
    file.open(metrics_out_file_path, std::ios::app);

    if (!file.is_open())
    {
      std::cerr << "Error: Could not open the file '" << metrics_out_file_path << "' - " << std::strerror(errno) << std::endl;
      throw std::runtime_error("Error in opening file: " + metrics_out_file_path);
    }

    file << output_json.dump(-1) << std::endl;

    file.close();
  }

} // namespace amracut_textbench
