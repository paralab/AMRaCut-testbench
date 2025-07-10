#include <vector>
#include <string>
#include <cinttypes>
#include <stdexcept>

#include "dtypes/dtypes.hpp"

namespace amracut_testbench
{

  void export_parititions_to_vtk(
      const std::vector<amracut_testbench::Element> &elements,
      const std::vector<uint64_t> &sfc_partition_labels,
      // const std::vector<uint64_t> &p2,
      // const std::vector<uint64_t> &p3,
      const std::string output_filename);
}