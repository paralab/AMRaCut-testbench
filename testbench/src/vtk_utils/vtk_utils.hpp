#ifndef AMRACUT_TESTBENCH_VTK_UTIS_H
#define AMRACUT_TESTBENCH_VTK_UTIS_H

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
      const std::vector<uint64_t> &amracut_partition_labels,
      const std::string output_filename);
}


#endif