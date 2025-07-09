#ifndef AMRACUT_TESTBENCH_MESH_UTILS_H
#define AMRACUT_TESTBENCH_MESH_UTILS_H

#include <vector>

namespace amracut_testbench
{
  ElementType GetElementType(const std::string &mesh_file_path, MPI_Comm comm);

  template <class T>
  void GetInitialElementsDistribution(const std::string &mesh_file_path, std::vector<T> &elements_out,
                                      ElementType element_type, MPI_Comm comm);
} // namespace amracut_testbench

#include "mesh_utils.tcc"

#endif