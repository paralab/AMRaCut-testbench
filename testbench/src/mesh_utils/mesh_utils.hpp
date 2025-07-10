#ifndef AMRACUT_TESTBENCH_MESH_UTILS_H
#define AMRACUT_TESTBENCH_MESH_UTILS_H

#include <vector>

namespace amracut_testbench
{
  ElementType GetElementType(const std::string &mesh_file_path, MPI_Comm comm);

  template <class T>
  void GetInitialElementsDistribution(const std::string &mesh_file_path, std::vector<T> &elements_out,
                                      ElementType element_type, MPI_Comm comm);

  template <class T>
  void ResolveLocalElementConnectivity(const std::vector<T> &elements, ElementType element_type,
                                       std::vector<std::pair<uint64_t, uint64_t>> &connected_element_idx_pairs_out,
                                       std::vector<ElementFace> &unconnected_elements_faces_out);

  void ResolveBoundaryElementConnectivity(std::vector<ElementFace> &unpaired_element_faces,
                                          std::vector<int> &proc_element_counts,
                                          std::vector<int> &proc_element_counts_scanned,
                                          std::vector<std::pair<uint64_t, uint64_t>> &boundary_connected_element_idx_pairs_out,
                                          MPI_Comm comm);                                                                             
} // namespace amracut_testbench

#include "mesh_utils/mesh_utils.tcc"

#endif