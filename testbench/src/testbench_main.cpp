#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include "mpi.h"

#include "dtypes/dtypes.hpp"
#include "mesh_utils/mesh_utils.hpp"
#include "utils/utils.hpp"
#include "usort/parUtils.h"
#include "usort/ompUtils.h"
#include "sfc/sfc.hpp"

#ifdef ENABLE_VTK_FEATURES
#include "vtk_utils/vtk_utils.hpp"
#endif




/**
 * steps
 *
 * 1.
 *
 */

template <class T>
amracut_testbench::PartitionStatus ReadAndDistributeSFC(std::string mesh_file_path, 
                               amracut_testbench::ElementType element_type, 
                               std::vector<T> &elements_out, MPI_Comm comm);


int main(int argc, char *argv[])
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int procs_n, my_rank, hostname_len;
  char hostname[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &procs_n);
  MPI_Comm_rank(comm, &my_rank);
  MPI_Get_processor_name(hostname, &hostname_len);



  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <mesh file path>" << std::endl;
    return 1;
  }
  const std::string mesh_file_path = argv[1];

  if (!my_rank)
  {
    amracut_testbench::print_log("running on", procs_n, "MPI processs");
    amracut_testbench::print_log("partitioning: ", mesh_file_path);
  }

  amracut_testbench::ElementType element_type = amracut_testbench::GetElementType(mesh_file_path, comm);

  std::vector<amracut_testbench::HexElement> local_sfc_elements_hex;
  std::vector<amracut_testbench::TetElement> local_sfc_elements_tet;
  std::vector<amracut_testbench::Element> local_sfc_elements;

  std::vector<amracut_testbench::Element> global_elements;


  std::vector<std::pair<uint64_t, uint64_t>> local_connected_element_pairs;
  std::vector<std::pair<uint64_t, uint64_t>> boundary_connected_element_pairs;

  std::vector<amracut_testbench::ElementFace> local_unconnected_elements_faces;


  int local_element_count;
  int global_element_count;

  std::vector<int> proc_element_counts(procs_n);
  std::vector<int> proc_element_counts_scanned(procs_n);

  amracut_testbench::PartitionStatus SFC_status;

  switch (element_type)
  {
  case amracut_testbench::ElementType::TET:
  {
    SFC_status = ReadAndDistributeSFC(mesh_file_path, element_type, local_sfc_elements_tet, comm);
    // amracut_testbench::print_log_mpi(my_rank, "elements:", local_sfc_elements_tet.size());
    local_element_count = static_cast<int>(local_sfc_elements_tet.size());
    
    break;                                                    
  }
  case amracut_testbench::ElementType::HEX:
  {
    SFC_status = ReadAndDistributeSFC(mesh_file_path, element_type, local_sfc_elements_hex, comm);
    // amracut_testbench::print_log_mpi(my_rank, "elements:", local_sfc_elements_hex.size()); 
    local_element_count = static_cast<int>(local_sfc_elements_hex.size());
    break;                                                     

  }
  default:
  {
    throw std::runtime_error("Unknown element type");
    break;
  }
  }

  MPI_Allgather(&local_element_count, 1, MPI_INT, proc_element_counts.data(),1,MPI_INT, comm);

  MPI_Allreduce(&local_element_count, &global_element_count, 1, MPI_INT, MPI_SUM, comm);

  omp_par::scan(&proc_element_counts[0], &proc_element_counts_scanned[0], procs_n);

  local_sfc_elements.resize(local_element_count);
  uint64_t global_idx_start = static_cast<uint64_t>(proc_element_counts_scanned[my_rank]);

  switch (element_type)
  {
  case amracut_testbench::ElementType::TET:
  {
    for (int local_elem_i = 0; local_elem_i < local_element_count; local_elem_i++)
    {
      local_sfc_elements_tet[local_elem_i].global_idx = global_idx_start + local_elem_i;
      
      local_sfc_elements[local_elem_i].global_idx = global_idx_start + local_elem_i;
      local_sfc_elements[local_elem_i].x = local_sfc_elements_tet[local_elem_i].x;
      local_sfc_elements[local_elem_i].y = local_sfc_elements_tet[local_elem_i].y;
      local_sfc_elements[local_elem_i].z = local_sfc_elements_tet[local_elem_i].z;
    }
    amracut_testbench::ResolveLocalElementConnectivity(local_sfc_elements_tet, amracut_testbench::ElementType::TET, 
                                                       local_connected_element_pairs, local_unconnected_elements_faces);
    break;
  }
  case amracut_testbench::ElementType::HEX:
  {

    for (int local_elem_i = 0; local_elem_i < local_element_count; local_elem_i++)
    {
      local_sfc_elements_hex[local_elem_i].global_idx = global_idx_start + local_elem_i;

      local_sfc_elements[local_elem_i].global_idx = global_idx_start + local_elem_i;
      local_sfc_elements[local_elem_i].x = local_sfc_elements_hex[local_elem_i].x;
      local_sfc_elements[local_elem_i].y = local_sfc_elements_hex[local_elem_i].y;
      local_sfc_elements[local_elem_i].z = local_sfc_elements_hex[local_elem_i].z;
    }
    amracut_testbench::ResolveLocalElementConnectivity(local_sfc_elements_hex, amracut_testbench::ElementType::HEX, 
                                                       local_connected_element_pairs, local_unconnected_elements_faces);
    break;
  }

  default:
  {
    throw std::runtime_error("Unknown element type");
    break;
  }
  }



  amracut_testbench::ResolveBoundaryElementConnectivity(local_unconnected_elements_faces,
                                                        proc_element_counts, proc_element_counts_scanned, 
                                                        boundary_connected_element_pairs, comm);







  std::vector<uint64_t> local_sfc_partition_labels(local_element_count, my_rank);

  #ifdef ENABLE_VTK_FEATURES
  std::vector<amracut_testbench::Element> global_sfc_elements;
  std::vector<uint64_t> global_sfc_partition_labels;

  if (!my_rank)
  {
    global_sfc_partition_labels.resize(global_element_count);
    global_sfc_elements.resize(global_element_count);
  }

  MPI_Gatherv(local_sfc_elements.data(), local_element_count, par::Mpi_datatype<amracut_testbench::Element>::value(), global_sfc_elements.data(),
              proc_element_counts.data(),proc_element_counts_scanned.data(),par::Mpi_datatype<amracut_testbench::Element>::value(), 0, comm);
  MPI_Barrier(comm);
  
  MPI_Gatherv(local_sfc_partition_labels.data(), local_element_count, MPI_UINT64_T, global_sfc_partition_labels.data(),
              proc_element_counts.data(),proc_element_counts_scanned.data(),MPI_UINT64_T, 0, comm);

  if (!my_rank)
  {
    amracut_testbench::export_parititions_to_vtk(global_sfc_elements, global_sfc_partition_labels, "partitions");
  }
  

  #endif


  MPI_Finalize();



  return 0;
}

template <class T>
amracut_testbench::PartitionStatus ReadAndDistributeSFC(std::string mesh_file_path, 
                               amracut_testbench::ElementType element_type, 
                               std::vector<T> &elements_out, MPI_Comm comm)
{

  int procs_n, my_rank;
  MPI_Comm_size(comm, &procs_n);
  MPI_Comm_rank(comm, &my_rank);

  std::vector<T> initial_elements; // before SFC
  amracut_testbench::GetInitialElementsDistribution(mesh_file_path, 
                                                    initial_elements, 
                                                    element_type, comm);
  amracut_testbench::SetMortonEncoding(initial_elements, element_type, comm);
  elements_out.resize(initial_elements.size());
  if (!my_rank)
  {
    amracut_testbench::print_log("starting sfc sort");
  }

  MPI_Barrier(comm);
  std::vector<T> initial_elements_cpy(initial_elements);
  par::sampleSort<T>(initial_elements_cpy,elements_out,comm);       //warmup?

  MPI_Barrier(comm);
  auto start = std::chrono::high_resolution_clock::now();
  par::sampleSort<T>(initial_elements,elements_out,comm);
  MPI_Barrier(comm);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  if (!my_rank)
  {
    amracut_testbench::print_log("global sfc sort done");
    amracut_testbench::print_log("SFC sort time:\t", duration.count(), " us");
  }

  amracut_testbench::PartitionStatus status;
  status.return_code = 0;
  status.time_us = duration.count();

  return status;
}