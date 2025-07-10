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
#include "dgraph/dgraph.hpp"


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



  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " <mesh file path> <metrics output json path>" << std::endl;
    return 1;
  }
  const std::string mesh_file_path = argv[1];
  const std::string output_file_path = argv[2];


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

  if (!my_rank)
  {
    amracut_testbench::print_log("starting elemenet connectivity");
  }

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

  if (!my_rank)
  {
    amracut_testbench::print_log("elemenet connectivity done");
    amracut_testbench::print_log("starting graph creation");
  }

  amracut_testbench::DGraph dist_graph(local_connected_element_pairs, boundary_connected_element_pairs,
                                       proc_element_counts, proc_element_counts_scanned ,DIST_GRAPH_UNWEIGHTED , comm);
  if (!my_rank)
  {
    amracut_testbench::print_log("graph created");
  }

  std::vector<uint64_t> local_amracut_partition_labels(local_element_count);

  if (!my_rank)
  {
    amracut_testbench::print_log("starting AMRaCut partitioning");
  }
  amracut_testbench::PartitionStatus amracut_status = dist_graph.PartitionAMRaCut(local_amracut_partition_labels, true);


  std::vector<uint64_t> local_sfc_partition_labels(local_element_count, my_rank);

  // collecting partitioning metrics (Only valid on MPI rank 0)

  std::vector<uint64_t> global_amracut_partition_sizes;
  std::vector<uint64_t> global_amracut_partition_boundaries;
  std::vector<uint64_t> global_amracut_partition_cuts;

  dist_graph.GetPartitionMetrics(local_amracut_partition_labels, global_amracut_partition_sizes,
                                 global_amracut_partition_boundaries, global_amracut_partition_cuts);

  std::vector<uint64_t> global_sfc_partition_sizes;
  std::vector<uint64_t> global_sfc_partition_boundaries;
  std::vector<uint64_t> global_sfc_partition_cuts;

  dist_graph.GetPartitionMetrics(local_sfc_partition_labels, global_sfc_partition_sizes,
                                 global_sfc_partition_boundaries, global_sfc_partition_cuts);

  if(!my_rank)
  {
    amracut_testbench::ExportMetricsToJson(mesh_file_path, procs_n, global_element_count, 
                                           SFC_status.time_us, global_sfc_partition_sizes, global_sfc_partition_boundaries, global_sfc_partition_cuts,
                                           amracut_status.time_us, global_amracut_partition_sizes, global_amracut_partition_boundaries, global_amracut_partition_cuts,
                                           output_file_path);
  }

#ifdef ENABLE_VTK_FEATURES
  std::vector<amracut_testbench::Element> global_sfc_elements;
  std::vector<uint64_t> global_sfc_partition_labels;
  std::vector<uint64_t> global_amracut_partition_labels;


  if (!my_rank)
  {
    amracut_testbench::print_log("starting VTK export");
    global_sfc_partition_labels.resize(global_element_count);
    global_amracut_partition_labels.resize(global_element_count);
    global_sfc_elements.resize(global_element_count);
  }

  MPI_Gatherv(local_sfc_elements.data(), local_element_count, par::Mpi_datatype<amracut_testbench::Element>::value(), global_sfc_elements.data(),
              proc_element_counts.data(),proc_element_counts_scanned.data(),par::Mpi_datatype<amracut_testbench::Element>::value(), 0, comm);
  MPI_Barrier(comm);
  
  MPI_Gatherv(local_sfc_partition_labels.data(), local_element_count, MPI_UINT64_T, global_sfc_partition_labels.data(),
              proc_element_counts.data(),proc_element_counts_scanned.data(),MPI_UINT64_T, 0, comm);

  MPI_Barrier(comm);
  
  MPI_Gatherv(local_amracut_partition_labels.data(), local_element_count, MPI_UINT64_T, global_amracut_partition_labels.data(),
              proc_element_counts.data(),proc_element_counts_scanned.data(),MPI_UINT64_T, 0, comm);

  if (!my_rank)
  {
    amracut_testbench::export_parititions_to_vtk(global_sfc_elements, 
                                                 global_sfc_partition_labels, global_amracut_partition_labels, 
                                                 "partitions");
    amracut_testbench::print_log("VTK export done");
                                                     
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