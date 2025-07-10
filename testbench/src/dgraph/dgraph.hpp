#ifndef AMRACUT_TESTBENCH_DGRAPH_H
#define AMRACUT_TESTBENCH_DGRAPH_H

#include <cinttypes>
#include <vector>

#include <mpi.h>

namespace amracut_testbench
{
  class DGraph
  {
  private:
    MPI_Comm comm;
    uint64_t global_count;
    uint64_t local_count;
    uint64_t ghost_count;
    uint64_t send_count;

    std::vector<uint64_t> local_degrees;
    std::vector<uint64_t> local_xdj;
    std::vector<uint64_t> local_adjncy;
    std::vector<uint64_t> dist_adjncy;

    std::vector<uint64_t> local_vertex_wgts;
    std::vector<uint64_t> dist_adjwgt;
    uint64_t wgt_flag;

    std::vector<int> vtx_dist;
    std::vector<int> vtx_counts;


    std::vector<int> ghost_counts;
    std::vector<int> ghost_counts_scanned;

    std::vector<uint64_t> sending_scatter_map;
    std::vector<int> send_counts;
    std::vector<int> send_counts_scanned;

    void ExtractGhosts(std::vector<std::pair<uint64_t, uint64_t>> &bdry_con_sorted,
                       std::vector<uint64_t> &unique_ghosts_out,
                       std::vector<int> &ghost_counts_per_proc_out);

  public:
    DGraph(
        const std::vector<std::pair<uint64_t, uint64_t>> &local_con,
        const std::vector<std::pair<uint64_t, uint64_t>> &bdry_con,
        const std::vector<int> &proc_element_counts,
        const std::vector<int> &proc_element_counts_scanned,
        const int wgt_flag,
        MPI_Comm comm);

    // ~DGraph();
    PartitionStatus PartitionAMRaCut(std::vector<uint64_t> &partition_labels_out, bool use_diffusion);

    void GetPartitionMetrics(std::vector<uint64_t> &local_partition_labels,
                             std::vector<uint64_t> &partition_sizes_out,
                             std::vector<uint64_t> &partition_boundaries_out,
                             std::vector<uint64_t> &partition_cuts_out);

  }; // class DGraph

} // namespace amracut_testbench

#endif