
#include <vector>
#include <cinttypes>
#include <algorithm>
#include <chrono>
#include <mpi.h>

#include "amracut.h"

#include "utils/utils.hpp"
#include "dtypes/dtypes.hpp"
#include "dgraph/dgraph.hpp"
#include "usort/ompUtils.h"


namespace amracut_testbench
{
  DGraph::DGraph(
      const std::vector<std::pair<uint64_t, uint64_t>> &local_con,
      const std::vector<std::pair<uint64_t, uint64_t>> &bdry_con,
      const std::vector<int> &proc_element_counts,
      const std::vector<int> &proc_element_counts_scanned,
      const int wgt_flag,
      MPI_Comm comm)

  {
    this->comm = comm;
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);

    this->vtx_counts.assign(proc_element_counts.begin(), proc_element_counts.end());
    this->vtx_dist.assign(proc_element_counts_scanned.begin(), proc_element_counts_scanned.end());


    this->vtx_dist.push_back(proc_element_counts_scanned[procs_n-1] + proc_element_counts[procs_n-1]);
    this->global_count = proc_element_counts_scanned[procs_n-1] + proc_element_counts[procs_n-1];

    this->local_count = proc_element_counts[my_rank];

    this->wgt_flag = wgt_flag;

    std::vector<std::pair<uint64_t, uint64_t>> bdry_con_sorted(bdry_con);
    std::sort(bdry_con_sorted.begin(), bdry_con_sorted.end(), [](const auto & a, const auto & b) {
      return a.second < b.second;
    });

    std::vector<uint64_t> unique_ghosts;
    std::vector<int> ghost_counts_per_proc(procs_n, 0);

    this->ExtractGhosts(bdry_con_sorted, unique_ghosts, ghost_counts_per_proc);

    this->ghost_count = unique_ghosts.size();

    this->ghost_counts = ghost_counts_per_proc;
    this->ghost_counts_scanned.resize(procs_n);
    omp_par::scan(&this->ghost_counts[0], &this->ghost_counts_scanned[0], procs_n);


    this->local_degrees.resize(this->local_count + this->ghost_count);
    std::fill(this->local_degrees.begin(), this->local_degrees.end(), 0);
    this->local_xdj.resize(this->local_count + this->ghost_count + 1);
    std::fill(this->local_xdj.begin(), this->local_xdj.end(), 0);

    this->local_adjncy.resize(2 * (local_con.size() + bdry_con.size()));
    this->dist_adjncy.resize(2*local_con.size()  + bdry_con.size());
    // this->dist_adjwgt.resize(2*local_con.size()  + boundary_connectivity.size());

    for (auto& edge : local_con) {
        this->local_degrees[edge.first - proc_element_counts_scanned[my_rank]]++;
        this->local_degrees[edge.second - proc_element_counts_scanned[my_rank]]++;
    }


    /**
     * Calculating graph degree for local and ghost elements
     */
    if (!bdry_con_sorted.empty())
    {
      uint64_t ghost_idx = this->local_count;
      // assumes first element belongs ro this process and second element belongs to another process
      this->local_degrees[bdry_con_sorted[0].first - proc_element_counts_scanned[my_rank]]++;
      this->local_degrees[ghost_idx]++;

      for (size_t boundary_edge_i = 1; boundary_edge_i < bdry_con_sorted.size(); boundary_edge_i++)
      {
        // since we sorted boundary edges by second element
        if (bdry_con_sorted[boundary_edge_i].second != bdry_con_sorted[boundary_edge_i - 1].second)
        {
          ghost_idx++;
        }

        // assumes first element belongs ro this process and second element belongs to another process
        this->local_degrees[bdry_con_sorted[boundary_edge_i].first - proc_element_counts_scanned[my_rank]]++;
        this->local_degrees[ghost_idx]++;
      }
    }

    omp_par::scan(&this->local_degrees[0], &this->local_xdj[0], this->local_count + this->ghost_count);
    //scan operation does not populate the last extra entry for CSR encoding, hence we have to manually populate the last element
    this->local_xdj[this->local_count + this->ghost_count] = 
            this->local_xdj[this->local_count + this->ghost_count -1 ] + this->local_degrees[this->local_count + this->ghost_count-1];

    /**
     * populating adjacency structure
     */
    std::vector<uint64_t> next_index;
    next_index.assign(this->local_xdj.begin(), this->local_xdj.end() - 1);

    for (size_t local_edge_i = 0; local_edge_i < local_con.size(); local_edge_i++)
    {
      auto edge = local_con[local_edge_i];
      auto local_index_1 = edge.first - proc_element_counts_scanned[my_rank];
      auto local_index_2 = edge.second - proc_element_counts_scanned[my_rank];

      this->local_adjncy[next_index[local_index_1]] = local_index_2;
      this->local_adjncy[next_index[local_index_2]] = local_index_1;

      this->dist_adjncy[next_index[local_index_1]] = edge.second;
      this->dist_adjncy[next_index[local_index_2]] = edge.first;

      // this->dist_adjwgt[next_index[local_index_1]] = local_conectivity_weights[local_edge_i];
      // this->dist_adjwgt[next_index[local_index_2]] = local_conectivity_weights[local_edge_i];

      next_index[local_index_1]++;
      next_index[local_index_2]++;
    }

    if (!bdry_con_sorted.empty())
    {
      uint64_t ghost_idx = this->local_count;
      // assumes first element belongs to this process and second element belongs to another process
      {
        auto local_index_1 = bdry_con_sorted[0].first - proc_element_counts_scanned[my_rank];
        auto local_index_2 = ghost_idx;
        this->local_adjncy[next_index[local_index_1]] = local_index_2;
        this->local_adjncy[next_index[local_index_2]] = local_index_1;

        this->dist_adjncy[next_index[local_index_1]] = bdry_con_sorted[0].second;
        // this->dist_adjwgt[next_index[local_index_1]] = boundary_connectivity_weights[0];

        next_index[local_index_1]++;
        next_index[local_index_2]++;
      }

      for (size_t boundary_edge_i = 1; boundary_edge_i < bdry_con_sorted.size(); boundary_edge_i++)
      {
        // since we sorted boundary edges by second element
        if (bdry_con_sorted[boundary_edge_i].second != bdry_con_sorted[boundary_edge_i - 1].second)
        {
          ghost_idx++;
        }

        // assumes first element belongs ro this process and second element belongs to another process
        auto local_index_1 = bdry_con_sorted[boundary_edge_i].first - proc_element_counts_scanned[my_rank];
        auto local_index_2 = ghost_idx;
        this->local_adjncy[next_index[local_index_1]] = local_index_2;
        this->local_adjncy[next_index[local_index_2]] = local_index_1;

        this->dist_adjncy[next_index[local_index_1]] = bdry_con_sorted[boundary_edge_i].second;
        // this->dist_adjwgt[next_index[local_index_1]] = boundary_connectivity_weights[boundary_edge_i];

        next_index[local_index_1]++;
        next_index[local_index_2]++;
      }
    }

    /**
     * building the scatter map
     * */
    this->sending_scatter_map.clear();
    this->send_counts.resize(procs_n);
    std::fill(this->send_counts.begin(), this->send_counts.end(), 0);
    this->send_counts_scanned.resize(procs_n);
    std::fill(this->send_counts_scanned.begin(), this->send_counts_scanned.end(), 0);
    if (!bdry_con_sorted.empty())
    {
      uint64_t current_other_proc = 0;

      std::vector<uint64_t> send_elements_with_dups;
      std::vector<uint64_t> send_elements_with_dups_counts(procs_n, 0);
      std::vector<uint64_t> send_elements_with_dups_counts_scanned(procs_n, 0);

      for (auto &edge : bdry_con_sorted)
      {
        if (edge.second >= proc_element_counts_scanned[current_other_proc] &&
            edge.second < proc_element_counts_scanned[current_other_proc] + proc_element_counts[current_other_proc])
        {
          send_elements_with_dups.push_back(edge.first);
          send_elements_with_dups_counts[current_other_proc]++;
        }
        else
        {
          while (1)
          {
            current_other_proc++;
            if (edge.second >= proc_element_counts_scanned[current_other_proc] &&
                edge.second < proc_element_counts_scanned[current_other_proc] + proc_element_counts[current_other_proc])
            {
              send_elements_with_dups.push_back(edge.first);
              send_elements_with_dups_counts[current_other_proc]++;
              break;
            }
          }
        }
      }
      omp_par::scan(&send_elements_with_dups_counts[0], &send_elements_with_dups_counts_scanned[0], procs_n);

      std::vector<uint64_t> send_elements;

      // sort each send elements per each process
      for (int proc_i = 0; proc_i < procs_n; proc_i++)
      {
        if (send_elements_with_dups_counts[proc_i])
        {
          omp_par::merge_sort(&send_elements_with_dups[send_elements_with_dups_counts_scanned[proc_i]],
                              &send_elements_with_dups[send_elements_with_dups_counts_scanned[proc_i] + send_elements_with_dups_counts[proc_i]]);

          send_elements.push_back(send_elements_with_dups[send_elements_with_dups_counts_scanned[proc_i]]);
          send_counts[proc_i]++;
          for (size_t elem_i = send_elements_with_dups_counts_scanned[proc_i] + 1; elem_i < (send_elements_with_dups_counts_scanned[proc_i] + send_elements_with_dups_counts[proc_i]); elem_i++)
          {
            if (send_elements_with_dups[elem_i] != send_elements_with_dups[elem_i - 1])
            {
              send_elements.push_back(send_elements_with_dups[elem_i]);
              send_counts[proc_i]++;
            }
          }
        }
      }

      omp_par::scan(&this->send_counts[0], &this->send_counts_scanned[0], procs_n);
      this->send_count = this->send_counts_scanned[procs_n - 1] + this->send_counts[procs_n - 1];
      this->sending_scatter_map.resize(this->send_count);


      for (size_t send_elem_i = 0; send_elem_i < send_elements.size(); send_elem_i++)
      {
        this->sending_scatter_map[send_elem_i] = send_elements[send_elem_i] - proc_element_counts_scanned[my_rank];
      }
    }

  }   // DGraph::DGraph




  void DGraph::ExtractGhosts(std::vector<std::pair<uint64_t, uint64_t>> &bdry_con_sorted,
                             std::vector<uint64_t> &unique_ghosts_out,
                             std::vector<int> &ghost_counts_per_proc_out)

  {
    unique_ghosts_out.clear();

    if (bdry_con_sorted.size() > 0)
    {
      unique_ghosts_out.push_back(bdry_con_sorted[0].second);
    }

    for (size_t bdry_edge_i = 1; bdry_edge_i < bdry_con_sorted.size(); bdry_edge_i++)
    {
      if (bdry_con_sorted[bdry_edge_i].second !=
          bdry_con_sorted[bdry_edge_i - 1].second)
      {
        unique_ghosts_out.push_back(bdry_con_sorted[bdry_edge_i].second);
      }
    }

    std::fill(ghost_counts_per_proc_out.begin(), ghost_counts_per_proc_out.end(), 0);

    uint64_t current_proc = 0;
    for (size_t g_i = 0; g_i < unique_ghosts_out.size(); g_i++)
    {
      if (unique_ghosts_out[g_i] >= this->vtx_dist[current_proc] &&
          unique_ghosts_out[g_i] <
              (this->vtx_dist[current_proc] + this->vtx_counts[current_proc]))
      {
        ghost_counts_per_proc_out[current_proc]++;
      }
      else
      {
        while (1)
        {
          current_proc++;
          if (unique_ghosts_out[g_i] >= this->vtx_dist[current_proc] &&
              unique_ghosts_out[g_i] <
                  (this->vtx_dist[current_proc] + this->vtx_counts[current_proc]))
          {
            ghost_counts_per_proc_out[current_proc]++;
            break;
          }
        }
      }
    }
  }

  PartitionStatus DGraph::PartitionAMRaCut(std::vector<uint64_t> &partition_labels_out, bool use_diffusion)
  {
    int procs_n, my_rank;
    MPI_Comm_size(this->comm, &procs_n);
    MPI_Comm_rank(this->comm, &my_rank);

    std::vector<amracut_uint_t> xadj__(this->local_xdj.begin(), this->local_xdj.begin() + (this->local_count + 1));
    std::vector<amracut_uint_t> vtxdist__(this->vtx_dist.begin(), this->vtx_dist.end());
    std::vector<amracut_uint_t> adjncy__(this->dist_adjncy.begin(), this->dist_adjncy.end());
    std::vector<amracut_uint_t> partitions_labels(this->local_count);
    // std::vector<amracut_uint_t> local_vertex_wgts__(this->local_vertex_wgts.begin(), this->local_vertex_wgts.end());
    // std::vector<amracut_uint_t> adjwgt__(this->dist_adjwgt.begin(), this->dist_adjwgt.end());

    amracut_ctrl ctrl;
    amracut_setup(&ctrl, vtxdist__.data(), xadj__.data(), adjncy__.data(), NULL, NULL, this->wgt_flag, &(this->comm));
    MPI_Barrier(comm);
    auto start__ = std::chrono::high_resolution_clock::now();
    amracut_partgraph(&ctrl, partitions_labels.data(), use_diffusion, 0);
    MPI_Barrier(comm);
    auto end__ = std::chrono::high_resolution_clock::now();
    auto duration__ = std::chrono::duration_cast<std::chrono::microseconds>(end__ - start__);
    amracut_destroy(&ctrl);

    if (!my_rank)
    {
      print_log("AMRaCut total time:\t", duration__.count(), " us");
    }
    partition_labels_out.assign(partitions_labels.begin(), partitions_labels.end());
    return {.return_code = 0, .time_us = static_cast<uint64_t>(duration__.count())};
  }

} // namespace amracut_testbench