#include <iostream>
#include <string>

#include "mpi.h"

#include "dtypes/dtypes.hpp"
#include "mesh_utils/mesh_utils.hpp"
#include "utils/utils.hpp"


/**
 * steps
 *
 * 1.
 *
 */

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

  amracut_testbench::ElementType elementType = amracut_testbench::GetElementType(mesh_file_path, comm);

  amracut_testbench::print_log("[", my_rank, "]", elementType);

  MPI_Finalize();



  return 0;
}
