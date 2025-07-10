#include <cinttypes>
#include <vector>
#include <mpi.h>

namespace amracut_testbench
{

  template <class T>
  void SetMortonEncoding(std::vector<T> &elements, ElementType element_type, MPI_Comm comm)
  {
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);

    double bounding_box_local_min = elements[0].x;
    double bounding_box_local_max = elements[0].x;
    for (size_t i = 0; i < elements.size(); i++)
    {
      bounding_box_local_min = std::min({bounding_box_local_min, elements[i].x, elements[i].y, elements[i].z});
      bounding_box_local_max = std::max({bounding_box_local_max, elements[i].x, elements[i].y, elements[i].z});
    }
    double bounding_box_global_min;
    double bounding_box_global_max;
    MPI_Allreduce(&bounding_box_local_min, &bounding_box_global_min, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&bounding_box_local_max, &bounding_box_global_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    uint64_t levels = 20;
    double spacing = (bounding_box_global_max - bounding_box_global_min) / ((double)(1 << levels));
    // print_log("[", my_rank, "]:", "bounding_box_global_min = ", bounding_box_global_min);
    // print_log("[", my_rank, "]:", "bounding_box_global_max = ", bounding_box_global_max);

    for (size_t element_i = 0; element_i < elements.size(); element_i++)
    {
      uint64_t x = (uint64_t)((elements[element_i].x - bounding_box_global_min) / spacing);
      uint64_t y = (uint64_t)((elements[element_i].y - bounding_box_global_min) / spacing);
      uint64_t z = (uint64_t)((elements[element_i].z - bounding_box_global_min) / spacing);

      elements[element_i].morton_encoding = mortonEncode_magicbits(x, y, z);
    }


  }

}
