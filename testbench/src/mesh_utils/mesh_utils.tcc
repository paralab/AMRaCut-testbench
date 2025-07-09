#include <string>
#include <vector>
#include <map>
#include <cassert>
#include "mpi.h"
#include "gmsh.h"

#include "dtypes/dtypes.hpp"
#include "usort/ompUtils.h"
#include "usort/dtypes.h"
#include "utils/utils.hpp"



namespace amracut_testbench
{

  template <class T>
  void GetInitialElementsDistribution(const std::string &mesh_file_path,
                                      std::vector<T> &local_elements_out,
                                      ElementType element_type, MPI_Comm comm)
  {
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);

    int local_element_count;
    std::vector<int> proc_element_counts(procs_n);         // populated only at root
    std::vector<int> proc_element_counts_scanned(procs_n); // populated only at root
    std::vector<T> all_elements;                           // populated only at root

    // root process will read the mesh and do a simple initial element scatter distribution
    if (!my_rank)
    {

      gmsh::initialize();
      // gmsh::option::setNumber("General.Terminal", 0);
      gmsh::open(mesh_file_path);

      int gmsh_element_type;
      int gmsh_face_type;
      int faces_per_element;
      int nodes_per_element;
      int nodes_per_face;

      switch (element_type)
      {
      case ElementType::TET: // linear tet
      {
        gmsh_element_type = 4;
        gmsh_face_type = 3; // triangle
        nodes_per_face = 3;
        faces_per_element = 4;
        nodes_per_element = 4;

        break;
      }
      case 5: // linear hexahedra
      {
        gmsh_element_type = 5;
        gmsh_face_type = 4; // quadtriangle
        nodes_per_face = 4;
        faces_per_element = 6;
        nodes_per_element = 8;

        break;
      }

      default:
      {
        throw std::invalid_argument("unknown element type\t exiting...");
        break;
      }
      }

      std::vector<uint64_t> element_node_tags;
      std::vector<uint64_t> element_tags;

      gmsh::model::mesh::getElementsByType(gmsh_element_type, element_tags, element_node_tags, -1);

      // print_log("[", my_rank, "] element_tags:", VectorToString(element_tags));

      size_t total_element_count = element_tags.size();
      assert(element_node_tags.size() == total_element_count * nodes_per_element);

      std::vector<uint64_t> allNodeTags;
      std::vector<double> allNodeCoords, allNodeParams;

      gmsh::model::mesh::getNodes(allNodeTags, allNodeCoords, allNodeParams, -1, -1, true, false);

      std::map<uint64_t, uint64_t> node_tag_to_index;
      for (size_t node_i = 0; node_i < allNodeTags.size(); node_i++)
      {
        node_tag_to_index[allNodeTags[node_i]] = node_i;
      }

      std::vector<double> element_coordinates(total_element_count * 3); // 3 for 3D - x,y,z

      for (size_t element_i = 0; element_i < total_element_count; element_i++)
      {
        double x = 0;
        double y = 0;
        double z = 0;

        for (size_t elem_node_i = 0; elem_node_i < nodes_per_element; elem_node_i++)
        {
          size_t nodeTag = element_node_tags[element_i * nodes_per_element + elem_node_i];

          x += allNodeCoords[node_tag_to_index[nodeTag] * 3];
          y += allNodeCoords[node_tag_to_index[nodeTag] * 3 + 1];
          z += allNodeCoords[node_tag_to_index[nodeTag] * 3 + 2];
        }
        x = x / nodes_per_element;
        y = y / nodes_per_element;
        z = z / nodes_per_element;

        element_coordinates[element_i * 3] = x;
        element_coordinates[element_i * 3 + 1] = y;
        element_coordinates[element_i * 3 + 2] = z;
      }

      std::vector<std::size_t> faceNodes;
      gmsh::model::mesh::getElementFaceNodes(gmsh_element_type, gmsh_face_type, faceNodes, -1, false);

      assert(faceNodes.size() == total_element_count * faces_per_element * nodes_per_face);

      gmsh::model::mesh::createFaces();

      std::vector<std::size_t> faceTags;
      std::vector<int> faceOrientations;
      gmsh::model::mesh::getFaces(gmsh_face_type, faceNodes, faceTags, faceOrientations);
      assert(faceTags.size() == (faces_per_element * total_element_count));
      gmsh::finalize();

      all_elements.resize(total_element_count);

      for (size_t element_i = 0; element_i < total_element_count; element_i++)
      {
        all_elements[element_i].element_tag = element_tags[element_i];
        all_elements[element_i].x = element_coordinates[element_i * 3];
        all_elements[element_i].y = element_coordinates[element_i * 3 + 1];
        all_elements[element_i].z = element_coordinates[element_i * 3 + 2];

        // for (size_t elem_node_i = 0; elem_node_i < nodes_per_element; elem_node_i++)
        // {
        //     size_t nodeTag = element_node_tags[element_i*nodes_per_element + elem_node_i];
        //     all_elements[element_i].node_tags[elem_node_i] = nodeTag;
        // }

        for (size_t elem_face_i = 0; elem_face_i < faces_per_element; elem_face_i++)
        {
          all_elements[element_i].face_tags[elem_face_i] = faceTags[element_i * faces_per_element + elem_face_i];
        }
      }
      std::fill(proc_element_counts.begin(), proc_element_counts.end(), total_element_count / procs_n);

      // Distribute the remainder element counts
      int rem = total_element_count % procs_n;
      for (int proc_i = 0; proc_i < rem; proc_i++)
      {
        proc_element_counts[proc_i]++;
      }
      omp_par::scan(&proc_element_counts[0], &proc_element_counts_scanned[0], procs_n);
      print_log("mesh reading done");
    } // mesh reading done

    MPI_Scatter(proc_element_counts.data(), 1, MPI_INT, &local_element_count, 1, MPI_INT, 0, comm);
    local_elements_out.resize(local_element_count);
    MPI_Scatterv(all_elements.data(),
                 proc_element_counts.data(), proc_element_counts_scanned.data(), par::Mpi_datatype<T>::value(),
                 local_elements_out.data(), local_element_count, par::Mpi_datatype<T>::value(),
                 0, comm);

    // print_log("[", my_rank, "] elements:", VectorToString(local_elements_out));
  }
}