#include <string>
#include <stdexcept>
#include <vector>
#include "gmsh.h"

#include "dtypes/dtypes.hpp"


namespace amracut_testbench
{
  ElementType GetElementType(const std::string &mesh_file_path, MPI_Comm comm)
  {
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);

    int gmsh_element_type;

    // only the root will check and read the mesh element type. no point of every process reading the file.
    if (!my_rank)
    {
      gmsh::initialize();
      gmsh::option::setNumber("General.Verbosity", 0);
      gmsh::open(mesh_file_path);

      std::vector<std::pair<int, int>> dimTags;
      gmsh::model::getEntities(dimTags);

      std::vector<int> elementTypes;
      gmsh::model::mesh::getElementTypes(elementTypes, 3);

      gmsh::finalize();

      if (elementTypes.size() == 0)
      {
        throw std::invalid_argument("no 3D elements\t exiting...");
      }

      if (elementTypes.size() > 1)
      {
        throw std::invalid_argument("more than 1 element type\t exiting...");
      }
      gmsh_element_type = elementTypes[0];
    }

    MPI_Bcast(&gmsh_element_type, 1, MPI_INT, 0, comm);

    ElementType element_type;

    switch (gmsh_element_type)
    {
    case 4: // linear tet
    {
      element_type = ElementType::TET;
      // std::cout << "linear tetrahedra mesh\n";
      break;
    }
    case 5: // linear hexahedra
    {
      element_type = ElementType::HEX;
      // std::cout << "linear hexahedra mesh\n";
      break;
    }

    default:
    {
      throw std::invalid_argument("unknown element type\t exiting...");
      break;
    }
    }

    return element_type;
  }
} // namespace amracut_testbench
