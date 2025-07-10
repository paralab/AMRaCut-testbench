#ifndef AMRACUT_TESTBENCH_DTYPES_H
#define AMRACUT_TESTBENCH_DTYPES_H

#include <cinttypes>
#include <iostream>

#include "usort/dtypes.h"


#define DIST_GRAPH_UNWEIGHTED           0
#define DIST_GRAPH_VTX_WEIGHTED         1
#define DIST_GRAPH_EDGE_WEIGHTED        2
#define DIST_GRAPH_VTX_EDGE_WEIGHTED    3

namespace amracut_testbench
{
  enum ElementType { TET=4, HEX=5 };

  struct TetElement
  {
    uint64_t element_tag;
    uint64_t global_idx; // populated after global morton sort
    double x;
    double y;
    double z;
    uint64_t morton_encoding;
    uint64_t face_tags[4];

    bool operator==(const TetElement &other) const
    {
      return morton_encoding == other.morton_encoding;
    }
    bool operator<(const TetElement &other) const
    {
      return morton_encoding < other.morton_encoding;
    }
    bool operator<=(const TetElement &other) const
    {
      return morton_encoding <= other.morton_encoding;
    }
    bool operator>=(const TetElement &other) const
    {
      return morton_encoding >= other.morton_encoding;
    }
    bool operator>(const TetElement &other) const
    {
      return morton_encoding > other.morton_encoding;
    }
  };
  // Overloading the << operator for TetElement
  std::ostream &operator<<(std::ostream &os, const TetElement &obj);

  struct HexElement
  {
    uint64_t element_tag;
    uint64_t global_idx; // populated after global morton sort
    double x;
    double y;
    double z;
    uint64_t morton_encoding;
    uint64_t face_tags[6];

    bool operator==(const HexElement &other) const
    {
      return morton_encoding == other.morton_encoding;
    }
    bool operator<(const HexElement &other) const
    {
      return morton_encoding < other.morton_encoding;
    }
    bool operator<=(const HexElement &other) const
    {
      return morton_encoding <= other.morton_encoding;
    }
    bool operator>=(const HexElement &other) const
    {
      return morton_encoding >= other.morton_encoding;
    }
    bool operator>(const HexElement &other) const
    {
      return morton_encoding > other.morton_encoding;
    }
  };
  // Overloading the << operator for HexElement
  std::ostream &operator<<(std::ostream &os, const HexElement &obj);

  struct ElementFace
  {
    uint64_t global_idx;
    uint64_t face_tag;
    bool operator==(const ElementFace &other) const
    {
      return face_tag == other.face_tag;
    }
    bool operator<(const ElementFace &other) const
    {
      return face_tag < other.face_tag;
    }
    bool operator<=(const ElementFace &other) const
    {
      return face_tag <= other.face_tag;
    }
    bool operator>=(const ElementFace &other) const
    {
      return face_tag >= other.face_tag;
    }
    bool operator>(const ElementFace &other) const
    {
      return face_tag > other.face_tag;
    }
  };

  struct Element
  {
    uint64_t global_idx;
    double x;
    double y;
    double z;
  };

  struct PartitionStatus
  {
    int return_code;
    uint64_t time_us;
  };

} // namespace amracut_testbench

namespace par
{
  template <>
  class Mpi_datatype<amracut_testbench::TetElement>
  {

    /**
    @return the MPI_Datatype for the C++ datatype "TetElement"
    **/
  public:
    static MPI_Datatype value()
    {
      static bool first = true;
      static MPI_Datatype custom_mpi_type;

      if (first)
      {
        first = false;
        int block_lengths[7] = {1, 1, 1, 1, 1, 1, 4};
        MPI_Datatype types[7] = {MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UINT64_T, MPI_UINT64_T};
        MPI_Aint offsets[7];
        offsets[0] = offsetof(amracut_testbench::TetElement, element_tag);
        offsets[1] = offsetof(amracut_testbench::TetElement, global_idx);
        offsets[2] = offsetof(amracut_testbench::TetElement, x);
        offsets[3] = offsetof(amracut_testbench::TetElement, y);
        offsets[4] = offsetof(amracut_testbench::TetElement, z);
        offsets[5] = offsetof(amracut_testbench::TetElement, morton_encoding);
        offsets[6] = offsetof(amracut_testbench::TetElement, face_tags);

        MPI_Type_create_struct(7, block_lengths, offsets, types, &custom_mpi_type);
        MPI_Type_commit(&custom_mpi_type);
      }

      return custom_mpi_type;
    }
  };

  template <>
  class Mpi_datatype<amracut_testbench::HexElement>
  {

    /**
    @return the MPI_Datatype for the C++ datatype "HexElement"
    **/
  public:
    static MPI_Datatype value()
    {
      static bool first = true;
      static MPI_Datatype custom_mpi_type;

      if (first)
      {
        first = false;
        int block_lengths[7] = {1, 1, 1, 1, 1, 1, 6};
        MPI_Datatype types[7] = {MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UINT64_T, MPI_UINT64_T};
        MPI_Aint offsets[7];
        offsets[0] = offsetof(amracut_testbench::HexElement, element_tag);
        offsets[1] = offsetof(amracut_testbench::HexElement, global_idx);
        offsets[2] = offsetof(amracut_testbench::HexElement, x);
        offsets[3] = offsetof(amracut_testbench::HexElement, y);
        offsets[4] = offsetof(amracut_testbench::HexElement, z);
        offsets[5] = offsetof(amracut_testbench::HexElement, morton_encoding);
        offsets[6] = offsetof(amracut_testbench::HexElement, face_tags);

        MPI_Type_create_struct(7, block_lengths, offsets, types, &custom_mpi_type);
        MPI_Type_commit(&custom_mpi_type);
      }

      return custom_mpi_type;
    }
  };

  template <>
  class Mpi_datatype<amracut_testbench::Element>
  {

  /**
  @return the MPI_Datatype for the C++ datatype "Element"
  **/
  public:
    static MPI_Datatype value()
    {
      static bool first = true;
      static MPI_Datatype custom_mpi_type;

      if (first)
      {
        first = false;
        int block_lengths[4] = {1, 1, 1, 1};
        MPI_Datatype types[4] = {MPI_UINT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
        MPI_Aint offsets[4];
        offsets[0] = offsetof(amracut_testbench::Element, global_idx);
        offsets[1] = offsetof(amracut_testbench::Element, x);
        offsets[2] = offsetof(amracut_testbench::Element, y);
        offsets[3] = offsetof(amracut_testbench::Element, z);

        MPI_Type_create_struct(4, block_lengths, offsets, types, &custom_mpi_type);
        MPI_Type_commit(&custom_mpi_type);
      }

      return custom_mpi_type;
    }
  };

  template <>
  class Mpi_datatype<amracut_testbench::ElementFace>
  {

  /**
  @return the MPI_Datatype for the C++ datatype "ElementFace"
  **/
  public:
    static MPI_Datatype value()
    {
      static bool first = true;
      static MPI_Datatype custom_mpi_type;

      if (first)
      {
        first = false;
        int block_lengths[2] = {1, 1};
        MPI_Datatype types[2] = {MPI_UINT64_T, MPI_UINT64_T};
        MPI_Aint offsets[2];
        offsets[0] = offsetof(amracut_testbench::ElementFace, global_idx);
        offsets[1] = offsetof(amracut_testbench::ElementFace, face_tag);

        MPI_Type_create_struct(2, block_lengths, offsets, types, &custom_mpi_type);
        MPI_Type_commit(&custom_mpi_type);
      }

      return custom_mpi_type;
    }
  };

  template <>
  class Mpi_pairtype<uint64_t, uint64_t>
  {

  /**
  @return the MPI_Datatype for the C++ datatype "std::pair<uint64_t,uint64_t>"
  **/
  public:
    static MPI_Datatype value()
    {
      static bool first = true;
      static MPI_Datatype custom_mpi_type;

      if (first)
      {

        first = false;
        MPI_Datatype inner_type = MPI_UINT64_T;

        int second_value_offset;
        MPI_Type_size(inner_type, &second_value_offset);
        int block_lengths[2] = {1, 1};
        MPI_Datatype types[2] = {inner_type, inner_type};
        MPI_Aint offsets[2];
        offsets[0] = 0;
        offsets[1] = static_cast<MPI_Aint>(second_value_offset);

        MPI_Type_create_struct(2, block_lengths, offsets, types, &custom_mpi_type);
        MPI_Type_commit(&custom_mpi_type);
      }

      return custom_mpi_type;
    }
  };

} // namespace par

#endif