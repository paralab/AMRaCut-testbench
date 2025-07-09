#ifndef AMRACUT_TESTBENCH_DTYPES_H
#define AMRACUT_TESTBENCH_DTYPES_H

#include <cinttypes>
#include <iostream>

#include "usort/dtypes.h"

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
    uint64_t element_tag;
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

} // namespace par

#endif