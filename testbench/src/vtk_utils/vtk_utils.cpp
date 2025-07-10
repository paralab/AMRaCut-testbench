#include <vector>
#include <string>
#include <cinttypes>
#include <stdexcept>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkVertex.h>
#include <vtkPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkXMLPolyDataWriter.h>

#include "dtypes/dtypes.hpp"

namespace amracut_testbench
{

  void export_parititions_to_vtk(
      const std::vector<amracut_testbench::Element> &elements,
      const std::vector<uint64_t> &sfc_partition_labels,
      const std::vector<uint64_t> &amracut_partition_labels,
      const std::string output_filename)
  {
    size_t n = elements.size();
    if (sfc_partition_labels.size() != n || amracut_partition_labels.size() != n)
    {
      throw std::runtime_error("Element and partition label vector size mismatch");
    }

    auto points = vtkSmartPointer<vtkPoints>::New();
    auto vertices = vtkSmartPointer<vtkCellArray>::New();
    auto polyData = vtkSmartPointer<vtkPolyData>::New();

    for (vtkIdType i = 0; i < static_cast<vtkIdType>(n); ++i)
    {
      points->InsertNextPoint(elements[i].x, elements[i].y, elements[i].z);
      vtkIdType pid = i;
      vertices->InsertNextCell(1, &pid); // each vertex is a cell
    }

    polyData->SetPoints(points);
    polyData->SetVerts(vertices); // Define vertices as cells

    auto array_sfc = vtkSmartPointer<vtkDoubleArray>::New();
    array_sfc->SetName("SFC");
    array_sfc->SetNumberOfComponents(1);
    array_sfc->SetNumberOfTuples(n);

    auto array_amracut = vtkSmartPointer<vtkDoubleArray>::New();
    array_amracut->SetName("AMRaCut");
    array_amracut->SetNumberOfComponents(1);
    array_amracut->SetNumberOfTuples(n);

    for (vtkIdType i = 0; i < static_cast<vtkIdType>(n); ++i)
    {
      array_sfc->SetValue(i, static_cast<double>(sfc_partition_labels[i]));
      array_amracut->SetValue(i, static_cast<double>(amracut_partition_labels[i]));
    }

    polyData->GetCellData()->AddArray(array_sfc);
    polyData->GetCellData()->AddArray(array_amracut);

    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName((output_filename + ".vtp").c_str());
    writer->SetInputData(polyData);
    writer->SetDataModeToBinary();     // Write in binary format
    writer->SetCompressorTypeToZLib(); // Optional: compress to reduce size
    writer->Write();
  }

} // namespace amracut_testbench