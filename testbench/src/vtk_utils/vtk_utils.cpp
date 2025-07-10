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
      // const std::vector<uint64_t> &p2,
      // const std::vector<uint64_t> &p3,
      const std::string output_filename)
  {
    size_t n = elements.size();
    if (sfc_partition_labels.size() != n /*|| p2.size() != n || p3.size() != n*/)
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

    // auto array_p2 = vtkSmartPointer<vtkDoubleArray>::New();
    // array_p2->SetName("p2");
    // array_p2->SetNumberOfComponents(1);
    // array_p2->SetNumberOfTuples(n);

    // auto array_p3 = vtkSmartPointer<vtkDoubleArray>::New();
    // array_p3->SetName("p3");
    // array_p3->SetNumberOfComponents(1);
    // array_p3->SetNumberOfTuples(n);

    for (vtkIdType i = 0; i < static_cast<vtkIdType>(n); ++i)
    {
      array_sfc->SetValue(i, static_cast<double>(sfc_partition_labels[i]));
      // array_p2->SetValue(i, static_cast<double>(p2[i]));
      // array_p3->SetValue(i, static_cast<double>(p3[i]));
    }

    polyData->GetCellData()->AddArray(array_sfc);
    // polyData->GetCellData()->AddArray(array_p2);
    // polyData->GetCellData()->AddArray(array_p3);

    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName((output_filename + ".vtp").c_str());
    writer->SetInputData(polyData);
    writer->SetDataModeToBinary();     // Write in binary format
    writer->SetCompressorTypeToZLib(); // Optional: compress to reduce size
    writer->Write();
  }

} // namespace amracut_testbench