#include "stkLibiglBoolean3DMesher.h"

//---------VTK----------------------------------
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkSmartPointer.h>

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

//---------libigl----------------------------------
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/copyleft/cgal/piecewise_constant_winding_number.h>

vtkStandardNewMacro(stkLibiglBoolean3DMesher);

// ----------------------------------------------------------------------------
int stkLibiglBoolean3DMesher::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector)
{
  //  Get the input and output data objects.
  //  Get the info objects
  vtkPolyData* inputMeshA = this->GetInputMeshA();
  vtkPolyData* inputMeshB = this->GetInputMeshB();

  if (inputMeshA == nullptr)
  {
    vtkErrorMacro("Input Mesh A is empty.");
    return 0;
  }

  if (inputMeshA->GetPolys() == nullptr)
  {
    vtkErrorMacro("Input Mesh A does not contain any cell structure.");
    return 0;
  }

  if (inputMeshA->GetNumberOfCells() == 0)
  {
    vtkErrorMacro("Input Mesh A contains no cells.");
    return 0;
  }

  if (inputMeshA->GetPoints() == nullptr)
  {
    vtkErrorMacro("Input Mesh A does not contain any point structure.");
    return 0;
  }

  if (inputMeshA->GetNumberOfPoints() == 0)
  {
    vtkErrorMacro("Input Mesh A contains no points.");
    return 0;
  }

  if (inputMeshB == nullptr)
  {
    vtkErrorMacro("Input Mesh B is empty.");
    return 0;
  }

  if (inputMeshB->GetPolys() == nullptr)
  {
    vtkErrorMacro("Input Mesh B does not contain any cell structure.");
    return 0;
  }

  if (inputMeshB->GetNumberOfCells() == 0)
  {
    vtkErrorMacro("Input Mesh B contains no cells.");
    return 0;
  }

  if (inputMeshB->GetPoints() == nullptr)
  {
    vtkErrorMacro("Input Mesh B does not contain any point structure.");
    return 0;
  }

  if (inputMeshB->GetNumberOfPoints() == 0)
  {
    vtkErrorMacro("Input Mesh B contains no points.");
    return 0;
  }

  vtkPolyData* output0 = vtkPolyData::GetData(outputVector->GetInformationObject(0));

  vtkIdType numMeshVerts = inputMeshA->GetNumberOfPoints();
  Eigen::MatrixXd inputMeshAVerts(numMeshVerts, 3);

  vtkIdType numMeshCells = inputMeshA->GetNumberOfCells();
  Eigen::MatrixXi inputMeshACells(numMeshCells, 3);

  ConvertVTKMeshToEigenVerts(inputMeshA, inputMeshAVerts);
  ConvertVTKMeshToEigenCells(inputMeshA, inputMeshACells);

  numMeshVerts = inputMeshB->GetNumberOfPoints();
  Eigen::MatrixXd inputMeshBVerts(numMeshVerts, 3);

  numMeshCells = inputMeshB->GetNumberOfCells();
  Eigen::MatrixXi inputMeshBCells(numMeshCells, 3);

  ConvertVTKMeshToEigenVerts(inputMeshB, inputMeshBVerts);
  ConvertVTKMeshToEigenCells(inputMeshB, inputMeshBCells);

  // Preconditions
  if (!this->SkipPreconditions)
  {
    if (!igl::copyleft::cgal::piecewise_constant_winding_number(inputMeshAVerts, inputMeshACells))
    {
      vtkErrorMacro("Input mesh A is not PWN. See "
                    "http://www.cs.columbia.edu/cg/mesh-arrangements/#definitions.");
      return 0;
    }

    if (!igl::copyleft::cgal::piecewise_constant_winding_number(inputMeshBVerts, inputMeshBCells))
    {
      vtkErrorMacro("Input mesh B is not PWN. See "
                    "http://www.cs.columbia.edu/cg/mesh-arrangements/#definitions.");
      return 0;
    }
  }

  Eigen::MatrixXd outputMeshVerts;
  Eigen::MatrixXi outputMeshCells;
  try
  {
    if (this->Mode == stkLibiglBoolean3DMesher::Modes::UNION)
    {
      igl::copyleft::cgal::mesh_boolean(inputMeshAVerts, inputMeshACells, inputMeshBVerts,
        inputMeshBCells, igl::MeshBooleanType::MESH_BOOLEAN_TYPE_UNION, outputMeshVerts,
        outputMeshCells);
    }
    else if (this->Mode == stkLibiglBoolean3DMesher::Modes::INTERSECTION)
    {
      igl::copyleft::cgal::mesh_boolean(inputMeshAVerts, inputMeshACells, inputMeshBVerts,
        inputMeshBCells, igl::MeshBooleanType::MESH_BOOLEAN_TYPE_INTERSECT, outputMeshVerts,
        outputMeshCells);
    }
    else if (this->Mode == stkLibiglBoolean3DMesher::Modes::DIFFERENCE1)
    {
      igl::copyleft::cgal::mesh_boolean(inputMeshAVerts, inputMeshACells, inputMeshBVerts,
        inputMeshBCells, igl::MeshBooleanType::MESH_BOOLEAN_TYPE_MINUS, outputMeshVerts,
        outputMeshCells);
    }
    else if (this->Mode == stkLibiglBoolean3DMesher::Modes::DIFFERENCE2)
    {
      igl::copyleft::cgal::mesh_boolean(inputMeshBVerts, inputMeshBCells, inputMeshAVerts,
        inputMeshACells, igl::MeshBooleanType::MESH_BOOLEAN_TYPE_MINUS, outputMeshVerts,
        outputMeshCells);
    }
  }
  catch (const std::exception& e)
  {
    vtkErrorMacro(<< "Error caught : " << e.what());
    return 0;
  }

  ConvertEigenToVTKMesh(outputMeshVerts, outputMeshCells, output0);

  return 1;
}

//----------------------------------------------------------------------------
int stkLibiglBoolean3DMesher::ConvertVTKMeshToEigenVerts(
  vtkDataSet* object, Eigen::MatrixXd& vertices)
{
  auto numVerts = object->GetNumberOfPoints();
  assert(vertices.cols() == 3);

  double tmpPt[3];
  for (decltype(numVerts) i = 0; i < numVerts; ++i)
  {
    object->GetPoint(i, tmpPt);
    vertices(i, 0) = tmpPt[0];
    vertices(i, 1) = tmpPt[1];
    vertices(i, 2) = tmpPt[2];
  }
  return 1;
}

//----------------------------------------------------------------------------
int stkLibiglBoolean3DMesher::ConvertVTKMeshToEigenCells(vtkDataSet* object, Eigen::MatrixXi& cells)
{
  auto numCells = object->GetNumberOfCells();
  auto numRows = cells.rows(); // number of cells
  auto numCols = cells.cols(); // number of vertices per cell
  assert(numRows == numCells);

  vtkSmartPointer<vtkIdList> cell_pts_ids = vtkSmartPointer<vtkIdList>::New();
  for (decltype(numCells) i = 0; i < numCells; ++i)
  {
    object->GetCellPoints(i, cell_pts_ids);
    auto numIds = cell_pts_ids->GetNumberOfIds();
    assert(numIds == numCols);

    for (decltype(numCols) j = 0; j < numCols; ++j)
    {
      cells(i, j) = cell_pts_ids->GetId(j);
    }
  }

  return 1;
}

//----------------------------------------------------------------------------
int stkLibiglBoolean3DMesher::ConvertEigenToVTKMesh(
  const Eigen::MatrixXd& verts, const Eigen::MatrixXi& cells, vtkPointSet* outObject)
{
  auto numVRows = verts.rows(); // number of vertices
  auto numVCols = verts.cols();
  assert(numVCols == 3);

  auto numCRows = cells.rows(); // number of cells
  auto numCCols = cells.cols(); // number of vertices per cell

  vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
  newPoints->SetNumberOfPoints(numVRows);

  vtkSmartPointer<vtkCellArray> newCells = vtkSmartPointer<vtkCellArray>::New();

  double ptTmp[3];
  for (int i = 0; i < numVRows; ++i)
  {
    ptTmp[0] = verts(i, 0);
    ptTmp[1] = verts(i, 1);
    ptTmp[2] = verts(i, 2);
    newPoints->SetPoint(i, ptTmp);
  }

  for (int i = 0; i < numCRows; i++)
  {
    newCells->InsertNextCell(cells.cols());
    for (int j = 0; j < numCCols; ++j)
    {
      newCells->InsertCellPoint(cells(i, j));
    }
  }

  outObject->SetPoints(newPoints);

  if (outObject->IsA("vtkPolyData"))
  {
    vtkPolyData* castedOutObject = vtkPolyData::SafeDownCast(outObject);
    castedOutObject->SetPolys(newCells);
  }
  else if (outObject->IsA("vtkUnstructuredGrid"))
  {
    vtkUnstructuredGrid* castedOutObject = vtkUnstructuredGrid::SafeDownCast(outObject);
    if (numCCols == 3)
    {
      castedOutObject->SetCells(VTK_TRIANGLE, newCells);
    }
    else if (numCCols == 4)
    {
      castedOutObject->SetCells(VTK_TETRA, newCells);
    }
  }

  return 1;
}
