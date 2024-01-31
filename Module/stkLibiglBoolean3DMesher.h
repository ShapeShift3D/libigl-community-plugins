/**
 * @class stkLibiglBoolean3DMesher
 * @brief Applies a boolean operation between two 3D surfaces.
 *
 * This filter takes two inputs, inputMeshA and inputMeshB, of type vtkPolyData and applies
 * one of the four boolean operations (Union, Intersection, Difference1 (A - B) and Difference2 (B -
 * A)) between them.
 *
 * Inputs: Surface Mesh A (port 0, vtkPolyData), Surface Mesh B (port 1, vtkPolyData)
 * Output: Surface Mesh Result (port 0, vtkPolyData)
 *
 * @sa
 * stkLibiglBoolean3DMesher
 */
#pragma once

#include <stkLibiglBoolean3DMesherInterface.h>
#include <stkLibiglCopyleftModule.h>

#include <Eigen/Dense>

/**
 * @ingroup stkLibiglCopyleftModule
 *
 */
class STKLIBIGLCOPYLEFT_EXPORT stkLibiglBoolean3DMesher : public stkLibiglBoolean3DMesherInterface
{
public:
  static stkLibiglBoolean3DMesher* New();
  vtkTypeMacro(stkLibiglBoolean3DMesher, stkLibiglBoolean3DMesherInterface);

protected:
  stkLibiglBoolean3DMesher() = default;
  ~stkLibiglBoolean3DMesher() = default;

  virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
  stkLibiglBoolean3DMesher(const stkLibiglBoolean3DMesher&) = delete;
  void operator=(const stkLibiglBoolean3DMesher&) = delete;

  int ConvertVTKMeshToEigenVerts(vtkDataSet* object, Eigen::MatrixXd& vertices);
  int ConvertVTKMeshToEigenCells(vtkDataSet* object, Eigen::MatrixXi& cells);
  int ConvertEigenToVTKMesh(
    const Eigen::MatrixXd& verts, const Eigen::MatrixXi& cells, vtkPointSet* outObject);
};
