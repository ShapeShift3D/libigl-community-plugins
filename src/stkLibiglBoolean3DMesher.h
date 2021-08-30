#pragma once

#include <stkLibiglBoolean3DMesherInterface.h>
#include <stkLibiglCopyleftModule.h>

#include <Eigen/Dense>

// Inherit from the desired filter
class STKLIBIGLCOPYLEFT_EXPORT stkLibiglBoolean3DMesher : public stkLibiglBoolean3DMesherInterface
{
public:
  // VTK requirements
  static stkLibiglBoolean3DMesher* New();
  vtkTypeMacro(stkLibiglBoolean3DMesher, stkLibiglBoolean3DMesherInterface);

protected:
  stkLibiglBoolean3DMesher() = default;
  ~stkLibiglBoolean3DMesher() = default;

  virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
  // needed but not implemented
  stkLibiglBoolean3DMesher(const stkLibiglBoolean3DMesher&) = delete;
  void operator=(const stkLibiglBoolean3DMesher&) = delete;

  int ConvertVTKMeshToEigenVerts(vtkDataSet* object, Eigen::MatrixXd& vertices);
  int ConvertVTKMeshToEigenCells(vtkDataSet* object, Eigen::MatrixXi& cells);
  int ConvertEigenToVTKMesh(
    const Eigen::MatrixXd& verts, const Eigen::MatrixXi& cells, vtkPointSet* outObject);
};
