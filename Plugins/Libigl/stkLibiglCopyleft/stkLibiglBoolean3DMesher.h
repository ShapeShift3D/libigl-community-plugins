#ifndef stkLibiglBoolean3DMesher_h
#define stkLibiglBoolean3DMesher_h
// Gives access to macros for communication with the UI
#include "vtkPolyDataAlgorithm.h" 
#include <stkLibiglCopyleftModule.h>

#include <Eigen/Dense>

// Inherit from the desired filter
class STKLIBIGLCOPYLEFT_EXPORT stkLibiglBoolean3DMesher : public vtkPolyDataAlgorithm
{
public:
  // VTK requirements
  static stkLibiglBoolean3DMesher* New();
  vtkTypeMacro(stkLibiglBoolean3DMesher, vtkPolyDataAlgorithm);
  
  vtkPolyData* GetInputMeshA();
  vtkPolyData* GetInputMeshB();

  enum Modes {
      UNION = 0,
      INTERSECTION,
      DIFFERENCE,
      DIFFERENCE2
  };

  vtkGetMacro(Mode, int);
  vtkSetMacro(Mode, int);

  void SetModeToUnion();
  void SetModeToIntersection();
  void SetModeToDifference();
  void SetModeToDifference2();

  //@{
  /**
  * Skip preconditions if true
  */
  vtkGetMacro(SkipPreconditions, bool);
  vtkSetMacro(SkipPreconditions, bool);
  vtkBooleanMacro(SkipPreconditions, bool);
  //@}

protected:
  stkLibiglBoolean3DMesher();
  ~stkLibiglBoolean3DMesher(){}

  int Mode;

  virtual int RequestData(vtkInformation*, 
                            vtkInformationVector**, 
                            vtkInformationVector*) override;

  bool SkipPreconditions;

private:
  // needed but not implemented
  stkLibiglBoolean3DMesher(const stkLibiglBoolean3DMesher&) = delete;
  void operator=(const stkLibiglBoolean3DMesher&) = delete;

  int ConvertVTKMeshToEigenVerts(vtkDataSet* object, Eigen::MatrixXd& vertices);
  int ConvertVTKMeshToEigenCells(vtkDataSet* object, Eigen::MatrixXi& cells);
  int ConvertEigenToVTKMesh(const Eigen::MatrixXd& verts, const Eigen::MatrixXi& cells,
      vtkPointSet* outObject);
};
#endif
