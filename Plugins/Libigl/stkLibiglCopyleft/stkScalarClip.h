#ifndef __stkScalarClip_h
#define __stkScalarClip_h

#include <stkLibiglCopyleftModule.h>
#include <vtkPolyDataAlgorithm.h>
#include "stkLibiglBoolean3DMesher.h"

#include <Eigen/dense>

/**
* @ingroup stkLibiglCopyleftModule
* 
*/
class STKLIBIGLCOPYLEFT_EXPORT stkScalarClip : public vtkPolyDataAlgorithm
{
public:
	static stkScalarClip* New();
	vtkTypeMacro(stkScalarClip, vtkPolyDataAlgorithm);
	
	vtkGetMacro(key, int);
	vtkSetMacro(key, int);
	
protected:
	stkScalarClip();
	~stkScalarClip() {};

	int RequestData(
	  vtkInformation* request,
	  vtkInformationVector** inputVector,
	  vtkInformationVector* outputVector
	) override;

private:
	stkScalarClip(const stkScalarClip&) =delete;
	void operator = (const stkScalarClip&) =delete;
	stkLibiglBoolean3DMesher lb3d;

	int key;
};

#endif // __stkScalarClip_h
