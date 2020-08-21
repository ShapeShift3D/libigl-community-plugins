#ifndef __stkLibiglCopyLeftScalarClip_h
#define __stkLibiglCopyLeftScalarClip_h

#include <stkLibiglCopyleftModule.h>
#include <vtkPolyDataAlgorithm.h>
#include "stkLibiglBoolean3DMesher.h"

#include <Eigen/dense>

/**
* @ingroup stkLibiglCopyleftModule
* 
*/
class STKLIBIGLCOPYLEFT_EXPORT stkLibiglCopyLeftScalarClip : public vtkPolyDataAlgorithm
{
public:
	static stkLibiglCopyLeftScalarClip* New();
	vtkTypeMacro(stkLibiglCopyLeftScalarClip, vtkPolyDataAlgorithm);
	
	vtkGetMacro(key, int);
	vtkSetMacro(key, int);
	
protected:
	stkLibiglCopyLeftScalarClip();
	~stkLibiglCopyLeftScalarClip() {};

	int RequestData(
	  vtkInformation* request,
	  vtkInformationVector** inputVector,
	  vtkInformationVector* outputVector
	) override;

private:
	stkLibiglCopyLeftScalarClip(const stkLibiglCopyLeftScalarClip&) =delete;
	void operator = (const stkLibiglCopyLeftScalarClip&) =delete;
	stkLibiglBoolean3DMesher lb3d;

	int key;
};

#endif // __stkLibiglCopyLeftScalarClip_h
