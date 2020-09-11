#include "stkLibiglCopyLeftScalarClip.h"
#include "stkLibiglBoolean3DMesher.h"
#include "vtkCellData.h"
#include <vtkCellArray.h>
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>

#include <Eigen/sparse>

#include <igl/copyleft/marching_cubes.h>
#include <igl/signed_distance.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Core>

vtkStandardNewMacro(stkLibiglCopyLeftScalarClip);

//----------------------------------------------------------------------------
stkLibiglCopyLeftScalarClip::stkLibiglCopyLeftScalarClip()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

}

//----------------------------------------------------------------------------
int stkLibiglCopyLeftScalarClip::RequestData(
	vtkInformation* vtkNotUsed(request),
	vtkInformationVector** inputVector,
	vtkInformationVector* outputVector)
{	


	vtkInformation* mesh_info = inputVector[0]->GetInformationObject(0);
	vtkPolyData* mesh = vtkPolyData::SafeDownCast(mesh_info->Get(vtkDataObject::DATA_OBJECT()));

	vtkPolyData* outputData = vtkPolyData::GetData(outputVector, 0);

	if (mesh == nullptr)
	{
		vtkErrorMacro("Mesh is empty.");
		return 0;
	}

	if (mesh->GetNumberOfPoints() == 0)
	{
		vtkErrorMacro("Mesh is empty.");
		return 0;
	}
	// -------------------------------

	Eigen::MatrixXi F;
	Eigen::MatrixXd V;

	// Construct V matrix
	auto numberOfMeshVerts = mesh->GetNumberOfPoints();
	Eigen::MatrixXd meshVerts(numberOfMeshVerts, 3);
	lb3d.ConvertVTKMeshToEigenVerts(mesh, meshVerts);

	V = meshVerts;

	// Construct F matrix
	auto numberOfMeshCells = mesh->GetNumberOfCells();
	Eigen::MatrixXi meshCells(numberOfMeshCells, 3);
	lb3d.ConvertVTKMeshToEigenCells(mesh, meshCells);

	F = meshCells;

	/// ----------------------------------------------------------------------

	 // number of vertices on the largest side
	const int s = 50;
	const Eigen::RowVector3d Vmin = V.colwise().minCoeff();
	const Eigen::RowVector3d Vmax = V.colwise().maxCoeff();
	const double h = (Vmax - Vmin).maxCoeff() / (double)s;
	const Eigen::RowVector3i res = (s * ((Vmax - Vmin) / (Vmax - Vmin).maxCoeff())).cast<int>();

	// create grid
	Eigen::MatrixXd GV(res(0) * res(1) * res(2), 3);
	for (int zi = 0; zi < res(2); zi++)
	{
		const auto lerp = [&](const int di, const int d)->double
		{return Vmin(d) + (double)di / (double)(res(d) - 1) * (Vmax(d) - Vmin(d)); };
		const double z = lerp(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double y = lerp(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double x = lerp(xi, 0);
				GV.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(x, y, z);
			}
		}
	}
	// compute values
	Eigen::VectorXd S, B;
	{
		Eigen::VectorXi I;
		Eigen::MatrixXd C, N;
		igl::signed_distance(GV, V, F, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, S, I, C, N);
		// Convert distances to binary inside-outside data --> aliasing artifacts
		B = S;
		std::for_each(B.data(), B.data() + B.size(), [](double& b) {b = (b > 0 ? 1 : (b < 0 ? -1 : 0)); });
	}
	Eigen::MatrixXd SV, BV;
	Eigen::MatrixXi SF, BF;
	igl::copyleft::marching_cubes(S, GV, res(0), res(1), res(2), SV, SF);
	igl::copyleft::marching_cubes(B, GV, res(0), res(1), res(2), BV, BF);


// '1'  Show original mesh.
// '2'  Show marching cubes contour of signed distance.
// '3'  Show marching cubes contour of indicator function.


	// Send new positions
	vtkNew<vtkPolyData> meshOutput;
	lb3d.ConvertEigenToVTKMesh(SV, SF, meshOutput);
	outputData->DeepCopy(meshOutput);

	switch (key)
	{
	default:
		return 0;
	case '1':
		// Send new positions
		lb3d.ConvertEigenToVTKMesh(V, F, meshOutput);
		outputData->DeepCopy(meshOutput);
		break;
	case '2':
		// Send new positions
		lb3d.ConvertEigenToVTKMesh(SV, SF, meshOutput);
		outputData->DeepCopy(meshOutput);
		break;
	case '3':
		// Send new positions
		lb3d.ConvertEigenToVTKMesh(BV, BF, meshOutput);
		outputData->DeepCopy(meshOutput);
		break;
	}


}
