#pragma once

#include "Utils/STL_Macros.h"
#include "Utils/GMM_Macros.h"

class SpaceDeformer2D : public MPxDeformerNode
{
public:
	SpaceDeformer2D();
	virtual ~SpaceDeformer2D();

	static void* creator();
	static MStatus initialize();

	virtual MStatus deform(MDataBlock& block, MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex);//deformed position of each internal point is computed on the fly
	MStatus updateCage(MFnMesh& cageMeshFn);
	MStatus doSetup(MItGeometry& iter, MFnMesh& cageMeshFn);//pre-proccessing , calculati the bary coordinates


public:
	const static MTypeId mTypeId;
	const static MString mTypeName;

protected:
	static MObject mCageAttr;
	static MObject mCoordinateTypeAttr;

protected:
	bool mIsFirstTime;

	GMMDenseComplexColMatrix mCageVertices; //this matrix is actually a column vector. dimensions are: n x 1
	GMMDenseComplexColMatrix mCauchyCoordinates; //this matrix is actually a column vector. dimensions are: m x n
	GMMDenseComplexColMatrix mInternalPoints; //this matrix is actually a column vector. dimensions are: m x 1

	GMMDenseComplexColMatrix mCauchyCoordsForInterpolation; //dimensions are: n x n
	GMMDenseComplexColMatrix mInterpolatedNewVertices_f; //dimensions are: n x 1
private:
	void SpaceDeformer2D::matlabCalcNewVerticesForInterpolation();

};
