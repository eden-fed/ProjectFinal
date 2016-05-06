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
	MStatus updateControlPoints(MFnMesh& cageMeshFn);

	MStatus doSetup(MItGeometry& iter, MFnMesh& cageMeshFn);//pre-proccessing , calculati the bary coordinates


public:
	const static MTypeId mTypeId;
	const static MString mTypeName;

protected:
	static MObject mCageAttr;
	static MObject mCageP2pAttr;
	static MObject mCoordinateTypeAttr;
	static MObject mNumOfSegmentsAttr;


protected:
	bool mIsFirstTime;
	int mNumOfSegments;

	GMMDenseComplexColMatrix mUserCageVertices; //this matrix is actually a column vector. dimensions are: n x 1
	GMMDenseComplexColMatrix mCauchyCoordinates; //dimensions are: m x n
	GMMDenseComplexColMatrix mInternalPoints; //this matrix is actually a column vector. dimensions are: m x 1

	GMMDenseComplexColMatrix mCauchyCoordsOfOriginalCageVertices; //dimensions are: n x n
	GMMDenseComplexColMatrix mInterpolationGenCage_f; //dimensions are: n x 1

	GMMDenseComplexColMatrix mUserP2P; //dimentions are k x 1
	GMMDenseComplexColMatrix mCauchyCoordsOfOriginalP2P; //dimensions are: k x nLarge
	GMMDenseComplexColMatrix mP2PGenCageVertices_f; //dimensions are: nLarge x 1

	//************
	MPointArray mIncreasedVerteciesInitial_n; //dimensions are: l x 1
	GMMDenseComplexColMatrix mCauchyCoordinatesIncForP2P; //dimensions are: m x nLarge
	//*************

	MPointArray mIncreasedVertecies_a; //dimensions are: l x 1
	GMMDenseComplexColMatrix mSecondDerOfIncCageVertexCoords; //dimensions are: l x nLarge

private:
	void matlabCalcNewVerticesForInterpolation();
	void matlabCalcNewVerticesForP2P();
	std::string RelativeToFullPath(char* relPath);

};
