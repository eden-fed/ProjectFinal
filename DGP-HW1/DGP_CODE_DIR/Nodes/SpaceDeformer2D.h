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
	MStatus getData(MDataBlock& block, MObject& cageMesh, MObject& p2pMesh);
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
	static MObject mNlargeAttr;
	static MObject mkAttr;
	static MObject mSigmaaAttr;
	static MObject msigmabAttr;
	static MObject mZ0Attr;
	static MObject mlambdaAttr;
	static MObject mIterAttr;

	bool mIsFirstTime;
	float env;
	long coordinateType;
	int mNumOfSegmentsA;
	int mNLarge;
	int mNumOfSegmentsAOld;
	int mNLargeOld;
	double k;
	double SigmaA;
	double sigmaB;
	double lambda;
	int iterationsNum;
	float3 mZ0NotOnMesh;
	int mZ0index;
	Complex mZ0onMesh;

	int mNumOfInternalPoints;
	int mNumOfCageVerticies;
	int mNumOfControlPoints;
	MPointArray mCartCageVerticesNos; //cartesian coordinates
	MPointArray mCartCageVerticesNos_sizeA; //cartesian coordinates
	Complex* mCompCageVerticesWos;
	MPointArray mInternalPoints_MPoint; //mpoint array for runtime dosetup
	MPointArray mInitialcontrolPoints;

	MObject mcageMesh;


	GMMDenseComplexColMatrix mUserCageVerticesNos; //this matrix is actually a column vector. dimensions are: n x 1
	GMMDenseComplexColMatrix mCauchyCoordinates; //dimensions are: m x n
	GMMDenseComplexColMatrix mInternalPoints; //this matrix is actually a column vector. dimensions are: m x 1

	GMMDenseComplexColMatrix mCauchyCoordsOfOriginalCageVertices; //dimensions are: n x n
	GMMDenseComplexColMatrix mInterpolationGenCage_f; //dimensions are: n x 1

	GMMDenseComplexColMatrix mUserP2P; //dimentions are k x 1
	GMMDenseComplexColMatrix mCauchyCoordsOfOriginalP2P; //dimensions are: k x nLarge
	GMMDenseComplexColMatrix mP2PGenCageVertices_f; //dimensions are: nLarge x 1

	//************
	GMMDenseComplexColMatrix mCauchyCoordinatesIncForP2P; //dimensions are: m x nLarge
	//*************
	GMMDenseComplexColMatrix mIncCageVertexCoords; //dimensions are: a x nLarge
	GMMDenseComplexColMatrix mSecondDerOfIncCageVertexCoords; //dimensions are: a x nLarge
	GMMDenseComplexColMatrix mFirstDerOfIncCageVertexCoords; //dimensions are: a x nLarge
	//******************************

	GMMDenseColMatrix mNumOfVerticesInEdgesSizeA;//dimensions are: numOfEdges x 1 - sum of vertices=A
	GMMDenseColMatrix mNumOfVerticesInEdgesSizeNlarge;//dimensions are: numOfEdges x 1 - sum of vertices=nLarge

private:
	void matlabCalcNewVerticesForInterpolation();
	void matlabCalcNewVerticesForP2P();
	void matlabCalcLforHprojection();
	void matlabCalcLforLvprojection();
	void matlabCalcLforLvprojectionAccel();
	std::string RelativeToFullPath(char* relPath);
	MStatus runTimeDoSetup();
	int findClosestInternalPointsToZ0();
	MStatus preprocessingIntegral(MFnMesh& inputMesh, MObject InputGeom);
	void IncreaseVertecies(Complex* OriginalCompCageVertecies, int OrigCageSize, Complex** IncreasedCompCageVertecies, int& numOfIncreasedCageVertecies);
	void IncreaseVertecies(MPointArray& OriginalCageVertecies, MPointArray& IncreasedCageVertecies, int numOfIncreasedCageVertecies, bool countNumOfVerticesInEdges);
	MStatus showIncVertecies(MPointArray& IncreasedCageVertecies);

};
