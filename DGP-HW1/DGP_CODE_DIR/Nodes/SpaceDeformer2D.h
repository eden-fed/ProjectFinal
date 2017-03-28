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
	static MObject mEpsilonAttr;

	bool mIsFirstTime;
	bool mNeedToCalcNewLine;
	float env;
	long coordinateType;
	int mNumOfSegmentsA;
	int mNLarge;
	int mNumOfSegmentsAOld;
	int mNLargeOld;
	int mCurrentNumOfSegmentsA;
	int mCurrentNLarge;
	double k;
	double SigmaA;
	double sigmaB;
	double lambda;
	int iterationsNum;
	double epsilon;
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

	std::vector<double> mXvaluesOfIntersections;
	std::vector<double> mYvaluesOfIntersections;

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

	GMMDenseComplexColMatrix mPinvOfIncCageVertexCoords; //dimensions are: nLarge x a
	GMMSparseComplexRowMatrix mLMatrixForLipmansMethod; //dimensions are: 2*nLarge x 2*nLarge
	GMMSparseComplexRowMatrix mUMatrixForLipmansMethod; //dimensions are: 2*nLarge x 2*nLarge
	GMMDenseColMatrix mPRowPermOfLUForLipmansMethod; //dimensions are: 1 x 2*nLarge
	GMMDenseComplexColMatrix mTtrasposeForLipmansMethod;//dimensions are: 2*nLarge x a

	GMMDenseComplexColMatrix mAOfLineSegmentInLvAccelerated; //dimensions are: 5 x 1 **unused**
	GMMDenseComplexColMatrix mBOfLineSegmentInLvAccelerated; //dimensions are: 5 x 1 **unused**
	/*GMMDenseColMatrix mSlopeOfLineApproxForCurveInLv;
	GMMDenseColMatrix mXcoordOfIntersectionPointForCurveLv;*/

	double mSlopeOfLineApproxForCurveInLv;
	double mXcoordOfIntersectionPointForCurveLv;
	double mYcoordOfIntersectionPointForCurveLv;

	GMMDenseComplexColMatrix mUserCageVerticesNos_sizeA;
	GMMDenseComplexColMatrix mUserCageVerticesNos_sizeNLarge;
	
	GMMDenseComplexColMatrix mCompCageVerticesNos_sizeA;
	GMMDenseComplexColMatrix mCompCageVerticesNos;

	GMMDenseColMatrix mEndIndicesForIntegral; //dimensions are: 1 X m-1
	GMMDenseColMatrix mStartIndicesForIntegral; //dimensions are: 1 X m-1
	GMMDenseComplexColMatrix mEdgeVectorsForIntegral; //dimensions are: m X 1

	GMMDenseComplexColMatrix mCauchyCoordsOfz0; //dimensions are: 1 x nLarge

private:
	void matlabCalcNewVerticesForInterpolation();
	void matlabCalcNewVerticesForP2P();
	void matlabCalcLforHprojection();
	void matlabCalcLforLvprojection_curve();
	void matlabCalcLforLvprojection();
	void matlabCalcLforLvprojectionConformalAccel();
	void matlabCalcLforLvprojectionConformal();
	void matlabCalcLforLvprojectionAccel();
	void matlabCalcLforLvprojectionLipman();
	void matlabCalcLforLvprojectionDykstra();
	std::string RelativeToFullPath(char* relPath);
	MStatus runTimeDoSetup();
	int findClosestInternalPointsToZ0();
	MStatus preprocessingIntegral(MFnMesh& inputMesh, MObject InputGeom);
	MStatus calcSegments();
	void IncreaseVertecies(Complex* OriginalCompCageVertecies, int OrigCageSize, Complex** IncreasedCompCageVertecies, int& numOfIncreasedCageVertecies);
	void IncreaseVertecies(MPointArray& OriginalCageVertecies, MPointArray& IncreasedCageVertecies, int numOfIncreasedCageVertecies, bool countNumOfVerticesInEdges);
	MStatus showIncVertecies(MPointArray& IncreasedCageVertecies);
	MStatus findLineApproximationForCurve();

	void calcLvprojectionLGcpu();
	void evalFzAndFzBar(GMMDenseComplexColMatrix& CageVerticesNos, GMMDenseComplexColMatrix& userCageVerticesNos, GMMDenseColMatrix& mNumOfVerticesInEdges, int numOfIncreasedCageVertecies, int numOfCageVertices, GMMDenseComplexColMatrix& fz, GMMDenseComplexColMatrix& fzBar);
	void IncreaseVerteciesAfterMap(GMMDenseComplexColMatrix& OriginalCageVertecies, GMMDenseComplexColMatrix& IncreasedCageVertecies, int numOfIncreasedCageVertecies, GMMDenseColMatrix& mNumOfVerticesInEdges);
	void logarithmExtraction(GMMDenseComplexColMatrix& cageVerticesNos_sizeA, GMMDenseComplexColMatrix& fz, GMMDenseComplexColMatrix& cageVerteciesAfterMapSizeA, int a, GMMDenseComplexColMatrix& log_fz);
	void find_nu_f(GMMDenseComplexColMatrix& fz, GMMDenseComplexColMatrix& fzBar, int a, GMMDenseComplexColMatrix& nu_f);
	int doLocalGlobalIterations(GMMDenseComplexColMatrix& log_fz, GMMDenseComplexColMatrix& nu_f, GMMDenseComplexColMatrix& l, GMMDenseComplexColMatrix& nu, void (SpaceDeformer2D::*projectionFunction)(double&, double&));
	bool localStep(GMMDenseComplexColMatrix& log_fz, GMMDenseComplexColMatrix& nu_f, void (SpaceDeformer2D::*projectionFunction)(double&, double&));
	bool checkIfInsidePolygon(double x, double y);
	void projectPointToPolygonMinSeg(double& x, double& y);
	void projectPointToPolygonWithK(double& x, double& y);
	void projectPointToPolygonNoK(double& x, double& y);
	void findPHI(GMMDenseComplexColMatrix& PHI, GMMDenseComplexColMatrix& PHItag, GMMDenseComplexColMatrix& LonInternalPoints);
	void findPSI(GMMDenseComplexColMatrix& PSI, GMMDenseComplexColMatrix& NUonInternalPoints, GMMDenseComplexColMatrix& PHItag);
	void calcIntegralUsingSpaningTree(GMMDenseComplexColMatrix& PHI, GMMDenseComplexColMatrix& PHItag,Complex phi_Z0);

};
