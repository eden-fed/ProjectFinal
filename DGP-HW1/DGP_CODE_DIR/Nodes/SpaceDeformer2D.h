#pragma once

#include "Utils/STL_Macros.h"
#include "Utils/GMM_Macros.h"
#include "Utils/FastMatrix.h"
#include "Utils/MatrixCPU.h"
#include "Utils/MatrixGPU.h"

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
	MAYA_float3 mZ0NotOnMesh;
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
	//GMMSparseComplexRowMatrix mLMatrixForLipmansMethod; //dimensions are: 2*nLarge x 2*nLarge
	//GMMSparseComplexRowMatrix mUMatrixForLipmansMethod; //dimensions are: 2*nLarge x 2*nLarge
	//GMMDenseComplexColMatrix mTtrasposeForLipmansMethod;//dimensions are: 2*nLarge x a
	//GMMDenseComplexColMatrix mTForLipmansMethod;//dimensions are: 2*nLarge x a
	GMMDenseComplexColMatrix mInvMtransCForLipmansMethod;//dimensions are: 2*nLarge x a

	GMMDenseComplexColMatrix mAOfLineSegmentInLvAccelerated; //dimensions are: 5 x 1 **unused**
	GMMDenseComplexColMatrix mBOfLineSegmentInLvAccelerated; //dimensions are: 5 x 1 **unused**
	/*GMMDenseColMatrix mSlopeOfLineApproxForCurveInLv;
	GMMDenseColMatrix mXcoordOfIntersectionPointForCurveLv;*/

	double mSlopeOfLineApproxForCurveInLv;
	double mXcoordOfIntersectionPointForCurveLv;
	double mYcoordOfIntersectionPointForCurveLv;

	ComplexDoubleCPUMatrix mUserCageVerticesNos_sizeA;
	ComplexDoubleCPUMatrix mUserCageVerticesNos_sizeNLarge;
	
	GMMDenseComplexColMatrix mCompCageVerticesNos_sizeA;
	GMMDenseComplexColMatrix mCompCageVerticesNos;

	GMMDenseColMatrix mEndIndicesForIntegral; //dimensions are: 1 X m-1
	GMMDenseColMatrix mStartIndicesForIntegral; //dimensions are: 1 X m-1
	GMMDenseComplexColMatrix mEdgeVectorsForIntegral; //dimensions are: m X 1

	GMMDenseComplexColMatrix mCauchyCoordsOfz0; //dimensions are: 1 x nLarge

	//***arrays for L_nu projection in c++***
	ComplexDoubleCPUMatrix mfz;
	ComplexDoubleCPUMatrix mfzBar;
	ComplexDoubleCPUMatrix mLog_fz;
	ComplexDoubleCPUMatrix mNu_f;
	ComplexDoubleGPUMatrix mLog_fz_gpu;
	ComplexDoubleGPUMatrix mNu_f_gpu;
	ComplexDoubleGPUMatrix mL_gpu;
	ComplexDoubleGPUMatrix mNu_gpu;

	std::vector<int> mNumOfVerticesInEdgesSizeA_stdVec;//dimensions are: numOfEdges x 1 - sum of vertices=A
	std::vector<int> mNumOfVerticesInEdgesSizeNlarge_stdVec;//dimensions are: numOfEdges x 1 - sum of vertices=nLarge

	ComplexDoubleGPUMatrix mPinvOfIncCageVertexCoords_gpuMat;
	ComplexDoubleGPUMatrix mIncCageVertexCoords_gpuMat;

	ComplexDoubleCPUMatrix mLonInternalPoints;
	ComplexDoubleCPUMatrix mNUonInternalPoints;
	ComplexDoubleGPUMatrix mLonInternalPoints_gpu;
	ComplexDoubleGPUMatrix mNUonInternalPoints_gpu;

	ComplexDoubleGPUMatrix mCauchyCoordinatesIncForP2P_gpuMat;
	ComplexDoubleCPUMatrix mPHI;
	ComplexDoubleCPUMatrix mPHItag;
	ComplexDoubleCPUMatrix mPSI;

	ComplexDoubleGPUMatrix mCauchyCoordsOfz0_gpuMat;

	DoubleGPUMatrix mXvaluesOfIntersections_gpu;//test
	DoubleGPUMatrix mYvaluesOfIntersections_gpu;//test


	ComplexDoubleGPUMatrix mX_gpu;
	ComplexDoubleGPUMatrix mX_local_gpu;
	ComplexDoubleGPUMatrix mn_0forStopCondition_gpu;
	ComplexDoubleGPUMatrix mn_0_nu_forLipmansMethod_gpu;
	ComplexDoubleGPUMatrix mn_0_l_forLipmansMethod_gpu;
	ComplexDoubleGPUMatrix mTempCalc_nu_forLipmansMethod_gpu;
	ComplexDoubleGPUMatrix mTempCalc_l_forLipmansMethod_gpu; 
	//ComplexDoubleGPUMatrix mc_0forLipmansMethod_gpu;
	//ComplexDoubleGPUMatrix my_cforLipmansMethod_gpu;
	ComplexDoubleGPUMatrix my_eta_l_forLipmansMethod_gpu;
	ComplexDoubleGPUMatrix my_eta_nu_forLipmansMethod_gpu;
	//ComplexDoubleGPUMatrix mTtrasposeForLipmansMethod_gpu;
	//ComplexDoubleGPUMatrix mTForLipmansMethod_gpu;
	ComplexDoubleGPUMatrix mInvMtransCForLipmansMethod_gpu;
	//ComplexDoubleGPUMatrix mLnu_forLipman_gpu;

	//****

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
	//MStatus showIncVertecies(MPointArray& IncreasedCageVertecies);
	MStatus findLineApproximationForCurve();

	int calcLvprojectionLGgpu();
	int calcLvprojectionHPgpu();
	void evalFzAndFzBar(GMMDenseComplexColMatrix& CageVerticesNos, GMMDenseComplexColMatrix& userCageVerticesNos, std::vector<int>& mNumOfVerticesInEdges, int numOfIncreasedCageVertecies, int numOfCageVertices, ComplexDoubleCPUMatrix& fz, ComplexDoubleCPUMatrix& fzBar);
	void IncreaseVerteciesAfterMap(GMMDenseComplexColMatrix& OriginalCageVertecies, ComplexDoubleCPUMatrix& IncreasedCageVertecies, int numOfIncreasedCageVertecies, std::vector<int>& mNumOfVerticesInEdges);
	void logarithmExtraction(GMMDenseComplexColMatrix& cageVerticesNos_sizeA, ComplexDoubleCPUMatrix& fz, ComplexDoubleCPUMatrix& cageVerteciesAfterMapSizeA, int a, ComplexDoubleCPUMatrix& log_fz);
	void find_nu_f(ComplexDoubleCPUMatrix& fz, ComplexDoubleCPUMatrix& fzBar, int a, ComplexDoubleCPUMatrix& nu_f);
	//int doLocalGlobalIterations(void (SpaceDeformer2D::*projectionFunction)(double&, double&));
	//int doLocalGlobalIterations(bool (SpaceDeformer2D::*localStepFunction)(ComplexDoubleGPUMatrix&, ComplexDoubleGPUMatrix&));
	int doLocalGlobalIterations();
	int doHyperPlaneIterations();
	bool localStep(ComplexDoubleCPUMatrix& log_fz, ComplexDoubleCPUMatrix& nu_f, void (SpaceDeformer2D::*projectionFunction)(double&, double&));
	bool checkIfInsidePolygon(double x, double y);
	void projectPointToPolygonMinSeg(double& x, double& y);
	void projectPointToPolygonWithK(double& x, double& y);
	void projectPointToPolygonNoK(double& x, double& y);
	void findPHI(ComplexDoubleCPUMatrix& PHI, ComplexDoubleCPUMatrix& PHItag, ComplexDoubleCPUMatrix& LonInternalPoints);
	void findPSI(ComplexDoubleCPUMatrix& PSI, ComplexDoubleCPUMatrix& NUonInternalPoints, ComplexDoubleCPUMatrix& PHItag);
	void calcIntegralUsingSpaningTree(ComplexDoubleCPUMatrix& f, ComplexDoubleCPUMatrix& f_tag, Complex phi_Z0);

	bool localStep_noK_gpu(ComplexDoubleGPUMatrix& log_fz, ComplexDoubleGPUMatrix& nu_f);
	bool localStep_withK_gpu(ComplexDoubleGPUMatrix& log_fz, ComplexDoubleGPUMatrix& nu_f);
	bool localStep_noK_HP_gpu(ComplexDoubleGPUMatrix& x_vec);
	bool localStep_withK_HP_gpu(ComplexDoubleGPUMatrix& x_vec);
	bool localStep_minSeg_gpu(ComplexDoubleGPUMatrix& log_fz, ComplexDoubleGPUMatrix& nu_f);
	bool localStep_minSeg_HP_gpu(ComplexDoubleGPUMatrix& x_vec);
};
