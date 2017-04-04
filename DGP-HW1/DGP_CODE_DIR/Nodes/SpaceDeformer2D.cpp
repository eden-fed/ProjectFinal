#include "stdafx.h"

#include "SpaceDeformer2D.h"
#include "Utils/STL_Macros.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "Utils/Timer.h"
#include "GPULocalStep.h"

#define EPS 0.0001
#define CVX_INTERPOLATION 

#define IS_NUMERIC_ZERO(a) (abs(a)<0.000001?1:0)

#define IN
#define OUT

//#define DEBUG
//#define DEBUG2
//#define COMPARE_TOTAL_TIMINGS
#define COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
//#define DEBUG_FINAL_RESULT


const MTypeId SpaceDeformer2D::mTypeId(0x6723c);
const MString SpaceDeformer2D::mTypeName("SpaceDeformer2D");

MObject SpaceDeformer2D::mCageAttr;
MObject SpaceDeformer2D::mCageP2pAttr;
MObject SpaceDeformer2D::mCoordinateTypeAttr;
MObject SpaceDeformer2D::mNumOfSegmentsAttr;
MObject SpaceDeformer2D::mNlargeAttr;
MObject SpaceDeformer2D::mkAttr;
MObject SpaceDeformer2D::mSigmaaAttr;
MObject SpaceDeformer2D::msigmabAttr;
MObject SpaceDeformer2D::mZ0Attr;
MObject SpaceDeformer2D::mlambdaAttr;
MObject SpaceDeformer2D::mIterAttr;
MObject SpaceDeformer2D::mEpsilonAttr;



SpaceDeformer2D::SpaceDeformer2D() : mIsFirstTime(true), mNeedToCalcNewLine(true)
{
	mCompCageVerticesWos = NULL;
	mNumOfSegmentsAOld = mNumOfSegmentsA = 0;
	mNLargeOld = mNLarge = 0;
}

SpaceDeformer2D::~SpaceDeformer2D()
{
	if (mCompCageVerticesWos!=NULL)
		delete[] mCompCageVerticesWos;
}

void* SpaceDeformer2D::creator()
{
	return new SpaceDeformer2D();
}

MStatus SpaceDeformer2D::initialize()
{
	MStatus stat;

	MFnTypedAttribute cageAttr;
	mCageAttr = cageAttr.create("cage", "cage", MFnData::kMesh, MObject::kNullObj, &stat);
	CHECK_MSTATUS(addAttribute(mCageAttr));
	CHECK_MSTATUS(attributeAffects(mCageAttr, outputGeom));

	MFnTypedAttribute p2pCage;
	mCageP2pAttr = p2pCage.create("p2pcage", "p2pcage", MFnData::kMesh, MObject::kNullObj, &stat);
	CHECK_MSTATUS(addAttribute(mCageP2pAttr));
	CHECK_MSTATUS(attributeAffects(mCageP2pAttr, outputGeom));

	MFnEnumAttribute coordinateTypeAttr;
	mCoordinateTypeAttr = coordinateTypeAttr.create("coordinateType", "coordinateType", 0, &stat);
	CHECK_MSTATUS(coordinateTypeAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mCoordinateTypeAttr));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Cauchy", 0));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Cauchy Interpolation", 1));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Point to point", 2));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Projection H Space", 3));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Projection Lv Space Conformal", 4));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Projection Lv Space Conformal Accel", 5));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Projection Lv Space, cvx curve", 6));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Projection Lv Space, cvx", 7));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Projection Lv Space, L/G", 8));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Projection Lv Space, Lipman", 9));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Projection Lv Space, Dykstra", 10));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Projection Lv Space, L/G gpu", 11));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Projection Lv Space, Lipman gpu", 12));


	CHECK_MSTATUS(attributeAffects(mCoordinateTypeAttr, outputGeom));

	MFnNumericAttribute numOfSegmentsAttr;
	mNumOfSegmentsAttr = numOfSegmentsAttr.create("numOfSegmentsA", "numOfSegmentsA", MFnNumericData::kInt, 1000, &stat);
	CHECK_MSTATUS(numOfSegmentsAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mNumOfSegmentsAttr));
	CHECK_MSTATUS(attributeAffects(mNumOfSegmentsAttr, outputGeom));

	MFnNumericAttribute nLargeAttr;
	mNlargeAttr = nLargeAttr.create("numOfSegmentsN", "nLarge", MFnNumericData::kInt, 150, &stat);
	CHECK_MSTATUS(nLargeAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mNlargeAttr));
	CHECK_MSTATUS(attributeAffects(mNlargeAttr, outputGeom));

	MFnNumericAttribute kAttr;
	mkAttr = kAttr.create("k", "k", MFnNumericData::kDouble, 0.4, &stat);
	CHECK_MSTATUS(kAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mkAttr));
	CHECK_MSTATUS(attributeAffects(mkAttr, outputGeom));

	MFnNumericAttribute sigmaaAttr;
	mSigmaaAttr = sigmaaAttr.create("Sigma(a)", "Sigma(a)", MFnNumericData::kDouble, 2, &stat);
	CHECK_MSTATUS(sigmaaAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mSigmaaAttr));
	CHECK_MSTATUS(attributeAffects(mSigmaaAttr, outputGeom));

	MFnNumericAttribute sigmabAttr;
	msigmabAttr = sigmabAttr.create("sigma(b)", "sigma(b)", MFnNumericData::kDouble, 0.5, &stat);
	CHECK_MSTATUS(sigmabAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(msigmabAttr));
	CHECK_MSTATUS(attributeAffects(msigmabAttr, outputGeom));

	MFnNumericAttribute z0Attr;
	mZ0Attr = z0Attr.create("z0", "z0", MFnNumericData::k3Float, 0.5, &stat);
	CHECK_MSTATUS(z0Attr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mZ0Attr));
	CHECK_MSTATUS(attributeAffects(mZ0Attr, outputGeom));

	MFnNumericAttribute lambdaAttr;
	mlambdaAttr = lambdaAttr.create("lambda", "lambda", MFnNumericData::kDouble, 1, &stat);
	CHECK_MSTATUS(lambdaAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mlambdaAttr));
	CHECK_MSTATUS(attributeAffects(mlambdaAttr, outputGeom));

	MFnNumericAttribute iterAttr;
	mIterAttr = iterAttr.create("maxIterations", "maxIterations", MFnNumericData::kInt, 2000, &stat);
	CHECK_MSTATUS(iterAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mIterAttr));
	CHECK_MSTATUS(attributeAffects(mIterAttr, outputGeom));

	MFnNumericAttribute epsAttr;
	mEpsilonAttr = epsAttr.create("epsilon", "epsilon", MFnNumericData::kDouble, 1e-4, &stat);
	CHECK_MSTATUS(epsAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mEpsilonAttr));
	CHECK_MSTATUS(attributeAffects(mEpsilonAttr, outputGeom));

	return MStatus::kSuccess;
}

void SpaceDeformer2D::matlabCalcNewVerticesForInterpolation() {
	MatlabGMMDataExchange::SetEngineDenseMatrix("C", mCauchyCoordsOfOriginalCageVertices);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("q", mUserCageVerticesNos);//send the matrix to matlab

//	int res = MatlabInterface::GetEngine().LoadAndRunScript("C:/Users/eden/Documents/MySWProjects/ProjectFinal/DGP-HW1/DGP_CODE_DIR/interpolatedCauchy.m");
	
	int res = MatlabInterface::GetEngine().LoadAndRunScript(RelativeToFullPath("\\matlab scripts\\interpolatedCauchy.m").c_str());
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'interpolatedCauchy.m' failed with error code " << res << std::endl;
	}
	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mInterpolationGenCage_f);//get the incersed matrix from matlab
}

void SpaceDeformer2D::matlabCalcNewVerticesForP2P() {
	MatlabGMMDataExchange::SetEngineDenseMatrix("C", mCauchyCoordsOfOriginalP2P);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("q", mUserP2P);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("D", mSecondDerOfIncCageVertexCoords);//send the matrix to matlab
	
	int res = MatlabInterface::GetEngine().LoadAndRunScript(RelativeToFullPath("\\matlab scripts\\P2P.m").c_str());
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'P2P.m' failed with error code " << res << std::endl;
	}
	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mP2PGenCageVertices_f);//get the incersed matrix from matlab
}
GMMDenseColMatrix doubleToGmmMat(const double value)
{
	GMMDenseColMatrix retVal(1, 1);
	retVal(0, 0) = value;
	return retVal;
}
GMMDenseComplexColMatrix compToGmmMat(const Complex value)
{
	GMMDenseComplexColMatrix retVal(1, 1);
	retVal(0, 0) = value;
	return retVal;
}
GMMDenseComplexColMatrix compPointArrayToGmmMat(const MPointArray array)
{
	int numV = array.length();
	GMMDenseComplexColMatrix retVal(numV, 1);
	for (int i = 0; i < numV; i++)
	{
		MPoint p = array[i];
		Complex c(p[0], p[1]);
		retVal(i, 0) = c;
	}
	return retVal;
}

void gmmToStdIntVector(const GMMDenseColMatrix& gmmVec, std::vector<int>& vec){
	vec.clear();
	for (int i = 0; i < gmmVec.size(); i++)
		vec.push_back(gmmVec(i, 0));
}

void gmmToCpuMatrix(const GMMDenseComplexColMatrix& gmmMat, ComplexDoubleCPUMatrix& cpuMat){
	for (int i = 0; i < gmmMat.nrows(); i++)
		for (int j = 0; j < gmmMat.ncols(); j++)
			cpuMat(i, j) = gmmMat(i, j);
}

void SpaceDeformer2D::matlabCalcLforHprojection()
{
	//MatlabInterface::GetEngine().Eval("clearvars -except edgeVectors startIndices endIndices");

	MatlabGMMDataExchange::SetEngineDenseMatrix("NumOfVerticesInEdges", mNumOfVerticesInEdgesSizeA);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("Ctag", mFirstDerOfIncCageVertexCoords);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("C_sizeM", mCauchyCoordinatesIncForP2P);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("C_sizeA", mIncCageVertexCoords);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("k", doubleToGmmMat(this->k));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("SIGMA", doubleToGmmMat(this->SigmaA));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("sigma", doubleToGmmMat(this->sigmaB));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("Z0index", doubleToGmmMat((this->mZ0index)+1));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("Z0", compToGmmMat(mZ0onMesh));
	MatlabGMMDataExchange::SetEngineDenseMatrix("lambda", doubleToGmmMat(this->lambda));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMap", mUserCageVerticesNos);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesB4Map", compPointArrayToGmmMat(mCartCageVerticesNos));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesB4Map_sizeA", mCompCageVerticesNos_sizeA);//send the matrix to matlab

	std::string res = MatlabInterface::GetEngine().LoadAndRunScriptToString(RelativeToFullPath("\\matlab scripts\\projectToH.m").c_str());
	std::cerr << res;

	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mInternalPoints);//get the map from matlab

	cout.flush();
}
void SpaceDeformer2D::matlabCalcLforLvprojection_curve()
{

/*	IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeA, mCurrentNumOfSegmentsA, mNumOfVerticesInEdgesSizeA);
	IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeNLarge, mCurrentNLarge, mNumOfVerticesInEdgesSizeNlarge);
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeA", mUserCageVerticesNos_sizeA);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeNLarge", mUserCageVerticesNos_sizeNLarge);//send the matrix to matlab*/

	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMap", mUserCageVerticesNos);//send the matrix to matlab

	std::string res = MatlabInterface::GetEngine().EvalToString("projectToLv_curve;");
	std::cerr << res;


	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mInternalPoints);//get the map from matlab
	cout.flush();

}


void SpaceDeformer2D::matlabCalcLforLvprojection()
{
	/*IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeA, mCurrentNumOfSegmentsA, mNumOfVerticesInEdgesSizeA);
	IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeNLarge, mCurrentNLarge, mNumOfVerticesInEdgesSizeNlarge);
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeA", mUserCageVerticesNos_sizeA);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeNLarge", mUserCageVerticesNos_sizeNLarge);//send the matrix to matlab*/

	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMap", mUserCageVerticesNos);//send the matrix to matlab

	std::string res = MatlabInterface::GetEngine().EvalToString("projectToLv;");
	std::cerr << res;

	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mInternalPoints);//get the map from matlab
	cout.flush();
}
void SpaceDeformer2D::matlabCalcLforLvprojectionConformalAccel()
{

	/*IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeA, mCurrentNumOfSegmentsA, mNumOfVerticesInEdgesSizeA);
	IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeNLarge, mCurrentNLarge, mNumOfVerticesInEdgesSizeNlarge);
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeA", mUserCageVerticesNos_sizeA);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeNLarge", mUserCageVerticesNos_sizeNLarge);//send the matrix to matlab*/

	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMap", mUserCageVerticesNos);//send the matrix to matlab

	std::string res = MatlabInterface::GetEngine().EvalToString("projectToLv_accel_conformal;");
	std::cerr << res;

	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mInternalPoints);//get the map from matlab
	cout.flush();

}
void SpaceDeformer2D::matlabCalcLforLvprojectionConformal()
{

	/*IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeA, mCurrentNumOfSegmentsA, mNumOfVerticesInEdgesSizeA);
	IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeNLarge, mCurrentNLarge, mNumOfVerticesInEdgesSizeNlarge);
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeA", mUserCageVerticesNos_sizeA);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeNLarge", mUserCageVerticesNos_sizeNLarge);//send the matrix to matlab*/

	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMap", mUserCageVerticesNos);//send the matrix to matlab

	std::string res = MatlabInterface::GetEngine().EvalToString("projectToLv_comformal;");
	std::cerr << res;

	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mInternalPoints);//get the map from matlab
	cout.flush();

}
void SpaceDeformer2D::matlabCalcLforLvprojectionAccel()
{

	/*IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeA, mCurrentNumOfSegmentsA, mNumOfVerticesInEdgesSizeA);
	IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeNLarge, mCurrentNLarge, mNumOfVerticesInEdgesSizeNlarge);

	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeA", mUserCageVerticesNos_sizeA);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeNLarge", mUserCageVerticesNos_sizeNLarge);//send the matrix to matlab*/

	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMap", mUserCageVerticesNos);//send the matrix to matlab

	std::string res = MatlabInterface::GetEngine().EvalToString("projectToLv_accel;");
	std::cerr << res;

	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mInternalPoints);//get the map from matlab


	cout.flush();
}
void SpaceDeformer2D::matlabCalcLforLvprojectionDykstra()
{

	/*IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeA, mCurrentNumOfSegmentsA, mNumOfVerticesInEdgesSizeA);
	IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeNLarge, mCurrentNLarge, mNumOfVerticesInEdgesSizeNlarge);
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeA", mUserCageVerticesNos_sizeA);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeNLarge", mUserCageVerticesNos_sizeNLarge);//send the matrix to matlab*/

	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMap", mUserCageVerticesNos);//send the matrix to matlab

	std::string res = MatlabInterface::GetEngine().EvalToString("projectToLv_dykstra;");
	std::cerr << res;

	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mInternalPoints);//get the map from matlab
	cout.flush();
}
void SpaceDeformer2D::matlabCalcLforLvprojectionLipman()
{

	/*IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeA, mCurrentNumOfSegmentsA, mNumOfVerticesInEdgesSizeA);
	IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeNLarge, mCurrentNLarge, mNumOfVerticesInEdgesSizeNlarge);
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeA", mUserCageVerticesNos_sizeA);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMapSizeNLarge", mUserCageVerticesNos_sizeNLarge);//send the matrix to matlab*/

	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesAfterMap", mUserCageVerticesNos);//send the matrix to matlab

	std::string res = MatlabInterface::GetEngine().EvalToString("projectToLv_Lipman;");
	std::cerr << res;

	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mInternalPoints);//get the map from matlab
	cout.flush();
}

void SpaceDeformer2D::IncreaseVerteciesAfterMap(GMMDenseComplexColMatrix& OriginalCageVertecies, ComplexDoubleCPUMatrix& IncreasedCageVertecies, int numOfIncreasedCageVertecies, std::vector<int>& mNumOfVerticesInEdges) {

	int numOfOriginalVertecies = OriginalCageVertecies.size();
	if (numOfOriginalVertecies >= numOfIncreasedCageVertecies) {
		gmmToCpuMatrix(OriginalCageVertecies, IncreasedCageVertecies);
		return;
	}

	//IncreasedCageVertecies.resize(numOfIncreasedCageVertecies, 1);//move to preprocesing

	int index = 0;
	for (int i = 0; i < numOfOriginalVertecies; i++) {
		//insert the original vertex
		IncreasedCageVertecies(index++,0) = OriginalCageVertecies(i,0);

		//find the number of new veritecies per edge
		Complex p1 = OriginalCageVertecies((i + 1) % numOfOriginalVertecies,0);
		Complex p2 = OriginalCageVertecies(i,0);
		Complex vec = p1 - p2;
		double edgeLength = abs(p1 - p2);
		int numOfSegmentsPerEdge = mNumOfVerticesInEdges[i];

		//find the size of a segment in this edge
		double segmentLengthInEdge = edgeLength / numOfSegmentsPerEdge;

		//create a vector the size of a segment in the edge direction
		vec = vec / edgeLength;//normalize
		vec = vec * segmentLengthInEdge;

		//create new points 
		for (int j = 1; j < numOfSegmentsPerEdge; j++) {
			IncreasedCageVertecies(index++,0) = OriginalCageVertecies(i,0) + vec*Complex(j);
		}
	}

}
void SpaceDeformer2D::evalFzAndFzBar(GMMDenseComplexColMatrix& CageVerticesNos, GMMDenseComplexColMatrix& userCageVerticesNos, std::vector<int>& mNumOfVerticesInEdges, int numOfIncreasedCageVertecies, int numOfCageVertices, ComplexDoubleCPUMatrix& fz, ComplexDoubleCPUMatrix& fzBar){

	int index = 0;
	for (int i = 0; i < numOfCageVertices; i++){
		Complex ds = CageVerticesNos((i + 1) % numOfCageVertices,0) - CageVerticesNos(i,0);
		Complex dd = userCageVerticesNos((i + 1) % numOfCageVertices,0) - userCageVerticesNos(i,0);

		Complex  fzCur = 0.5*(abs(dd) + abs(ds))*dd / (abs(dd) * ds);
		Complex fzBarCur = 0.5*(abs(dd) - abs(ds))*dd / (abs(dd) * conj(ds));

		int numOfSegmentsPerEdge = mNumOfVerticesInEdges[i];
		for (int j = 0; j < numOfSegmentsPerEdge; j++){
			fz(index,0) = fzCur;
			fzBar(index++,0) = fzBarCur;
		}

	}
}
void SpaceDeformer2D::logarithmExtraction(GMMDenseComplexColMatrix& cageVerticesNos_sizeA, ComplexDoubleCPUMatrix& fz, ComplexDoubleCPUMatrix& cageVerteciesAfterMapSizeA, int a, ComplexDoubleCPUMatrix& log_fz){
#ifdef DEBUG
	std::vector<double> cornerAngleChanges; //debug
#endif // DEBUG

	double arg_fz;
	double ln_abs_fz;
	for (int i = 0; i < a; i++){
		Complex edgeBeforeMap = cageVerticesNos_sizeA((i + 1) % a,0) - cageVerticesNos_sizeA(i,0);
		Complex edgeAfterMap = cageVerteciesAfterMapSizeA((i + 1) % a,0) - cageVerteciesAfterMapSizeA(i,0);
		if (i == 0){
			arg_fz = arg(edgeAfterMap / edgeBeforeMap);
			//debug:
#ifdef DEBUG
			Complex prevEdgeBeforeMap = cageVerticesNos_sizeA(0,0) - cageVerticesNos_sizeA(a - 1,0);
			double cornerSourceAngle = arg(edgeBeforeMap / prevEdgeBeforeMap) + M_PI;
			double cornerTargetAngle = arg(edgeBeforeMap * fz(0,0) / (prevEdgeBeforeMap * fz(a - 1,0))) + M_PI;
			double cornerAngleChange = cornerTargetAngle - cornerSourceAngle;
			cornerAngleChanges.push_back(cornerAngleChange);
			//end debug
#endif // DEBUG

		}
		else{
			Complex prevEdgeBeforeMap = cageVerticesNos_sizeA(i,0) - cageVerticesNos_sizeA(i - 1,0);
			Complex prevfz = fz(i - 1,0);
			double cornerSourceAngle = arg(edgeBeforeMap / prevEdgeBeforeMap) + M_PI;
			double cornerTargetAngle = arg(edgeBeforeMap * fz(i,0) / (prevEdgeBeforeMap * prevfz)) + M_PI;
			double cornerAngleChange = cornerTargetAngle - cornerSourceAngle;
#ifdef DEBUG
			cornerAngleChanges.push_back(cornerAngleChange);//debug
#endif // DEBUG

			arg_fz = arg_fz + cornerAngleChange;
		}
		ln_abs_fz = log(abs(fz(i,0)));
		log_fz(i,0) = ln_abs_fz + Complex(0, 1)*arg_fz;
	}
#ifdef DEBUG
	//check turning number
	if (abs(std::accumulate(cornerAngleChanges.begin(), cornerAngleChanges.end(), 0.0))>1e-6){
		std::cerr << "The turning number of target polygon is invalid!\n";
	}
	//check log
	int max = 0;
	for (int i = 0; i < a; i++){
		double value = abs(exp(log_fz(i,0)) - exp(log(fz(i,0))));
		if (value > max)
			max = value;
	}
	if (max>1e-6)
		std::cerr << "The sanity check for log of derivative failed!\n";
#endif // DEBUG

}
void SpaceDeformer2D::find_nu_f(ComplexDoubleCPUMatrix& fz, ComplexDoubleCPUMatrix& fzBar, int a, ComplexDoubleCPUMatrix& nu_f){
	//GMMDenseComplexColMatrix nu_f(a, 1);//return this
	for (int i = 0; i < a; i++){
		nu_f(i,0) = std::conj(fzBar(i,0)) / fz(i,0);
	}
}
bool SpaceDeformer2D::checkIfInsidePolygon(double x,double y){
	/*//using cross product
	for (int i = 0; i < mXvaluesOfIntersections.size()-1; i++){
		double crossSegValue = (mXvaluesOfIntersections[i + 1] - mXvaluesOfIntersections[i])*(y - mYvaluesOfIntersections[i]) - (mYvaluesOfIntersections[i + 1] - mYvaluesOfIntersections[i])*(x - mXvaluesOfIntersections[i]);
		if (crossSegValue > 0-EPS){
			return false;
		}
	}
	return true;*/

	//using line equation
	return ((x <= k + epsilon) && (x <= log(SigmaA) - y + epsilon) && (log(sigmaB) + mSlopeOfLineApproxForCurveInLv*x <= y + epsilon));

}
void SpaceDeformer2D::projectPointToPolygonMinSeg(double& x, double& y){//x and y are changed in this function
	std::vector<double> closestPointsXvalues;
	std::vector<double> closestPointsYvalues;
	std::vector<double> minDistances;

	//using cross product
	for (int i = 0; i < mXvaluesOfIntersections.size() - 1; i++){
		double dot = (x - mXvaluesOfIntersections[i])*(mXvaluesOfIntersections[i + 1] - mXvaluesOfIntersections[i]) + (y - mYvaluesOfIntersections[i])*(mYvaluesOfIntersections[i + 1] - mYvaluesOfIntersections[i]);
		double projectionOnLine = dot / (pow((mXvaluesOfIntersections[i + 1] - mXvaluesOfIntersections[i]),2) + pow((mYvaluesOfIntersections[i + 1] - mYvaluesOfIntersections[i]),2));
		double t = std::max(0.0, min(1.0, projectionOnLine));
		double closestPointXvalue = mXvaluesOfIntersections[i] + t*(mXvaluesOfIntersections[i + 1] - mXvaluesOfIntersections[i]);
		double closestPointYvalue = mYvaluesOfIntersections[i] + t*(mYvaluesOfIntersections[i + 1] - mYvaluesOfIntersections[i]);
		closestPointsXvalues.push_back(closestPointXvalue);
		closestPointsYvalues.push_back(closestPointYvalue);
		minDistances.push_back(pow((x - closestPointXvalue),2) + pow((y - closestPointYvalue),2));
	}
	int indexOfMin = distance(minDistances.begin(), min_element(minDistances.begin(), minDistances.end()));
	x = closestPointsXvalues[indexOfMin];
	y = closestPointsYvalues[indexOfMin];

}

void SpaceDeformer2D::projectPointToPolygonWithK(double& x, double& y){

	if (y >= x + log(SigmaA)){
		x = 0;
		y = log(SigmaA);
	}
	else if (y > x + log(SigmaA) - 2 * k){
		double prevX = x;
		x = (x - y + log(SigmaA)) / 2;
		y = (y - prevX + log(SigmaA)) / 2;
	}
	else if (y >= log(SigmaA) - k){
		x = k;
		y = log(SigmaA);
	}
	else if (y > log(sigmaB / (1 - k))){
		x = k;
	}
	else if (y >= -x / mSlopeOfLineApproxForCurveInLv + log(sigmaB / (1 - k)) + k / mSlopeOfLineApproxForCurveInLv){
		x = k;
		y = log(sigmaB / (1 - k));
	}
	else if (y > -x / mSlopeOfLineApproxForCurveInLv + log(sigmaB)){
		double prevX = x;
		x = (x + mSlopeOfLineApproxForCurveInLv*y - mSlopeOfLineApproxForCurveInLv*log(sigmaB)) / (pow(mSlopeOfLineApproxForCurveInLv, 2) + 1);
		y = mSlopeOfLineApproxForCurveInLv*((prevX + mSlopeOfLineApproxForCurveInLv*y - mSlopeOfLineApproxForCurveInLv*log(sigmaB)) / (pow(mSlopeOfLineApproxForCurveInLv, 2) + 1)) + log(sigmaB);
	}
	else{
		x = 0;
		y = log(sigmaB);
	}

}

void SpaceDeformer2D::projectPointToPolygonNoK(double& x, double& y){

	if (y >= x + log(SigmaA)){
		x = 0;
		y = log(SigmaA);
	}
	else if (y > x + mYcoordOfIntersectionPointForCurveLv - mXcoordOfIntersectionPointForCurveLv){
		double prevX = x;
		x = (x - y + log(SigmaA)) / 2;
		y = (y - prevX + log(SigmaA)) / 2;
	}
	else if (y >= -x / mSlopeOfLineApproxForCurveInLv + mYcoordOfIntersectionPointForCurveLv + mXcoordOfIntersectionPointForCurveLv / mSlopeOfLineApproxForCurveInLv){
		x = mXcoordOfIntersectionPointForCurveLv;
		y = mYcoordOfIntersectionPointForCurveLv;
	}
	else if (y > -x / mSlopeOfLineApproxForCurveInLv + log(sigmaB)){
		double prevX = x;
		x = (x+mSlopeOfLineApproxForCurveInLv*y - mSlopeOfLineApproxForCurveInLv*log(sigmaB)) / (pow(mSlopeOfLineApproxForCurveInLv, 2) + 1);
		y = mSlopeOfLineApproxForCurveInLv*((prevX + mSlopeOfLineApproxForCurveInLv*y - mSlopeOfLineApproxForCurveInLv*log(sigmaB)) / (pow(mSlopeOfLineApproxForCurveInLv, 2) + 1)) + log(sigmaB);
	}
	else{
		x = 0;
		y = log(sigmaB);
	}
}
bool SpaceDeformer2D::localStep_noK_gpu(ComplexDoubleGPUMatrix& log_fz, ComplexDoubleGPUMatrix& nu_f){
	cuProjectPointsToPolygonNoK(log_fz.nPitchedRows()*log_fz.nCols(), (std::complex<double>*)log_fz.getData(), (std::complex<double>*)nu_f.getData(),log(SigmaA),sigmaB,k,mXcoordOfIntersectionPointForCurveLv,epsilon,mSlopeOfLineApproxForCurveInLv);
	//cuProjectPointsToPolygonNoK(mCurrentNumOfSegmentsA, (std::complex<double>*)log_fz.getData(), (std::complex<double>*)nu_f.getData(), log(SigmaA), sigmaB, k, mXcoordOfIntersectionPointForCurveLv, epsilon, mSlopeOfLineApproxForCurveInLv);

	return false;//fix
}
bool SpaceDeformer2D::localStep_withK_gpu(ComplexDoubleGPUMatrix& log_fz, ComplexDoubleGPUMatrix& nu_f){
	cuProjectPointsToPolygonWithK(log_fz.nPitchedRows()*log_fz.nCols(), (std::complex<double>*)log_fz.getData(), (std::complex<double>*)nu_f.getData(), log(SigmaA), sigmaB, k, epsilon, mSlopeOfLineApproxForCurveInLv);
	//cuProjectPointsToPolygonWithK(mCurrentNumOfSegmentsA, (std::complex<double>*)log_fz.getData(), (std::complex<double>*)nu_f.getData(), log(SigmaA), sigmaB, k, epsilon, mSlopeOfLineApproxForCurveInLv);
	return false;//fix
}
bool SpaceDeformer2D::localStep_noK_HP_gpu(ComplexDoubleGPUMatrix& x_vec){
	cuProjectPointsToPolygonNoK_HP(mCurrentNumOfSegmentsA, (std::complex<double>*)x_vec.getData(), log(SigmaA), sigmaB, k, mXcoordOfIntersectionPointForCurveLv, epsilon, mSlopeOfLineApproxForCurveInLv);//check
	return false;//fix
}
bool SpaceDeformer2D::localStep_withK_HP_gpu(ComplexDoubleGPUMatrix& x_vec){
	cuProjectPointsToPolygonWithK_HP(mCurrentNumOfSegmentsA, (std::complex<double>*)x_vec.getData(), log(SigmaA), sigmaB, k, epsilon, mSlopeOfLineApproxForCurveInLv);//check
	return false;//fix
}
bool SpaceDeformer2D::localStep_minSeg_gpu(ComplexDoubleGPUMatrix& log_fz, ComplexDoubleGPUMatrix& nu_f){
	cuProjectPointToPolygonMinSeg(log_fz.nPitchedRows()*log_fz.nCols(), (std::complex<double>*)log_fz.getData(), (std::complex<double>*)nu_f.getData(), mXvaluesOfIntersections_gpu.getData(), mYvaluesOfIntersections_gpu.getData(), mYvaluesOfIntersections_gpu.nRows() - 1,epsilon);
	return false;//fix
}
bool SpaceDeformer2D::localStep_minSeg_HP_gpu(ComplexDoubleGPUMatrix& x_vec){
	cuProjectPointToPolygonMinSeg_HP(mCurrentNumOfSegmentsA, (std::complex<double>*)x_vec.getData(), mXvaluesOfIntersections_gpu.getData(), mYvaluesOfIntersections_gpu.getData(), mYvaluesOfIntersections_gpu.nRows() - 1, epsilon);
	return false;//fix
}
bool SpaceDeformer2D::localStep(ComplexDoubleCPUMatrix& log_fz, ComplexDoubleCPUMatrix& nu_f, void (SpaceDeformer2D::*projectionFunction)(double&, double&)){

	bool allPointsInPolygon = true;
	for (int i = 0; i < log_fz.nRows(); i++){
		double x = abs(nu_f(i,0));
		double y = log_fz(i,0).real();
		if (checkIfInsidePolygon(x, y))
			continue;
		else{
			allPointsInPolygon = false;
			//projectPointToPolygonMinSeg(x, y);
			(this->*projectionFunction)(x, y);
			nu_f(i, 0) = x*exp(Complex(0, arg(nu_f(i, 0))));
			log_fz(i, 0) = Complex(y, log_fz(i, 0).imag());
		}
	}
	return allPointsInPolygon;
}
int SpaceDeformer2D::doLocalGlobalIterations(bool (SpaceDeformer2D::*localStepFunction)(ComplexDoubleGPUMatrix&, ComplexDoubleGPUMatrix&)){
#ifdef COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
	CPUTimer timer;
	double duration1 = 0.0, duration2 = 0.0;
#endif //COMPARE_LOCAL_VS_GLOBAL_TIMIMGS

	int i;
	for (i = 0; i < iterationsNum; i++){
#ifdef COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
		timer.tic();
#endif //COMPARE_LOCAL_VS_GLOBAL_TIMIMGS

		//local step
		//bool allPointsInPolygon = localStep(mLog_fz, mNu_f, projectionFunction);

		//bool allPointsInPolygon = (this->*localStepFunction)(mLog_fz_gpu, mNu_f_gpu);
		bool allPointsInPolygon = localStep_minSeg_gpu(mLog_fz_gpu, mNu_f_gpu);
#ifdef COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
		duration1 = timer.toc();
		cout << "local step took " << duration1 << endl;
#endif //COMPARE_LOCAL_VS_GLOBAL_TIMIMGS

		//stop condition -fix
		if (allPointsInPolygon)
			break;

#ifdef COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
		timer.tic();
#endif //COMPARE_LOCAL_VS_GLOBAL_TIMIMGS

		//global step
		mL_gpu.mult(mPinvOfIncCageVertexCoords_gpuMat, mLog_fz_gpu);
		mNu_gpu.mult(mPinvOfIncCageVertexCoords_gpuMat, mNu_f_gpu);

		mLog_fz_gpu.mult(mIncCageVertexCoords_gpuMat, mL_gpu);
		mNu_f_gpu.mult(mIncCageVertexCoords_gpuMat, mNu_gpu);

#ifdef COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
		duration2 = timer.toc();
		cout << "global step took " << duration2 << endl;
#endif //COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
	}

#ifdef DEBUG
	ComplexDoubleCPUMatrix test = mL_gpu;
	GMMDenseComplexColMatrix matlab0(test.nRows(), 1);
	for (int j = 0; j < test.nRows(); ++j)
		matlab0(j, 0) = test(j, 0);
	MatlabGMMDataExchange::SetEngineDenseMatrix("mL_gpu", matlab0);
	test = mNu_gpu;
	GMMDenseComplexColMatrix matlab1(test.nRows(), 1);
	for (int j = 0; j < test.nRows(); ++j)
		matlab1(j, 0) = test(j, 0);
	MatlabGMMDataExchange::SetEngineDenseMatrix("mNu_gpu", matlab1);
#endif
	/*cout << "cout: local-global iterations took " << durationn3 << endl;*/
	return i;
}

int SpaceDeformer2D::doHyperPlaneIterations(bool (SpaceDeformer2D::*localStepFunction)(ComplexDoubleGPUMatrix&)){
#ifdef COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
	CPUTimer timer;
	double duration1 = 0.0, duration2 = 0.0;
#endif //COMPARE_LOCAL_VS_GLOBAL_TIMIMGS

	double gpuTime;

	int i;
	Complex scalar;
	double norm_n0;
	for (i = 0; i < iterationsNum; i++){
		mX_local_gpu = mX_gpu;

		//local step
#ifdef COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
		timer.tic();
#endif //COMPARE_LOCAL_VS_GLOBAL_TIMIMGS

		//(this->*localStepFunction)(mX_local_gpu);
		localStep_minSeg_HP_gpu(mX_local_gpu);

#ifdef COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
		duration1 = timer.toc();
		cout << "local step took " << duration1 << endl;
		timer.tic();
#endif //COMPARE_LOCAL_VS_GLOBAL_TIMIMGS


		//global step
		mn_0forLipmansMethod_gpu.sub(mX_gpu, mX_local_gpu, &gpuTime);
		norm_n0 = mn_0forLipmansMethod_gpu.norm();
		//stop condition
		if (norm_n0 < epsilon)
			break;

#ifdef COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
		timer.tic();
#endif //COMPARE_LOCAL_VS_GLOBAL_TIMIMGS

		/*Edwad formula:*/
		mEta_0forLipmansMethod_gpu.mult(mTtrasposeForLipmansMethod_gpu, mn_0forLipmansMethod_gpu);
		my_etaforLipmansMethod_gpu.mult(mInvMForLipmansMethod_gpu, mEta_0forLipmansMethod_gpu);

		scalar=mEta_0forLipmansMethod_gpu.dotColVec(my_etaforLipmansMethod_gpu);
		scalar=pow(norm_n0,2)/scalar;
		my_etaforLipmansMethod_gpu.scale(scalar);

		mLnu_forLipman_gpu.sub(mLnu_forLipman_gpu, my_etaforLipmansMethod_gpu);

		mX_gpu.mult(mTForLipmansMethod_gpu, mLnu_forLipman_gpu);

#ifdef COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
		duration2 = timer.toc();
		cout << "global step took " << duration2 << endl;
#endif //COMPARE_LOCAL_VS_GLOBAL_TIMIMGS
	}
	return i;
}

void SpaceDeformer2D::calcIntegralUsingSpaningTree(ComplexDoubleCPUMatrix& f, ComplexDoubleCPUMatrix& f_tag, Complex phi_Z0){
	f(mZ0index,0) = phi_Z0;
	for (int i = 0; i < mEndIndicesForIntegral.size(); i++){
		Complex diffPHItagOnEdge = f_tag(mEndIndicesForIntegral[i] - 1, 0) + f_tag(mStartIndicesForIntegral[i] - 1, 0);
		Complex integralOnEdge = diffPHItagOnEdge*mEdgeVectorsForIntegral[i];
		f(mEndIndicesForIntegral[i] - 1,0) = f(mStartIndicesForIntegral[i] - 1,0) + integralOnEdge;
	}
}

void SpaceDeformer2D::findPHI(ComplexDoubleCPUMatrix& PHI, ComplexDoubleCPUMatrix& PHItag, ComplexDoubleCPUMatrix& LonInternalPoints){

	/*for (int i = 0; i < LonInternalPoints.nRows(); i++){
		PHItag(i, 0) = exp(LonInternalPoints(i, 0));
	}*/
	PHItag = LonInternalPoints;
	PHItag.exponent();

	IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeNLarge, mCurrentNLarge, mNumOfVerticesInEdgesSizeNlarge_stdVec);
	ComplexDoubleGPUMatrix phi_Z0_gpu(1, 1);
	phi_Z0_gpu.mult(mCauchyCoordsOfz0_gpuMat, mUserCageVerticesNos_sizeNLarge);
	ComplexDoubleCPUMatrix phi_Z0(phi_Z0_gpu);
	calcIntegralUsingSpaningTree(PHI, PHItag, phi_Z0(0, 0));
}

void SpaceDeformer2D::findPSI(ComplexDoubleCPUMatrix& PSI, ComplexDoubleCPUMatrix& NUonInternalPoints, ComplexDoubleCPUMatrix& PHItag){

	ComplexDoubleCPUMatrix PSItag(mNumOfInternalPoints, 1);

	for (int i = 0; i < NUonInternalPoints.nRows(); i++){
		PSItag(i, 0) = NUonInternalPoints(i, 0)*PHItag(i, 0);
	}

	Complex psi_Z0 = 0;
	calcIntegralUsingSpaningTree(PSI, PSItag, psi_Z0);
}

void SpaceDeformer2D::calcLvprojectionLGgpu(){

	//eval fz and fzBar
	IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeA, mCurrentNumOfSegmentsA, mNumOfVerticesInEdgesSizeA_stdVec);

	evalFzAndFzBar(mCompCageVerticesNos, mUserCageVerticesNos, mNumOfVerticesInEdgesSizeA_stdVec, mCurrentNumOfSegmentsA, mUserCageVerticesNos.size(), mfz, mfzBar);

	//extract argument from gz, and evaluate log(gz), Vg on A
	logarithmExtraction(mCompCageVerticesNos_sizeA, mfz, mUserCageVerticesNos_sizeA, mCurrentNumOfSegmentsA, mLog_fz);
	find_nu_f(mfz, mfzBar, mCurrentNumOfSegmentsA, mNu_f);

	mLog_fz_gpu = mLog_fz;
	mNu_f_gpu = mNu_f;

	mL_gpu.mult(mPinvOfIncCageVertexCoords_gpuMat, mLog_fz_gpu);
	mNu_gpu.mult(mPinvOfIncCageVertexCoords_gpuMat, mNu_f_gpu);
	
	mLog_fz_gpu.mult(mIncCageVertexCoords_gpuMat, mL_gpu);
	mNu_f_gpu.mult(mIncCageVertexCoords_gpuMat, mNu_gpu);

	bool (SpaceDeformer2D::*localStepFunction)(ComplexDoubleGPUMatrix&, ComplexDoubleGPUMatrix&) = NULL;
	if (k == mXcoordOfIntersectionPointForCurveLv)
		localStepFunction = &SpaceDeformer2D::localStep_withK_gpu;
	else
		localStepFunction = &SpaceDeformer2D::localStep_noK_gpu;

	doLocalGlobalIterations(localStepFunction);

	mLonInternalPoints_gpu.mult(mCauchyCoordinatesIncForP2P_gpuMat, mL_gpu);
	mLonInternalPoints = mLonInternalPoints_gpu;

	mNUonInternalPoints_gpu.mult(mCauchyCoordinatesIncForP2P_gpuMat, mNu_gpu);
	mNUonInternalPoints = mNUonInternalPoints_gpu;

	//find phi(z) and psi(z)

	findPHI(mPHI, mPHItag, mLonInternalPoints);

	findPSI(mPSI, mNUonInternalPoints, mPHItag);

	for (int i = 0; i < mInternalPoints.size(); i++){
		mInternalPoints[i] = mPHI(i,0) + std::conj(mPSI(i,0));
	}
#ifdef DEBUG_FINAL_RESULT
	MatlabGMMDataExchange::SetEngineDenseMatrix("mInternalPoints", mInternalPoints);
#endif
}

void SpaceDeformer2D::calcLvprojectionHPgpu(){

	//eval fz and fzBar
	IncreaseVerteciesAfterMap(mUserCageVerticesNos, mUserCageVerticesNos_sizeA, mCurrentNumOfSegmentsA, mNumOfVerticesInEdgesSizeA_stdVec);

	evalFzAndFzBar(mCompCageVerticesNos, mUserCageVerticesNos, mNumOfVerticesInEdgesSizeA_stdVec, mCurrentNumOfSegmentsA, mUserCageVerticesNos.size(), mfz, mfzBar);

	//extract argument from gz, and evaluate log(gz), Vg on A
	logarithmExtraction(mCompCageVerticesNos_sizeA, mfz, mUserCageVerticesNos_sizeA, mCurrentNumOfSegmentsA, mLog_fz);
	find_nu_f(mfz, mfzBar, mCurrentNumOfSegmentsA, mNu_f);

	mLog_fz_gpu = mLog_fz;
	mNu_f_gpu = mNu_f;

	mL_gpu.mult(mPinvOfIncCageVertexCoords_gpuMat, mLog_fz_gpu);
	mNu_gpu.mult(mPinvOfIncCageVertexCoords_gpuMat, mNu_f_gpu);

	mLog_fz_gpu.mult(mIncCageVertexCoords_gpuMat, mL_gpu);
	mNu_f_gpu.mult(mIncCageVertexCoords_gpuMat, mNu_gpu);

	mLnu_forLipman_gpu.concatenate(mL_gpu, mNu_gpu);
	mX_gpu.concatenate(mLog_fz_gpu, mNu_f_gpu);

#ifdef DEBUG2
	ComplexDoubleCPUMatrix test = mX_gpu;//debug
#endif

	bool (SpaceDeformer2D::*localStepFunction)(ComplexDoubleGPUMatrix&) = NULL;
	if (k == mXcoordOfIntersectionPointForCurveLv)
		localStepFunction = &SpaceDeformer2D::localStep_withK_HP_gpu;
	else
		localStepFunction = &SpaceDeformer2D::localStep_noK_HP_gpu;

	doHyperPlaneIterations(localStepFunction);

#ifdef DEBUG2
	test = mLnu_forLipman_gpu;//debug
	GMMDenseComplexColMatrix matlab1(test.nRows(), 1);
	for (int i = 0; i < test.nRows(); ++i)
		matlab1(i, 0) = test(i, 0);
	MatlabGMMDataExchange::SetEngineDenseMatrix("mLnu_forLipman_gpu", matlab1);
#endif

	mLnu_forLipman_gpu.splitInMiddle(mL_gpu, mNu_gpu);

#ifdef DEBUG2
	ComplexDoubleCPUMatrix mL = mL_gpu;//debug
	ComplexDoubleCPUMatrix mNu = mNu_gpu;//debug
#endif
	mLonInternalPoints_gpu.mult(mCauchyCoordinatesIncForP2P_gpuMat, mL_gpu);
	mLonInternalPoints = mLonInternalPoints_gpu;

	mNUonInternalPoints_gpu.mult(mCauchyCoordinatesIncForP2P_gpuMat, mNu_gpu);
	mNUonInternalPoints = mNUonInternalPoints_gpu;

	//find phi(z) and psi(z)

	findPHI(mPHI, mPHItag, mLonInternalPoints);

	findPSI(mPSI, mNUonInternalPoints, mPHItag);

	for (int i = 0; i < mInternalPoints.size(); i++){
		mInternalPoints[i] = mPHI(i, 0) + std::conj(mPSI(i, 0));
	}
#ifdef DEBUG_FINAL_RESULT
	MatlabGMMDataExchange::SetEngineDenseMatrix("mInternalPoints", mInternalPoints);
#endif
}

MStatus SpaceDeformer2D::showIncVertecies(MPointArray& IncreasedCageVertecies) {
	int facesNum = 1;
	int verticesNum = IncreasedCageVertecies.length();
	MIntArray polygonCounts(facesNum, verticesNum);

	MIntArray polygonConnects(verticesNum);//array of vertex connections for each polygon
	for (int i = 0; i < verticesNum; i++)
		polygonConnects[i] = i;

	MFnMesh reconstructed(mcageMesh);
	reconstructed.create(verticesNum, facesNum, IncreasedCageVertecies, polygonCounts, polygonConnects);

	//Update the reconstructed mesh to appear in the maya GUI
	MCHECKERROR(reconstructed.updateSurface(), "faild to update the surface.");

}


std::string SpaceDeformer2D::RelativeToFullPath(char* relPath) {
	//please enter path like this "\\matlab scripts\\inverse.m"
	char *libvar;
	libvar = getenv("DGP_CODE_DIR");
	std::string str1(libvar);
	std::string str2(relPath);
	return str1 + str2;
}
MStatus SpaceDeformer2D::getData(IN MDataBlock& block,OUT MObject& cageMesh,OUT MObject& p2pMesh) {
	MStatus stat;
	//mNeedToCalcNewLine = false;
	MDataHandle envData = block.inputValue(envelope, &stat);
	if (MS::kSuccess != stat) return stat;

	env = envData.asFloat();

	MDataHandle coordHandle = block.inputValue(mCoordinateTypeAttr, &stat);
	coordinateType = coordHandle.asLong();

	MDataHandle segAHandle = block.inputValue(mNumOfSegmentsAttr, &stat);
	mNumOfSegmentsA = segAHandle.asInt();//we can only change this value in the first time (doSetup)

	MDataHandle segNHandle = block.inputValue(mNlargeAttr, &stat);
	mNLarge = segNHandle.asInt();//we can only change this value in the first time (doSetup)

	double test;
	MDataHandle kHandle = block.inputValue(mkAttr, &stat);
	test = kHandle.asDouble();
	if (k != test) mNeedToCalcNewLine = true;
	k = test;

	MDataHandle sigmaAHandle = block.inputValue(mSigmaaAttr, &stat);
	test = sigmaAHandle.asDouble();
	if (SigmaA != test) mNeedToCalcNewLine = true;
	SigmaA = test;

	MDataHandle sigmaBHandle = block.inputValue(msigmabAttr, &stat);
	test = sigmaBHandle.asDouble();
	if (sigmaB != test) mNeedToCalcNewLine = true;
	sigmaB = test;

	MDataHandle z0Handle = block.inputValue(mZ0Attr, &stat);
	MAYA_float3& z0= z0Handle.asFloat3();
	mZ0NotOnMesh[0] = z0[0]; mZ0NotOnMesh[1] = z0[1]; mZ0NotOnMesh[2] = 0;

	MDataHandle lambdaHandle = block.inputValue(mlambdaAttr, &stat);
	test = lambdaHandle.asDouble();
	if (lambda != test) {
		lambda = test;
		MatlabGMMDataExchange::SetEngineDenseMatrix("lambda", doubleToGmmMat(this->lambda));//send value to matlab
	}

	MDataHandle iterHandle = block.inputValue(mIterAttr, &stat);
	test = iterHandle.asInt();
	if (iterationsNum != test) {
		iterationsNum = test;
		MatlabGMMDataExchange::SetEngineDenseMatrix("max_iterations", doubleToGmmMat(this->iterationsNum));//send value to matlab
	}

	MDataHandle epsHandle = block.inputValue(mEpsilonAttr, &stat);
	test = epsHandle.asDouble();
	if (epsilon != test) {
		epsilon = test;
		MatlabGMMDataExchange::SetEngineDenseMatrix("epsilon", doubleToGmmMat(this->epsilon));//send value to matlab
	}

	MDataHandle handle = block.inputValue(mCageAttr, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	cageMesh = handle.asMesh();

	MDataHandle p2phandle = block.inputValue(mCageP2pAttr, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	p2pMesh = p2phandle.asMesh();

	/*if (mNeedToCalcNewLine){
		//calcSegments();
		findLineApproximationForCurve();
	}*/

	return stat;
}


MStatus SpaceDeformer2D::deform(MDataBlock& block, MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex)
{
	CPUTimer timer, timer2;
	double duration, duration2;
	MStatus stat;
	//MObject cageMesh;
	MObject p2pMesh;
	stat = getData(block, mcageMesh, p2pMesh);
	MCHECKERROR(stat, "could not get the data block");
	if (mcageMesh.isNull() || p2pMesh.isNull()) {
		return stat;
	}
	//*****************
	//retriving the mesh from the data block
	MArrayDataHandle hInput = block.outputArrayValue(input, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat)
		stat = hInput.jumpToElement(multiIndex);
	CHECK_MSTATUS_AND_RETURN_IT(stat)
		MObject oInputGeom = hInput.outputValue().child(inputGeom).asMesh();
	MFnMesh fnInputMesh(oInputGeom);

	//*****************


	MFnMesh cageMeshFn(mcageMesh, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	MFnMesh p2pMeshFn(p2pMesh, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	updateCage(cageMeshFn);//populates mUserCageVertices
	updateControlPoints(p2pMeshFn);//populates mUserP2P

	if (mIsFirstTime)
	{
		//	MatlabInterface::GetEngine().Eval("clear");
		stat = doSetup(iter, cageMeshFn);
		runTimeDoSetup();
		preprocessingIntegral(fnInputMesh, oInputGeom);
		CHECK_MSTATUS_AND_RETURN_IT(stat);
		mIsFirstTime = false;
	}

	//compute the deformation of all the internal points. This is done by simply multiplying the coordinate matrix by the cage vertices vector
	switch (coordinateType)
	{
	case 0://Approximating cage-based complex Cauchy coordinates
		gmm::mult(mCauchyCoordinates, mUserCageVerticesNos, mInternalPoints);//multiply cauchy coordinates with the new cage vertices and insert it to the internal points
		break;
	case 1://Interpolating complex Cauchy coordinates
#ifdef CVX_INTERPOLATION //calculate using cvx
		matlabCalcNewVerticesForInterpolation();
#else//calculate using the inverse matrix
		gmm::mult(mCauchyCoordsOfOriginalCageVertices, mUserCageVerticesNos, mInterpolationGenCage_f);
#endif
		gmm::mult(mCauchyCoordinates, mInterpolationGenCage_f, mInternalPoints);//get the new internal points 
		break;
	case 2://Point2Point coordinates
		//runtime
		runTimeDoSetup();
		matlabCalcNewVerticesForP2P();
		gmm::mult(mCauchyCoordinatesIncForP2P, mP2PGenCageVertices_f, mInternalPoints);//get the new internal points 
		break;
	case 3://projection via H space 
		//runtime
		runTimeDoSetup();
		matlabCalcLforHprojection();
		break;
	case 4://projection via Lv space conformal
		runTimeDoSetup();
		matlabCalcLforLvprojectionConformal();
		break;
	case 5://projection via Lv space conformal accelerated
		runTimeDoSetup();
		matlabCalcLforLvprojectionConformalAccel();
		break;
	case 6://projection via Lv space 
		runTimeDoSetup();
		matlabCalcLforLvprojection_curve();
		break;
	case 7://projection via Lv space 
		matlabCalcLforLvprojection();
		break;
	case 8://projection via Lv space using local global
		runTimeDoSetup();
#ifdef COMPARE_TOTAL_TIMINGS
		timer.tic();
#endif //COMPARE_TOTAL_TIMINGS
		matlabCalcLforLvprojectionAccel();
#ifdef COMPARE_TOTAL_TIMINGS
		duration = timer.toc();
		cout << "local-global on matlab took " << duration << endl;
#endif //COMPARE_TOTAL_TIMINGS
		break;
	case 9://projection via Lv space using lipman's methos
		runTimeDoSetup();
#ifdef COMPARE_TOTAL_TIMINGS
		timer.tic();
		runTimeDoSetup();
#endif //COMPARE_TOTAL_TIMINGS
		matlabCalcLforLvprojectionLipman();
#ifdef COMPARE_TOTAL_TIMINGS
		duration = timer.toc();
		cout << "Hyperplane on matlab took " << duration << endl;
#endif //COMPARE_TOTAL_TIMINGS
		break;
	case 10://projection via Lv space using Dykstra
		runTimeDoSetup();
		matlabCalcLforLvprojectionDykstra();
		break;
	case 11://projection via Lv space using local global on gpu
		runTimeDoSetup();
#ifdef COMPARE_TOTAL_TIMINGS
		timer2.tic();
#endif //COMPARE_TOTAL_TIMINGS
		calcLvprojectionLGgpu();
#ifdef COMPARE_TOTAL_TIMINGS
		duration2 = timer2.toc();
		cout << "local-global on gpu took " << duration2 << endl;
#endif //COMPARE_TOTAL_TIMINGS
		break;
	case 12://projection via Lv space using lipman's method on gpu
		runTimeDoSetup();
#ifdef COMPARE_TOTAL_TIMINGS
		timer2.tic();
#endif //COMPARE_TOTAL_TIMINGS
		calcLvprojectionHPgpu();
#ifdef COMPARE_TOTAL_TIMINGS
		duration2 = timer2.toc();
		cout << "Hyperplane on gpu took " << duration2 << endl;
#endif //COMPARE_TOTAL_TIMINGS
	}
	///////////////////////////////
	///////////////////////////////

	//update the new deformed position of all the internal vertices
	for (iter.reset(); !iter.isDone(); iter.next())
	{
		int i = iter.index();

		MPoint pt = iter.position();


		///// add your code here //////
		///////////////////////////////
		//update c to be the deformed position of the i'th vertex
		Complex c = mInternalPoints[i];

		///////////////////////////////
		///////////////////////////////

		iter.setPosition(MPoint(c.real(), c.imag(), 0.0));
	}

	return stat;
}


MStatus SpaceDeformer2D::updateCage(MFnMesh& cageMeshFn)
{
	MStatus stat;

	int numFaces = cageMeshFn.numPolygons(&stat);

	assert(numFaces == 1);

	MIntArray vertexIndices;
	cageMeshFn.getPolygonVertices(0, vertexIndices);
	int numV = vertexIndices.length();
	assert(numV >= 3);

	gmm::clear(mUserCageVerticesNos);
	gmm::resize(mUserCageVerticesNos, numV, 1);

	MPointArray vertexArray;
	stat = cageMeshFn.getPoints(vertexArray);

	assert(numV == vertexArray.length());

	for (int i = 0; i < numV; i++)
	{
		MPoint p = vertexArray[i];
		Complex c(p[0], p[1]);
		mUserCageVerticesNos(i, 0) = c;
	}
	return MS::kSuccess;
}


MStatus SpaceDeformer2D::updateControlPoints(MFnMesh& cageMeshFn)
{
	MStatus stat;

	int numFaces = cageMeshFn.numPolygons(&stat);

	assert(numFaces == 1);

	MIntArray vertexIndices;
	cageMeshFn.getPolygonVertices(0, vertexIndices);
	int numV = vertexIndices.length();
//	assert(numV >= 3);

	gmm::clear(mUserP2P);
	gmm::resize(mUserP2P, numV, 1);

	MPointArray vertexArray;
	stat = cageMeshFn.getPoints(vertexArray);

	assert(numV == vertexArray.length());

	for (int i = 0; i < numV; i++)
	{
		MPoint p = vertexArray[i];
		Complex c(p[0], p[1]);
		mUserP2P(i, 0) = c;
	}
	return MS::kSuccess;
}

void populateC(GMMDenseComplexColMatrix &C, Complex* cage, int n, MPointArray& constraints, int m) {
	for (int i = 0; i < m; i++) {
		MPoint pt = constraints[i];
		Complex z(pt[0], pt[1]); //cage point
		for (int j = 0; j < n; j++) {
			Complex K(0.0, 0.0);
			//update K to be value of the j'th coordinate at the i'th cage point
			int next = j + 1;
			int prev = j - 1;
			if (j == n - 1)
				next = 0;
			if (j == 0)
				prev = n - 1;

			Complex Zj = cage[j];
			Complex Zjnext = cage[next];
			Complex Zjprev = cage[prev];

			//the cauchy-green complex barycentric coordinates equesion

			Complex multiplicand1(0, (1 / (-2 * M_PI)));
			Complex multiplicand2 = ((Zjnext - z) / (Zjnext - Zj))*log((Zjnext - z) / (Zj - z)) -
				((Zjprev - z) / (Zj - Zjprev))*log((Zj - z) / (Zjprev - z));
			K = multiplicand1*multiplicand2;

			C(i, j) = K;
		}
	}
}

void populateD(GMMDenseComplexColMatrix &D, Complex* cage, int n, MPointArray& constraints, int m) {
	for (int i = 0; i < m; i++) {
		MPoint pt = constraints[i];
		Complex z(pt[0], pt[1]); //cage point
		for (int j = 0; j < n; j++) {
			Complex K(0.0, 0.0);
			//update K to be value of the j'th coordinate at the i'th cage point
			int next = j + 1;
			int prev = j - 1;
			if (j == n - 1)
				next = 0;
			if (j == 0)
				prev = n - 1;

			Complex Zj = cage[j];
			Complex Zjnext = cage[next];
			Complex Zjprev = cage[prev];

			//the cauchy-green complex barycentric coordinates equesion
			Complex multiplicand1(0, (1 / (-2 * M_PI)));
			Complex multiplicand2 = (((Complex)1 / ((Zjprev - z)*(Zj - z))) - ((Complex)1 / ((Zj - z)*(Zjnext - z))));


			K = multiplicand1*multiplicand2;

			D(i, j) = K;
		}
	}
}

void populateCtag(GMMDenseComplexColMatrix &D, Complex* cage, int n, MPointArray& constraints, int m) {
	for (int i = 0; i < m; i++) {
		MPoint pt = constraints[i];
		Complex z(pt[0], pt[1]); //cage point
		for (int j = 0; j < n; j++) {
			Complex K(0.0, 0.0);
			//update K to be value of the j'th coordinate at the i'th cage point
			int next = j + 1;
			int prev = j - 1;
			if (j == n - 1)
				next = 0;
			if (j == 0)
				prev = n - 1;

			Complex Zj = cage[j];
			Complex Zjnext = cage[next];
			Complex Zjprev = cage[prev];

			//the cauchy-green complex barycentric coordinates equesion
			Complex multiplicand1(0, (1 / (-2 * M_PI)));
			Complex multiplicand2 = (((Complex)1 / (Zjnext - Zj))*log((Zj - z) / (Zjnext - z)) +
									((Complex)1 / (Zj - Zjprev))*log((Zj - z) / (Zjprev - z)));


			K = multiplicand1*multiplicand2;

			D(i, j) = K;
		}
	}
}

void SpaceDeformer2D::IncreaseVertecies(MPointArray& OriginalCageVertecies, MPointArray& IncreasedCageVertecies,int numOfIncreasedCageVertecies, bool sizeAorNlarge) {

	//find the circumference of the cage polygon
	double circumference=0;
	int numOfOriginalVertecies = OriginalCageVertecies.length();
	IncreasedCageVertecies.clear();
	if (numOfOriginalVertecies >= numOfIncreasedCageVertecies) {
		IncreasedCageVertecies = OriginalCageVertecies;
		return;
	}

	for (int i = 0; i < numOfOriginalVertecies; i++) {

		circumference += (OriginalCageVertecies[(i + 1)% numOfOriginalVertecies].distanceTo(OriginalCageVertecies[i]));
	}

	//find the segment length
	double segmentLength = circumference / (double)numOfIncreasedCageVertecies;

	for (int i = 0; i < numOfOriginalVertecies; i++) {
		//insert the original vertex
		IncreasedCageVertecies.append(OriginalCageVertecies[i]);

		//find the number of new veritecies per edge
		double edgeLength = (OriginalCageVertecies[(i + 1) % numOfOriginalVertecies].distanceTo(OriginalCageVertecies[i]));
		int numOfSegmentsPerEdge = round(edgeLength / segmentLength);//the number of vertices in the edge is #seg-1
		numOfSegmentsPerEdge = std::max(1, numOfSegmentsPerEdge);

		//find the size of a segment in this edge
		double segmentLengthInEdge = edgeLength / numOfSegmentsPerEdge;

		//create a vector the size of a segment in the edge direction
		MVector vec = OriginalCageVertecies[(i + 1) % numOfOriginalVertecies] - (OriginalCageVertecies[i]);
		vec.normalize();
		vec = vec * segmentLengthInEdge;

		//create new points 
		for (int j = 1; j < numOfSegmentsPerEdge; j++) {
			IncreasedCageVertecies.append((MPoint)((OriginalCageVertecies[i])+ j*vec));
		}

		if (sizeAorNlarge) {
			mNumOfVerticesInEdgesSizeA[i] = numOfSegmentsPerEdge;
		}
		else
			mNumOfVerticesInEdgesSizeNlarge[i] = numOfSegmentsPerEdge;
	}

}

void SpaceDeformer2D::IncreaseVertecies(Complex* OriginalCompCageVertecies, int OrigCageSize, Complex** IncreasedCompCageVertecies, int& numOfIncreasedCageVertecies) {
	MPointArray OriginalCartCageVertecies;
	MPointArray IncreasedCageVertecies;

	for (int i = 0; i < OrigCageSize; i++) {
		Complex c = OriginalCompCageVertecies[i];
		MPoint p(c.real(), c.imag());
		OriginalCartCageVertecies.append(p);
	}

	IncreaseVertecies(OriginalCartCageVertecies, IncreasedCageVertecies, numOfIncreasedCageVertecies, false);

	int length1 = OriginalCartCageVertecies.length();
	numOfIncreasedCageVertecies = IncreasedCageVertecies.length();

	*IncreasedCompCageVertecies = new Complex[numOfIncreasedCageVertecies];

	for (int i = 0; i < numOfIncreasedCageVertecies; i++) {
		Complex c(IncreasedCageVertecies[i].x, IncreasedCageVertecies[i].y);
		(*IncreasedCompCageVertecies)[i] = c;
	}

	//*****temp for test*****
	showIncVertecies(IncreasedCageVertecies);

}

MStatus SpaceDeformer2D::doSetup(MItGeometry& iter, MFnMesh& cageMeshFn)
{
	MStatus stat;
//	MatlabInterface::GetEngine().EvalToString("addpath(genpath('C:/Users/Ben-PC/Documents/MySWprojects/ProjectFinal/DGP-HW1/DGP_CODE_DIR'))");

	mNumOfInternalPoints = iter.count(&stat); //num internal points (point of the triangulated cage)
	mNumOfCageVerticies = mUserCageVerticesNos.nrows(); //num of cage vertices
	mNumOfControlPoints = mUserP2P.nrows(); //num of control points



	gmm::clear(mCauchyCoordinates);
	gmm::resize(mCauchyCoordinates, mNumOfInternalPoints, mNumOfCageVerticies);

	gmm::clear(mInternalPoints);
	gmm::resize(mInternalPoints, mNumOfInternalPoints, 1);

	gmm::clear(mCauchyCoordsOfOriginalCageVertices);
	gmm::resize(mCauchyCoordsOfOriginalCageVertices, mNumOfCageVerticies, mNumOfCageVerticies);

	gmm::clear(mInterpolationGenCage_f);
	gmm::resize(mInterpolationGenCage_f, mNumOfCageVerticies, 1);


	//offset the cage by epsilon in the direction of the normal such that the triangle mesh (the image) is strictly inside the new cage
	cageMeshFn.getPoints(mCartCageVerticesNos);//extract the cage vertecies from the mesh object
	mCompCageVerticesNos = compPointArrayToGmmMat(mCartCageVerticesNos);
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesB4Map", mCompCageVerticesNos);//send the matrix to matlab
	this->mCompCageVerticesWos = new Complex[mNumOfCageVerticies];//complex coordinates

	//iterate through the cage vertecies and move all vertecies eps outwards
	for (int i = 0; i < mNumOfCageVerticies; i++)
	{
		int next = i + 1;
		int prev = i - 1;
		if (i == mNumOfCageVerticies - 1)
			next = 0;
		if (i == 0)
			prev = mNumOfCageVerticies - 1;


		MVector cageVertexNormal;//the normal to the vertex

		//Creating the neighbor edges to the cage vertex
		MVector Edge1 = mCartCageVerticesNos[i] - mCartCageVerticesNos[prev];
		MVector Edge2 = mCartCageVerticesNos[i] - mCartCageVerticesNos[next];
		Edge1.normalize();
		Edge2.normalize();

		//the notmal to the vertes is the angle bisector of the two edges
		float V = (Edge1^Edge2).z;
		if (IS_NUMERIC_ZERO(V)) {//no turn
			cageVertexNormal = MVector(-Edge2.y, Edge2.x);
		}else if (V>0) {//turn right
			cageVertexNormal = -(Edge1 + Edge2);
		}else {//turn left
			cageVertexNormal = Edge1 + Edge2;
		}
		
		cageVertexNormal.normalize();
		cageVertexNormal *= EPS; //epsilon to prevent points to be exactly on cage

		//move the vertex,the size of epsilon in the direction of the normal
		MVector offset = mCartCageVerticesNos[i] + cageVertexNormal;
		mCompCageVerticesWos[i] = Complex(offset.x, offset.y);
	}

	//calculate the interpolation coordinates to find the new vertices
	
	populateC(mCauchyCoordsOfOriginalCageVertices, mCompCageVerticesWos, mNumOfCageVerticies, mCartCageVerticesNos, mNumOfCageVerticies);

#ifndef CVX_INTERPOLATION

	//find the coordinates that will give us the new vertices so that the real cage point will look like interpolation
	//we will do it using cvx
	MatlabGMMDataExchange::SetEngineDenseMatrix("toInverse", mCauchyCoordsOfOriginalCageVertices);//send the matrix to matlab
	//load the matlab script
//	int res = MatlabInterface::GetEngine().LoadAndRunScript("%DGP_CODE_DIR%/matlab scripts/interpolatedCauchy.m"); -----not working
	
	int res = MatlabInterface::GetEngine().LoadAndRunScript(RelativeToFullPath("\\matlab scripts\\inverse.m").c_str());
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'interpolatedCauchy.m' failed with error code " << res << std::endl;
	}
	MatlabGMMDataExchange::GetEngineDenseMatrix("toInverse", mCauchyCoordsOfOriginalCageVertices);//get the incersed matrix from matlab
#endif


	//calculate Cj(z) the Cauchy-Green complex barycentric coordinates
	for (iter.reset(); !iter.isDone(); iter.next())
	{
		this->mInternalPoints_MPoint.append(iter.position());//fill the MPoint array
		//fill the complex array
		Complex c(iter.position().x, iter.position().y);
		int ind = iter.index();
		mInternalPoints(iter.index(), 0) = c;

	}

	populateC(mCauchyCoordinates, mCompCageVerticesWos, mNumOfCageVerticies, this->mInternalPoints_MPoint, mNumOfInternalPoints);

	//runTimeDoSetup();
	for (int i = 0; i < mNumOfControlPoints; i++) {
		Complex c = mUserP2P(i, 0);
		MPoint p(c.real(), c.imag());
		mInitialcontrolPoints.append(p);
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix("lambda", doubleToGmmMat(this->lambda));//send value to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("max_iterations", doubleToGmmMat(this->iterationsNum));//send value to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("epsilon", doubleToGmmMat(this->epsilon));//send value to matlab

	return MS::kSuccess;
}



int SpaceDeformer2D::findClosestInternalPointsToZ0() {
	MPoint Z0(mZ0NotOnMesh[0], mZ0NotOnMesh[1], mZ0NotOnMesh[2],0);
	double minDist = DBL_MAX;
	MPoint z0Point;
	int index=0;

	for (int i = 0; i < mInternalPoints_MPoint.length(); i++ ) {
		MPoint tempPoint = mInternalPoints_MPoint[i];
		double tempDist = sqrt(pow((Z0.x - tempPoint.x) ,2) + pow((Z0.y - tempPoint.y) ,2));

		if (tempDist < minDist) {
			z0Point = tempPoint;
			minDist = tempDist;
			index = i;
		}
	}
	mZ0onMesh = Complex(z0Point.x, z0Point.y);
	return index;
}

MStatus SpaceDeformer2D::runTimeDoSetup() {
	if (mNeedToCalcNewLine){
		findLineApproximationForCurve();
	}
	if (mNumOfSegmentsAOld == mNumOfSegmentsA && mNLargeOld == mNLarge) {
		return MS::kSuccess;
	}
	mNumOfSegmentsAOld = mNumOfSegmentsA;
	mNLargeOld = mNLarge;

	mZ0index = findClosestInternalPointsToZ0();

	Complex* IncreasedCompCageVertecies = NULL;// = new Complex[nLarge];

	gmm::clear(mNumOfVerticesInEdgesSizeNlarge);
	gmm::resize(mNumOfVerticesInEdgesSizeNlarge, mCartCageVerticesNos.length(), 1);

	IncreaseVertecies(IN mCompCageVerticesWos, IN mNumOfCageVerticies, OUT &IncreasedCompCageVertecies, mNLarge);//nLarge might change in this func
	mCurrentNLarge = mNLarge;

	MatlabGMMDataExchange::SetEngineDenseMatrix("NumOfVerticesInEdgesSizeNlarge", mNumOfVerticesInEdgesSizeNlarge);//send the matrix to matlab

	gmmToStdIntVector(mNumOfVerticesInEdgesSizeNlarge, mNumOfVerticesInEdgesSizeNlarge_stdVec);


	gmm::clear(mCauchyCoordinatesIncForP2P);
	gmm::resize(mCauchyCoordinatesIncForP2P, mNumOfInternalPoints, mNLarge);

	populateC(OUT mCauchyCoordinatesIncForP2P, IncreasedCompCageVertecies, mNLarge, this->mInternalPoints_MPoint, mNumOfInternalPoints);

	MatlabGMMDataExchange::SetEngineDenseMatrix("C_sizeM", mCauchyCoordinatesIncForP2P);//send the matrix to matlab

	ComplexDoubleCPUMatrix mCauchyCoordinatesIncForP2P_cpuMat(mCauchyCoordinatesIncForP2P.nrows(), mCauchyCoordinatesIncForP2P.ncols());
	gmmToCpuMatrix(mCauchyCoordinatesIncForP2P, mCauchyCoordinatesIncForP2P_cpuMat);
	mCauchyCoordinatesIncForP2P_gpuMat = mCauchyCoordinatesIncForP2P_cpuMat;

	MatlabGMMDataExchange::SetEngineDenseMatrix("Z0index", doubleToGmmMat((this->mZ0index) + 1));//send the matrix to matlab
	gmm::clear(mCauchyCoordsOfz0);
	gmm::resize(mCauchyCoordsOfz0, 1, mNLarge);
	std::string res = MatlabInterface::GetEngine().EvalToString("Cz0=C_sizeM(Z0index,:);");
	std::cerr << res << std::endl;

	MatlabGMMDataExchange::GetEngineDenseMatrix("Cz0", mCauchyCoordsOfz0);

	ComplexDoubleCPUMatrix mCauchyCoordsOfz0_cpuMat(mCauchyCoordsOfz0.nrows(), mCauchyCoordsOfz0.ncols());
	gmmToCpuMatrix(mCauchyCoordsOfz0, mCauchyCoordsOfz0_cpuMat);
	mCauchyCoordsOfz0_gpuMat = mCauchyCoordsOfz0_cpuMat;

	//*******************************************************************************
	//calculate D & C matrices for P2P

	gmm::clear(mCauchyCoordsOfOriginalP2P);
	gmm::resize(mCauchyCoordsOfOriginalP2P, mNumOfControlPoints, mNLarge);
	populateC(mCauchyCoordsOfOriginalP2P, IncreasedCompCageVertecies, mNLarge, mInitialcontrolPoints, mNumOfControlPoints);
	
	gmm::clear(mNumOfVerticesInEdgesSizeA);
	gmm::resize(mNumOfVerticesInEdgesSizeA, mCartCageVerticesNos.length(), 1);

	//MPointArray mCartCageVerticesNos_sizeA; //dimensions are: a x 1
	int x = mCartCageVerticesNos_sizeA.length();
	IncreaseVertecies(IN mCartCageVerticesNos, OUT mCartCageVerticesNos_sizeA, mNumOfSegmentsA,true);

	MatlabGMMDataExchange::SetEngineDenseMatrix("NumOfVerticesInEdgesSizeA", mNumOfVerticesInEdgesSizeA);//send the matrix to matlab
	mCompCageVerticesNos_sizeA = compPointArrayToGmmMat(mCartCageVerticesNos_sizeA);

	MatlabGMMDataExchange::SetEngineDenseMatrix("cageVerteciesB4Map_sizeA", mCompCageVerticesNos_sizeA);//send the matrix to matlab

	mCurrentNumOfSegmentsA = mCartCageVerticesNos_sizeA.length();

	gmmToStdIntVector(mNumOfVerticesInEdgesSizeA, mNumOfVerticesInEdgesSizeA_stdVec);

	gmm::clear(mSecondDerOfIncCageVertexCoords);
	gmm::resize(mSecondDerOfIncCageVertexCoords, mCurrentNumOfSegmentsA, mNLarge);

	populateD(OUT mSecondDerOfIncCageVertexCoords, IncreasedCompCageVertecies, mNLarge, mCartCageVerticesNos_sizeA, mCurrentNumOfSegmentsA);

	//********************************************************************************
	//calculate the first derivative of the cauchy grenn coords for the projection to H 
	//the matrix Ctag need to be a x nLarge

	gmm::clear(mFirstDerOfIncCageVertexCoords);
	gmm::resize(mFirstDerOfIncCageVertexCoords, mCurrentNumOfSegmentsA, mNLarge);

	gmm::clear(mIncCageVertexCoords);
	gmm::resize(mIncCageVertexCoords, mCurrentNumOfSegmentsA, mNLarge);

	populateC(OUT mIncCageVertexCoords, IncreasedCompCageVertecies, mNLarge, mCartCageVerticesNos_sizeA, mCurrentNumOfSegmentsA);
	MatlabGMMDataExchange::SetEngineDenseMatrix("C_sizeA", mIncCageVertexCoords);//send the matrix to matlab
	populateCtag(OUT mFirstDerOfIncCageVertexCoords, IncreasedCompCageVertecies, mNLarge, mCartCageVerticesNos_sizeA, mCurrentNumOfSegmentsA);

	ComplexDoubleCPUMatrix mIncCageVertexCoords_cpuMat(mIncCageVertexCoords.nrows(), mIncCageVertexCoords.ncols());
	gmmToCpuMatrix(mIncCageVertexCoords, mIncCageVertexCoords_cpuMat);
	mIncCageVertexCoords_gpuMat = mIncCageVertexCoords_cpuMat;
	
	//**********************************************
	// calc pinv for project to Lv accel
	gmm::clear(mPinvOfIncCageVertexCoords);
	gmm::resize(mPinvOfIncCageVertexCoords, mNLarge, mCurrentNumOfSegmentsA);
	
	res = MatlabInterface::GetEngine().EvalToString("p_inv=pinv(C_sizeA);");
	std::cerr << res << std::endl;

	MatlabGMMDataExchange::GetEngineDenseMatrix("p_inv", mPinvOfIncCageVertexCoords);

	ComplexDoubleCPUMatrix mPinvOfIncCageVertexCoords_cpuMat(mPinvOfIncCageVertexCoords.nrows(), mPinvOfIncCageVertexCoords.ncols());
	gmmToCpuMatrix(mPinvOfIncCageVertexCoords, mPinvOfIncCageVertexCoords_cpuMat);
	mPinvOfIncCageVertexCoords_gpuMat = mPinvOfIncCageVertexCoords_cpuMat;


	//**********************************************
	// calc LU prefactorization of matrix for lipman's method
	gmm::clear(mLMatrixForLipmansMethod); 
	gmm::resize(mLMatrixForLipmansMethod, 2 * mNLarge, 2 * mNLarge);
	gmm::clear(mUMatrixForLipmansMethod);
	gmm::resize(mUMatrixForLipmansMethod, 2 * mNLarge, 2 * mNLarge);
	gmm::clear(mTtrasposeForLipmansMethod);
	gmm::resize(mTtrasposeForLipmansMethod, 2 * mNLarge, 2 * mCurrentNumOfSegmentsA);
	gmm::clear(mTForLipmansMethod);
	gmm::resize(mTForLipmansMethod, 2 * mCurrentNumOfSegmentsA, 2 * mNLarge);

	res = MatlabInterface::GetEngine().EvalToString("T=blkdiag(C_sizeA,C_sizeA);T_trans=T';M=T_trans*T;[M_L, M_U] = lu(M); ");
	std::cerr << res << std::endl;

	MatlabGMMDataExchange::GetEngineSparseMatrix("M_L", mLMatrixForLipmansMethod);
	MatlabGMMDataExchange::GetEngineSparseMatrix("M_U", mUMatrixForLipmansMethod);
	MatlabGMMDataExchange::GetEngineDenseMatrix("T_trans", mTtrasposeForLipmansMethod);
	MatlabGMMDataExchange::GetEngineDenseMatrix("T", mTForLipmansMethod);

	ComplexDoubleCPUMatrix mTtrasposeForLipmansMethod_cpuMat(mTtrasposeForLipmansMethod.nrows(), mTtrasposeForLipmansMethod.ncols());
	gmmToCpuMatrix(mTtrasposeForLipmansMethod, mTtrasposeForLipmansMethod_cpuMat);
	mTtrasposeForLipmansMethod_gpu = mTtrasposeForLipmansMethod_cpuMat;

	ComplexDoubleCPUMatrix mTForLipmansMethod_cpuMat(mTForLipmansMethod.nrows(), mTForLipmansMethod.ncols());
	gmmToCpuMatrix(mTForLipmansMethod, mTForLipmansMethod_cpuMat);
	mTForLipmansMethod_gpu = mTForLipmansMethod_cpuMat;

	//calc the inverse matrix of M - need to decide if to use inv, chol or lu
	gmm::clear(mInvMForLipmansMethod);
	gmm::resize(mInvMForLipmansMethod, 2 * mNLarge, 2 * mNLarge);

	res = MatlabInterface::GetEngine().EvalToString("M_inv=inv(M); ");
	std::cerr << res << std::endl;

	MatlabGMMDataExchange::GetEngineDenseMatrix("M_inv", mInvMForLipmansMethod);

	ComplexDoubleCPUMatrix mInvMForLipmansMethod_cpuMat(mInvMForLipmansMethod.nrows(), mInvMForLipmansMethod.ncols());
	gmmToCpuMatrix(mInvMForLipmansMethod, mInvMForLipmansMethod_cpuMat);
	mInvMForLipmansMethod_gpu = mInvMForLipmansMethod_cpuMat;

	//***************************
	//initialize array sizes
	mUserCageVerticesNos_sizeA.resize(mCurrentNumOfSegmentsA, 1);
	mUserCageVerticesNos_sizeNLarge.resize(mCurrentNLarge, 1);
	mfz.resize(mCurrentNumOfSegmentsA, 1);
	mfzBar.resize(mCurrentNumOfSegmentsA, 1);
	mLog_fz.resize(mCurrentNumOfSegmentsA, 1);
	mNu_f.resize(mCurrentNumOfSegmentsA, 1);
	mLog_fz_gpu.resize(mCurrentNumOfSegmentsA, 1);
	mNu_f_gpu.resize(mCurrentNumOfSegmentsA, 1);
	mL_gpu.resize(mCurrentNLarge, 1);
	mNu_gpu.resize(mCurrentNLarge, 1);

	mLonInternalPoints.resize(mNumOfInternalPoints, 1);
	mNUonInternalPoints.resize(mNumOfInternalPoints, 1);
	mLonInternalPoints_gpu.resize(mNumOfInternalPoints, 1);
	mNUonInternalPoints_gpu.resize(mNumOfInternalPoints, 1);

	mPHI.resize(mNumOfInternalPoints, 1);
	mPHItag.resize(mNumOfInternalPoints, 1);
	mPSI.resize(mNumOfInternalPoints, 1);

	mX_gpu.resize(2*mCurrentNumOfSegmentsA, 1);
	mX_local_gpu.resize(2 * mCurrentNumOfSegmentsA, 1);
	mn_0forLipmansMethod_gpu.resize(2 * mCurrentNumOfSegmentsA, 1);
	mEta_0forLipmansMethod_gpu.resize(2 * mCurrentNLarge, 1);
	mc_0forLipmansMethod_gpu.resize(2 * mCurrentNLarge, 1);
	my_cforLipmansMethod_gpu.resize(2 * mCurrentNLarge, 1);
	my_etaforLipmansMethod_gpu.resize(2 * mCurrentNLarge, 1);
	
	mLnu_forLipman_gpu.resize(2 * mCurrentNLarge, 1);
	//*********************

	cout.flush();

	delete[] IncreasedCompCageVertecies;
	return MS::kSuccess;
}
MStatus SpaceDeformer2D::preprocessingIntegral(MFnMesh& inputMesh, MObject InputGeom) {
	MStatus stat;
	int numEdges1 = inputMesh.numEdges();
	int numVertices1 = inputMesh.numVertices();
	
	//create adjacency matrix
	GMMSparseRowMatrix adjacencyMatrix(mNumOfInternalPoints, mNumOfInternalPoints);
	MItMeshEdge edgesItr(InputGeom);//edge iterator for the mesh
	int numEdges = edgesItr.count();

	while (!edgesItr.isDone()) {
		MAYA_int2 connectedvertices;
		inputMesh.getEdgeVertices(edgesItr.index(), connectedvertices);

		adjacencyMatrix(connectedvertices[0], connectedvertices[1]) = 1;
		adjacencyMatrix(connectedvertices[1], connectedvertices[0]) = 1;
		edgesItr.next();

	}

	//send values to matlab
	MatlabGMMDataExchange::SetEngineSparseMatrix("adjacencyGraph", adjacencyMatrix);
	MatlabGMMDataExchange::SetEngineDenseMatrix("vertices", mInternalPoints);
	MatlabGMMDataExchange::SetEngineDenseMatrix("Z0index", doubleToGmmMat((this->mZ0index)+1));//send the matrix to matlab

	int res = MatlabInterface::GetEngine().LoadAndRunScript(RelativeToFullPath("\\matlab scripts\\preprocessingIntegral.m").c_str());
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'preprocessingIntegral.m' failed with error code " << res << std::endl;
	}

	gmm::clear(mEndIndicesForIntegral);
	gmm::resize(mEndIndicesForIntegral, mNumOfInternalPoints-1, 1);
	gmm::clear(mStartIndicesForIntegral);
	gmm::resize(mStartIndicesForIntegral, mNumOfInternalPoints - 1, 1);	
	gmm::clear(mEdgeVectorsForIntegral);
	gmm::resize(mEdgeVectorsForIntegral, mNumOfInternalPoints - 1, 1);

	MatlabGMMDataExchange::GetEngineDenseMatrix("endIndices_double", mEndIndicesForIntegral);//get the map from matlab
	MatlabGMMDataExchange::GetEngineDenseMatrix("startIndices_double", mStartIndicesForIntegral);//get the map from matlab
	MatlabGMMDataExchange::GetEngineDenseMatrix("edgeVectors", mEdgeVectorsForIntegral);//get the map from matlab


	//MatlabGMMDataExchange::GetEngineDenseMatrix("adjacencyGraph", mCauchyCoordsOfOriginalCageVertices);//get the incersed matrix from matlab
	return stat;
}
//unused for now- maybe change to one line calculation
MStatus SpaceDeformer2D::calcSegments(){
	const int numOfSegments = 50;
	gmm::clear(mAOfLineSegmentInLvAccelerated);
	gmm::resize(mAOfLineSegmentInLvAccelerated, numOfSegments, 1);
	gmm::clear(mBOfLineSegmentInLvAccelerated);
	gmm::resize(mBOfLineSegmentInLvAccelerated, numOfSegments, 1);

	MatlabGMMDataExchange::SetEngineDenseMatrix("k", doubleToGmmMat(this->k));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("SIGMA", doubleToGmmMat(this->SigmaA));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("sigma", doubleToGmmMat(this->sigmaB));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("numOfSegments", doubleToGmmMat(numOfSegments));//send the matrix to matlab

	std::string res = MatlabInterface::GetEngine().EvalToString("[A,B]=createLineSegments(sigma,SIGMA,k,numOfSegments+1);");
	std::cout << res << std::endl;
	std::cerr << res << std::endl;

	MatlabGMMDataExchange::GetEngineDenseMatrix("A", mAOfLineSegmentInLvAccelerated);//get the map from matlab
	MatlabGMMDataExchange::GetEngineDenseMatrix("B", mBOfLineSegmentInLvAccelerated);//get the map from matlab

	cout.flush();
	return MS::kSuccess;
}

MStatus SpaceDeformer2D::findLineApproximationForCurve(){
	GMMDenseColMatrix slopeOfLineApproxForCurveInLv(1, 1);
	GMMDenseColMatrix xCoordOfIntersectionPointForCurveLv(1, 1);

	MatlabGMMDataExchange::SetEngineDenseMatrix("k", doubleToGmmMat(this->k));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("SIGMA", doubleToGmmMat(this->SigmaA));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("sigma", doubleToGmmMat(this->sigmaB));//send the matrix to matlab

	std::string res = MatlabInterface::GetEngine().EvalToString("[m,intersectionX]=findLineApproxForCurve(sigma,SIGMA,k);");
	std::cout << res << std::endl;
	std::cerr << res << std::endl;

	MatlabGMMDataExchange::GetEngineDenseMatrix("m", slopeOfLineApproxForCurveInLv);//** not needed **
	MatlabGMMDataExchange::GetEngineDenseMatrix("intersectionX", xCoordOfIntersectionPointForCurveLv);//** not needed **

	mXcoordOfIntersectionPointForCurveLv = xCoordOfIntersectionPointForCurveLv(0, 0);
	mSlopeOfLineApproxForCurveInLv = slopeOfLineApproxForCurveInLv(0, 0);
	mYcoordOfIntersectionPointForCurveLv = log(sigmaB / (1 - mXcoordOfIntersectionPointForCurveLv));

	//use only with min seg
	if (k <= mXcoordOfIntersectionPointForCurveLv){
		mXvaluesOfIntersections = { 0, k, k, 0 };
		mYvaluesOfIntersections = { log(SigmaA), log(SigmaA) - k, log(sigmaB / (1 - k)), log(sigmaB) };
	}
	else{
		mXvaluesOfIntersections = { 0, mXcoordOfIntersectionPointForCurveLv, 0 };
		mYvaluesOfIntersections = { log(SigmaA), log(sigmaB / (1 - mXcoordOfIntersectionPointForCurveLv)), log(sigmaB) };
	}

	//gpu
	DoubleCPUMatrix temp_mXvaluesOfIntersections(mXvaluesOfIntersections.size(),1);
	DoubleCPUMatrix temp_mYvaluesOfIntersections(mYvaluesOfIntersections.size(),1);

	for (int i = 0; i < mXvaluesOfIntersections.size(); i++){
		temp_mXvaluesOfIntersections(i, 0) = mXvaluesOfIntersections[i];
		temp_mYvaluesOfIntersections(i, 0) = mYvaluesOfIntersections[i];

	}
	mXvaluesOfIntersections_gpu = temp_mXvaluesOfIntersections;
	mYvaluesOfIntersections_gpu = temp_mYvaluesOfIntersections;
	//***

	mNeedToCalcNewLine = false;
	cout.flush();
	return MS::kSuccess;
}


