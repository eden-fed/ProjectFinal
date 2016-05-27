#include "stdafx.h"

#include "SpaceDeformer2D.h"
#include "Utils/STL_Macros.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"

#define EPS 0.001
#define CVX_INTERPOLATION 

#define IN
#define OUT

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




SpaceDeformer2D::SpaceDeformer2D() : mIsFirstTime(true)
{
	mCompCageVertices = NULL;
	mNumOfSegmentsAOld = mNumOfSegmentsA = 0;
	mNLargeOld = mNLarge = 0;
}

SpaceDeformer2D::~SpaceDeformer2D()
{
	if (mCompCageVertices!=NULL)
		delete[] mCompCageVertices;
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
	CHECK_MSTATUS(attributeAffects(mCoordinateTypeAttr, outputGeom));

	MFnNumericAttribute numOfSegmentsAttr;
	mNumOfSegmentsAttr = numOfSegmentsAttr.create("numOfSegmentsA", "numOfSegmentsA", MFnNumericData::kInt, 500, &stat);
	CHECK_MSTATUS(numOfSegmentsAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mNumOfSegmentsAttr));
	CHECK_MSTATUS(attributeAffects(mNumOfSegmentsAttr, outputGeom));

	MFnNumericAttribute nLargeAttr;
	mNlargeAttr = nLargeAttr.create("numOfSegmentsN", "nLarge", MFnNumericData::kInt, 50, &stat);
	CHECK_MSTATUS(nLargeAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mNlargeAttr));
	CHECK_MSTATUS(attributeAffects(mNlargeAttr, outputGeom));

	MFnNumericAttribute kAttr;
	mkAttr = kAttr.create("k", "k", MFnNumericData::kDouble, 0.6, &stat);
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
	return MStatus::kSuccess;
}

void SpaceDeformer2D::matlabCalcNewVerticesForInterpolation() {
	MatlabGMMDataExchange::SetEngineDenseMatrix("C", mCauchyCoordsOfOriginalCageVertices);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("q", mUserCageVertices);//send the matrix to matlab

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

void SpaceDeformer2D::matlabCalcLforHprojection()
{
	MatlabGMMDataExchange::SetEngineDenseMatrix("Ctag", mFirstDerOfIncCageVertexCoords);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("k", doubleToGmmMat(this->k));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("SIGMA", doubleToGmmMat(this->SigmaA));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("sigma", doubleToGmmMat(this->sigmaB));//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("Z0index", doubleToGmmMat(this->mZ0index));//send the matrix to matlab


	int res = MatlabInterface::GetEngine().LoadAndRunScript(RelativeToFullPath("\\matlab scripts\\projectToH.m").c_str());
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'projectToH.m' failed with error code " << res << std::endl;
	}

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

	MDataHandle envData = block.inputValue(envelope, &stat);
	if (MS::kSuccess != stat) return stat;

	env = envData.asFloat();

	MDataHandle coordHandle = block.inputValue(mCoordinateTypeAttr, &stat);
	coordinateType = coordHandle.asLong();

	MDataHandle segAHandle = block.inputValue(mNumOfSegmentsAttr, &stat);
	mNumOfSegmentsA = segAHandle.asInt();//we can only change this value in the first time (doSetup)

	MDataHandle segNHandle = block.inputValue(mNlargeAttr, &stat);
	mNLarge = segNHandle.asInt();//we can only change this value in the first time (doSetup)

	MDataHandle kHandle = block.inputValue(mkAttr, &stat);
	k = kHandle.asDouble();

	MDataHandle sigmaAHandle = block.inputValue(mSigmaaAttr, &stat);
	SigmaA = sigmaAHandle.asDouble();

	MDataHandle sigmaBHandle = block.inputValue(msigmabAttr, &stat);
	sigmaB = sigmaBHandle.asDouble();
	
	MDataHandle z0Handle = block.inputValue(mZ0Attr, &stat);
	float3& z0= z0Handle.asFloat3();
	mz0[0] = z0[0]; mz0[1] = z0[1]; mz0[2] = 0;

	MDataHandle handle = block.inputValue(mCageAttr, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	cageMesh = handle.asMesh();

	MDataHandle p2phandle = block.inputValue(mCageP2pAttr, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	p2pMesh = p2phandle.asMesh();

	return stat;
}


MStatus SpaceDeformer2D::deform(MDataBlock& block, MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex)
{
		
	MStatus stat;
	MObject cageMesh;
	MObject p2pMesh;
	stat = getData(block, cageMesh, p2pMesh);

	MCHECKERROR(stat, "couls not get the data block");

	MFnMesh cageMeshFn(cageMesh, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	MFnMesh p2pMeshFn(p2pMesh, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	
	updateCage(cageMeshFn);//populates mUserCageVertices
	updateControlPoints(p2pMeshFn);//populates mUserP2P

	if (mIsFirstTime)
	{
		stat = doSetup(iter, cageMeshFn);
		CHECK_MSTATUS_AND_RETURN_IT(stat);
		mIsFirstTime = false;
	}

	//compute the deformation of all the internal points. This is done by simply multiplying the coordinate matrix by the cage vertices vector
	MatlabInterface::GetEngine().Eval("clear");
	switch (coordinateType)
	{
	case 0://Approximating cage-based complex Cauchy coordinates
		gmm::mult(mCauchyCoordinates, mUserCageVertices, mInternalPoints);//multiply cauchy coordinates with the new cage vertices and insert it to the internal points
		break;
	case 1://Interpolating complex Cauchy coordinates
#ifdef CVX_INTERPOLATION //calculate using cvx
		matlabCalcNewVerticesForInterpolation();
#else//calculate using the inverse matrix
		gmm::mult(mCauchyCoordsOfOriginalCageVertices, mUserCageVertices, mInterpolationGenCage_f);
#endif
		gmm::mult(mCauchyCoordinates, mInterpolationGenCage_f, mInternalPoints);//get the new internal points 
		break;
	case 2://Point2Point coordinates
		//runtime
		runTimeDoSetup();
		matlabCalcNewVerticesForP2P();
		gmm::mult(mCauchyCoordinatesIncForP2P, mP2PGenCageVertices_f, mInternalPoints);//get the new internal points 
		break;
	case 3:
		//runtime
		runTimeDoSetup();
		matlabCalcLforHprojection();
		break;
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

	gmm::clear(mUserCageVertices);
	gmm::resize(mUserCageVertices, numV, 1);

	MPointArray vertexArray;
	stat = cageMeshFn.getPoints(vertexArray);

	assert(numV == vertexArray.length());

	for (int i = 0; i < numV; i++)
	{
		MPoint p = vertexArray[i];
		Complex c(p[0], p[1]);
		mUserCageVertices(i, 0) = c;
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

void IncreaseVertecies(MPointArray& OriginalCageVertecies, MPointArray& IncreasedCageVertecies,int numOfIncreasedCageVertecies) {
	//find the circumference of the cgae polygon
	double circumference=0;
	int numOfOriginalVertecies = OriginalCageVertecies.length();

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
		double numOfSegmentsPerEdge = round(edgeLength / segmentLength);//the number of vertices in the edge is #seg-1

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
	}

}


void IncreaseVertecies(Complex* OriginalCompCageVertecies, int OrigCageSize, Complex** IncreasedCompCageVertecies, int& numOfIncreasedCageVertecies) {
	MPointArray OriginalCartCageVertecies;
	MPointArray IncreasedCageVertecies;

	for (int i = 0; i < OrigCageSize; i++) {
		Complex c = OriginalCompCageVertecies[i];
		MPoint p(c.real(), c.imag());
		OriginalCartCageVertecies.append(p);
	}

	IncreaseVertecies(OriginalCartCageVertecies, IncreasedCageVertecies, numOfIncreasedCageVertecies);

	int length1 = OriginalCartCageVertecies.length();
	numOfIncreasedCageVertecies = IncreasedCageVertecies.length();

	*IncreasedCompCageVertecies = new Complex[numOfIncreasedCageVertecies];

	for (int i = 0; i < numOfIncreasedCageVertecies; i++) {
		Complex c(IncreasedCageVertecies[i].x, IncreasedCageVertecies[i].y);
		(*IncreasedCompCageVertecies)[i] = c;
	}
}

MStatus SpaceDeformer2D::doSetup(MItGeometry& iter, MFnMesh& cageMeshFn)
{
	MStatus stat;

	mNumOfInternalPoints = iter.count(&stat); //num internal points (point of the triangulated cage)
	mNumOfCageVerticies = mUserCageVertices.nrows(); //num of cage vertices
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
	cageMeshFn.getPoints(mCartCageVertices);//extract the cage vertecies from the mesh object
	this->mCompCageVertices = new Complex[mNumOfCageVerticies];//complex coordinates

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
		MVector Edge1 = mCartCageVertices[i] - mCartCageVertices[prev];
		MVector Edge2 = mCartCageVertices[i] - mCartCageVertices[next];
		Edge1.normalize();
		Edge2.normalize();

		//the notmal to the vertes is the angle bisector of the two edges
		float V = (Edge1^Edge2).z;
		if (V>=0) {//turn tight
			cageVertexNormal = -(Edge1 + Edge2);
		}else {//turn left
			cageVertexNormal = Edge1 + Edge2;
		}
		
		cageVertexNormal.normalize();
		cageVertexNormal *= EPS; //epsilon to prevent points to be exactly on cage

		//move the vertex,the size of epsilon in the direction of the normal
		MVector offset = mCartCageVertices[i] + cageVertexNormal;
		mCompCageVertices[i] = Complex(offset.x, offset.y);
	}

	//calculate the interpolation coordinates to find the new vertices
	
	populateC(mCauchyCoordsOfOriginalCageVertices, mCompCageVertices, mNumOfCageVerticies, mCartCageVertices, mNumOfCageVerticies);

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
		this->mInternalPoints_MPoint.append(iter.position());
	}

	populateC(mCauchyCoordinates, mCompCageVertices, mNumOfCageVerticies, this->mInternalPoints_MPoint, mNumOfInternalPoints);


	
	//runTimeDoSetup();
	for (int i = 0; i < mNumOfControlPoints; i++) {
		Complex c = mUserP2P(i, 0);
		MPoint p(c.real(), c.imag());
		mInitialcontrolPoints.append(p);
	}
	return MS::kSuccess;
}



int SpaceDeformer2D::findClosestInternalPointsToZ0() {
	MPoint Z0(mz0[0], mz0[1], mz0[2],0);
	double minDist = DBL_MAX;
	int index=0;

	for (int i = 0; i < mInternalPoints_MPoint.length(); i++ ) {
		MPoint tempPoint = mInternalPoints_MPoint[i];
		double tempDist = sqrt((Z0.x - tempPoint.x) *(Z0.x - tempPoint.x) + (Z0.y - tempPoint.y) *(Z0.y - tempPoint.y));

		if (tempDist < minDist) {
			minDist = tempDist;
			index = i;
		}
	}
	return index;
}


MStatus SpaceDeformer2D::runTimeDoSetup() {
	if (mNumOfSegmentsAOld == mNumOfSegmentsA && mNLargeOld == mNLarge) {
		return MS::kSuccess;
	}
	mNumOfSegmentsAOld = mNumOfSegmentsA;
	mNLargeOld = mNLarge;

	mZ0index = findClosestInternalPointsToZ0();

	Complex* IncreasedCompCageVertecies = NULL;// = new Complex[nLarge];

	IncreaseVertecies(IN mCompCageVertices, IN mNumOfCageVerticies, OUT &IncreasedCompCageVertecies, mNLarge);//nLarge might change in this func

	gmm::clear(mCauchyCoordinatesIncForP2P);
	gmm::resize(mCauchyCoordinatesIncForP2P, mNumOfInternalPoints, mNLarge);


	populateC(OUT mCauchyCoordinatesIncForP2P, IncreasedCompCageVertecies, mNLarge, this->mInternalPoints_MPoint, mNumOfInternalPoints);

	//*******************************************************************************
	//calculate D & C matrices for P2P

	gmm::clear(mCauchyCoordsOfOriginalP2P);
	gmm::resize(mCauchyCoordsOfOriginalP2P, mNumOfControlPoints, mNLarge);
	populateC(mCauchyCoordsOfOriginalP2P, IncreasedCompCageVertecies, mNLarge, mInitialcontrolPoints, mNumOfControlPoints);

	MPointArray increasedVertecies_a; //dimensions are: a x 1
	
	IncreaseVertecies(IN mCartCageVertices, OUT increasedVertecies_a, mNumOfSegmentsA);

	int a = increasedVertecies_a.length();
	gmm::clear(mSecondDerOfIncCageVertexCoords);
	gmm::resize(mSecondDerOfIncCageVertexCoords, a, mNLarge);

	populateD(OUT mSecondDerOfIncCageVertexCoords, IncreasedCompCageVertecies, mNLarge, increasedVertecies_a, a);

	//********************************************************************************
	//calculate the first derivative of the cauchy grenn coords for the projection to H 
	//the matrix Ctag need to be a x nLarge

	gmm::clear(mFirstDerOfIncCageVertexCoords);
	gmm::resize(mFirstDerOfIncCageVertexCoords, a, mNLarge);

	populateCtag(OUT mFirstDerOfIncCageVertexCoords, IncreasedCompCageVertecies, mNLarge, increasedVertecies_a, a);

	//********************************************************************************
	cout.flush();
	delete[] IncreasedCompCageVertecies;
	return MS::kSuccess;
}