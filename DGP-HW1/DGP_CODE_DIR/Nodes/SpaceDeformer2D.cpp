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





SpaceDeformer2D::SpaceDeformer2D() : mIsFirstTime(true)
{

}

SpaceDeformer2D::~SpaceDeformer2D()
{

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
	CHECK_MSTATUS(attributeAffects(mCoordinateTypeAttr, outputGeom));

	return MStatus::kSuccess;
}

void SpaceDeformer2D::matlabCalcNewVerticesForInterpolation() {
	MatlabGMMDataExchange::SetEngineDenseMatrix("C", mCauchyCoordsOfOriginalCageVertices);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("q", mUserCageVertices);//send the matrix to matlab

//	int res = MatlabInterface::GetEngine().LoadAndRunScript("C:/Users/eden/Documents/MySWProjects/ProjectFinal/DGP-HW1/DGP_CODE_DIR/interpolatedCauchy.m");
	int res = MatlabInterface::GetEngine().LoadAndRunScript("C:/Users/Ben-PC/Documents/MySWprojects/ProjectFinal/DGP-HW1/DGP_CODE_DIR/matlab scripts/interpolatedCauchy.m");
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'interpolatedCauchy.m' failed with error code " << res << std::endl;
	}
	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mInterpolationGenCage_f);//get the incersed matrix from matlab
}

void SpaceDeformer2D::matlabCalcNewVerticesForP2P() {
	MatlabGMMDataExchange::SetEngineDenseMatrix("C", mCauchyCoordsOfOriginalP2P);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("q", mUserP2P);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("D", mSecondDerOfIncCageVertexCoords);//send the matrix to matlab

	int res = MatlabInterface::GetEngine().LoadAndRunScript("C:/Users/Ben-PC/Documents/MySWprojects/ProjectFinal/DGP-HW1/DGP_CODE_DIR/matlab scripts/P2P.m");
	//int res = MatlabInterface::GetEngine().LoadAndRunScript("C:/Users/eden/Documents/MySWProjects/ProjectFinal/DGP-HW1/DGP_CODE_DIR/matlab scripts/P2P.m");
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'interpolatedCauchy.m' failed with error code " << res << std::endl;
	}
	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mP2PGenCageVertices_f);//get the incersed matrix from matlab
}


MStatus SpaceDeformer2D::deform(MDataBlock& block, MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex)
{
	MStatus stat;

	MDataHandle envData = block.inputValue(envelope, &stat);
	if (MS::kSuccess != stat) return stat;

	float env = envData.asFloat();

	MDataHandle coordHandle = block.inputValue(mCoordinateTypeAttr, &stat);
	long coordinateType = coordHandle.asLong();

	MDataHandle handle = block.inputValue(mCageAttr, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	MObject cageMesh = handle.asMesh();

	MFnMesh cageMeshFn(cageMesh, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	MDataHandle p2phandle = block.inputValue(mCageP2pAttr, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	MObject p2pMesh = p2phandle.asMesh();

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

	///// add your code here //////
	///////////////////////////////
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
		matlabCalcNewVerticesForP2P();
		gmm::mult(mCauchyCoordinates, mP2PGenCageVertices_f, mInternalPoints);//get the new internal points 


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
		double numOfNewVeriteciesPerEdge = round((OriginalCageVertecies[(i + 1) % numOfOriginalVertecies].distanceTo(OriginalCageVertecies[i]))/ segmentLength)-1;

		//create a vect the size of a segment in the edge direction
		MVector vec = OriginalCageVertecies[(i + 1) % numOfOriginalVertecies] - (OriginalCageVertecies[i]);
		vec.normalize();
		vec = vec * (vec.length()/ (numOfNewVeriteciesPerEdge+1));

		//create new points
		for (int j = 0; j <= numOfNewVeriteciesPerEdge; j++) {
			IncreasedCageVertecies.append((MPoint)((OriginalCageVertecies[i])+ vec));
		}
	}

}

MStatus SpaceDeformer2D::doSetup(MItGeometry& iter, MFnMesh& cageMeshFn)
{
	MStatus stat;

	int m = iter.count(&stat); //num internal points (point of the triangulated cage)
	int n = mUserCageVertices.nrows(); //num of cage vertices
	int k = mUserP2P.nrows(); //num of control points



	gmm::clear(mCauchyCoordinates);
	gmm::resize(mCauchyCoordinates, m, n);

	gmm::clear(mInternalPoints);
	gmm::resize(mInternalPoints, m, 1);

	gmm::clear(mCauchyCoordsOfOriginalCageVertices);
	gmm::resize(mCauchyCoordsOfOriginalCageVertices, n, n);

	gmm::clear(mInterpolationGenCage_f);
	gmm::resize(mInterpolationGenCage_f, n, 1);

	gmm::clear(mCauchyCoordsOfOriginalP2P);
	gmm::resize(mCauchyCoordsOfOriginalP2P, k, n);

	//offset the cage by epsilon in the direction of the normal such that the triangle mesh (the image) is strictly inside the new cage
	MPointArray cartCageVertices; //cartesian coordinates
	cageMeshFn.getPoints(cartCageVertices);//extract the cage vertecies from the mesh object
	Complex* compCageVertices = new Complex[n];//complex coordinates

	//iterate through the cage vertecies and move all vertecies eps outwards
	for (int i = 0; i < n; i++)
	{
		int next = i + 1;
		int prev = i - 1;
		if (i == n - 1)
			next = 0;
		if (i == 0)
			prev = n - 1;


		MVector cageVertexNormal;//the normal to the vertex

		//Creating the neighbor edges to the cage vertex
		MVector Edge1 = cartCageVertices[i] - cartCageVertices[prev];
		MVector Edge2 = cartCageVertices[i] - cartCageVertices[next];
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
		MVector offset = cartCageVertices[i] + cageVertexNormal;
		compCageVertices[i] = Complex(offset.x, offset.y);
	}

	//calculate the interpolation coordinates to find the new vertices

	populateC(mCauchyCoordsOfOriginalCageVertices, compCageVertices, n, cartCageVertices, n);

#ifndef CVX_INTERPOLATION

	//find the coordinates that will give us the new vertices so that the real cage point will look like interpolation
	//we will do it using cvx
	MatlabGMMDataExchange::SetEngineDenseMatrix("toInverse", mCauchyCoordsOfOriginalCageVertices);//send the matrix to matlab
	//load the matlab script
//	int res = MatlabInterface::GetEngine().LoadAndRunScript("%DGP_CODE_DIR%/matlab scripts/interpolatedCauchy.m"); -----not working
	int res = MatlabInterface::GetEngine().LoadAndRunScript("C:/Users/Ben-PC/Documents/MySWprojects/ProjectFinal/DGP-HW1/DGP_CODE_DIR/matlab scripts/inverse.m");
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'interpolatedCauchy.m' failed with error code " << res << std::endl;
	}
	MatlabGMMDataExchange::GetEngineDenseMatrix("toInverse", mCauchyCoordsOfOriginalCageVertices);//get the incersed matrix from matlab
#endif



	//calculate Cj(z) the Cauchy-Green complex barycentric coordinates
	MPointArray internalPoints;
	for (iter.reset(); !iter.isDone(); iter.next())
	{
		internalPoints.append(iter.position());
	}

	populateC(mCauchyCoordinates, compCageVertices, n, internalPoints, m);
	
	//********************************************************************************
	MPointArray controlPoints;

	for (int i = 0; i < k; i++) {
		Complex c = mUserP2P(i, 0);
		MPoint p(c.real(), c.imag());
		controlPoints.append(p);
	}

	populateC(mCauchyCoordsOfOriginalP2P, compCageVertices, n, controlPoints, k);



	IncreaseVertecies(IN cartCageVertices, OUT mIncreasedVertecies, 500);

	int l = mIncreasedVertecies.length();
	gmm::clear(mSecondDerOfIncCageVertexCoords);
	gmm::resize(mSecondDerOfIncCageVertexCoords, l, n);

	populateD(OUT mSecondDerOfIncCageVertexCoords, compCageVertices, n, mIncreasedVertecies, l);
	//********************************************************************************
	cout.flush();

	return MS::kSuccess;
}



