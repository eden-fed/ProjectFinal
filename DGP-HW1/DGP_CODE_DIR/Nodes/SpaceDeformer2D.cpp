#include "stdafx.h"

#include "SpaceDeformer2D.h"
#include "Utils/STL_Macros.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"

#define EPS 0.001
#define CVX_INTERPOLATION 

const MTypeId SpaceDeformer2D::mTypeId(0x6723c);
const MString SpaceDeformer2D::mTypeName("SpaceDeformer2D");

MObject SpaceDeformer2D::mCageAttr;
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
	MatlabGMMDataExchange::SetEngineDenseMatrix("A", mCauchyCoordsForInterpolation);//send the matrix to matlab
	MatlabGMMDataExchange::SetEngineDenseMatrix("q", mCageVertices);//send the matrix to matlab

	int res = MatlabInterface::GetEngine().LoadAndRunScript("C:/Users/Ben-PC/Documents/MySWprojects/ProjectFinal/DGP-HW1/DGP_CODE_DIR/matlab scripts/interpolatedCauchy.m");
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'interpolatedCauchy.m' failed with error code " << res << std::endl;
	}
	MatlabGMMDataExchange::GetEngineDenseMatrix("f", mInterpolatedNewVertices_f);//get the incersed matrix from matlab
}

MStatus SpaceDeformer2D::deform(MDataBlock& block, MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex)
{
	MStatus stat;

	MDataHandle envData = block.inputValue(envelope, &stat);
	if (MS::kSuccess != stat) return stat;

	float env = envData.asFloat();

	MDataHandle coordHandle = block.inputValue(mCoordinateTypeAttr, &stat);
	short coordinateType = coordHandle.asShort();

	MDataHandle handle = block.inputValue(mCageAttr, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	MObject cageMesh = handle.asMesh();

	MFnMesh cageMeshFn(cageMesh, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	updateCage(cageMeshFn);

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
		gmm::mult(mCauchyCoordinates, mCageVertices, mInternalPoints);//multiply cauchy coordinates with the new cage vertices and insert it to the internal points
		break;
	case 1://Interpolating complex Cauchy coordinates
		#ifdef CVX_INTERPOLATION //calculate using cvx
		matlabCalcNewVerticesForInterpolation();
		#else//calculate using the inverse matrix
		gmm::mult(mCauchyCoordsForInterpolation, mCageVertices, mInterpolatedNewVertices_f);
		#endif
		gmm::mult(mCauchyCoordinates, mInterpolatedNewVertices_f, mInternalPoints);//get the new internal points 
		break;
	case 2://Point2Point coordinates
		MGlobal::displayError("Point 2 Point deformation hasnt been implemented yet .");
		return MS::kFailure;
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

	gmm::clear(mCageVertices);
	gmm::resize(mCageVertices, numV, 1);

	MPointArray vertexArray;
	stat = cageMeshFn.getPoints(vertexArray);

	assert(numV == vertexArray.length());

	for (int i = 0; i < numV; i++)
	{
		MPoint p = vertexArray[i];
		Complex c(p[0], p[1]);
		mCageVertices(i, 0) = c;
	}
	return MS::kSuccess;
}


MStatus SpaceDeformer2D::doSetup(MItGeometry& iter, MFnMesh& cageMeshFn)
{
	MStatus stat;

	int m = iter.count(&stat); //num internal points (point of the triangulated cage)
	int n = mCageVertices.nrows(); //num of cage vertices

	gmm::clear(mCauchyCoordinates);
	gmm::resize(mCauchyCoordinates, m, n);

	gmm::clear(mInternalPoints);
	gmm::resize(mInternalPoints, m, 1);

	gmm::clear(mCauchyCoordsForInterpolation);
	gmm::resize(mCauchyCoordsForInterpolation, n, n);

	gmm::clear(mInterpolatedNewVertices_f);
	gmm::resize(mInterpolatedNewVertices_f, n, 1);

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
	for (int i = 0; i < n; i++) {
		MPoint pt = cartCageVertices[i];
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

			Complex Zj = compCageVertices[j];
			Complex Zjnext = compCageVertices[next];
			Complex Zjprev = compCageVertices[prev];

			//the cauchy-green complex barycentric coordinates equesion

			Complex multiplicand1(0, (1 / (-2 * M_PI)));
			Complex multiplicand2 = ((Zjnext - z) / (Zjnext - Zj))*log((Zjnext - z) / (Zj - z)) -
				((Zjprev - z) / (Zj - Zjprev))*log((Zj - z) / (Zjprev - z));
			K = multiplicand1*multiplicand2;

			mCauchyCoordsForInterpolation(i, j) = K;
		}
	}

	#ifndef CVX_INTERPOLATION

	//find the coordinates that will give us the new vertices so that the real cage point will look like interpolation
	//we will do it using cvx
	MatlabGMMDataExchange::SetEngineDenseMatrix("toInverse", mCauchyCoordsForInterpolation);//send the matrix to matlab
	//load the matlab script
//	int res = MatlabInterface::GetEngine().LoadAndRunScript("%DGP_CODE_DIR%/matlab scripts/interpolatedCauchy.m"); -----not working
	int res = MatlabInterface::GetEngine().LoadAndRunScript("C:/Users/Ben-PC/Documents/MySWprojects/ProjectFinal/DGP-HW1/DGP_CODE_DIR/matlab scripts/inverse.m");
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'interpolatedCauchy.m' failed with error code " << res << std::endl;
	}
	MatlabGMMDataExchange::GetEngineDenseMatrix("toInverse", mCauchyCoordsForInterpolation);//get the incersed matrix from matlab
	#endif
	//////////////////////////////////////////////////
	//calculate Cj(z) the Cauchy-Green complex barycentric coordinates
	for (iter.reset(); !iter.isDone(); iter.next())
	{
		int i = iter.index();
		MPoint pt = iter.position();

		Complex z(pt[0], pt[1]); //internal point

		mInternalPoints(i, 0) = z;

		//for each  internal point i, we would like to calculat the n Cauchy coordinates (1,2,.....,n)
		for (int j = 0; j < n; j++)
		{
			Complex K(0.0, 0.0);

			///// add your code here //////
			///////////////////////////////
			//update K to be value of the j'th coordinate at the i'th internal point

			int next = j + 1;
			int prev = j - 1;
			if (j == n - 1)
				next = 0;
			if (j == 0)
				prev = n - 1;

			Complex Zj = compCageVertices[j];
			Complex Zjnext = compCageVertices[next];
			Complex Zjprev = compCageVertices[prev];

			//the cauchy-green complex barycentric coordinates equesion

			Complex multiplicand1(0, (1 / (-2 * M_PI)));
			Complex multiplicand2 = ((Zjnext - z) / (Zjnext - Zj))*log((Zjnext - z) / (Zj - z)) -
				((Zjprev - z) / (Zj - Zjprev))*log((Zj - z) / (Zjprev - z));
			K = multiplicand1*multiplicand2;

			///////////////////////////////
			mCauchyCoordinates(i, j) = K;
		}
	}
	cout.flush();

	return MS::kSuccess;
}

