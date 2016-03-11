#include "stdafx.h"

#include "SpaceDeformer2D.h"
#include "Utils/STL_Macros.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"

#define EPS 0.001

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
		gmm::mult(mCauchyCoordinates, mCageVertices, mInternalPoints);//multiply coucy coordinates with the new cage vertices and insert it to the internal points
		break;
	case 1://Interpolating complex Cauchy coordinates
		gmm::mult(mCauchyCoordsForInterpolation, mCageVertices, mInterpolatedNewVertices_f);
		gmm::mult(mCauchyCoordinates, mInterpolatedNewVertices_f, mInternalPoints);
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


		//Complex c(pt[0], pt[1]);

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
//***********************************************************************************

double max(double x, double y) {
	if (x > y) return x;
	else return y;
}
double min(double x, double y) {
	if (x < y) return x;
	else return y;
}

int orientation(MPoint p0, MPoint p1, MPoint p2)
{
	//calculate cross product of (p1-p0)X (p2-p0)

	int x1 = p1.x - p0.x;
	int y1 = p1.y - p0.y;
	int x2 = p2.x - p0.x;
	int y2 = p2.y - p0.y;

	int determinant = x1*y2 - x2*y1;

	if (determinant == 0)
		return 0;

	if (determinant>0)
		return 1;

	if (determinant<0)
		return 2;

}

int onsegment(MPoint p0, MPoint p1, MPoint p2)
{
	if (p2.x >= min(p0.x, p1.x) && p2.x <= max(p0.x, p1.x))
		return 1;

	return 0;
}

int doIntersect(MPoint p1, MPoint q1, MPoint p2, MPoint q2)
{
	int o1 = orientation(p1, q1, q2);
	int o2 = orientation(p1, q1, p2);
	int o3 = orientation(p2, q2, p1);
	int o4 = orientation(p2, q2, q1);

	if (o1 != o2&&o3 != o4)  //handles general cases
		return 1;

	if (o1 == 0 && o2 == 0 && o3 == 0 && o4 == 0)  //handles special cases when all four points are collinear
	{
		if (onsegment(p1, q1, p2) || onsegment(p1, q1, q2))
			return 1;

	}
	return 0;

}

//***********************************************************************************



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
		cageVertexNormal = Edge1 + Edge2;
		cageVertexNormal.normalize();
		cageVertexNormal *= EPS; //epsilon to prevent points to be exactly on cage

		//move the vertex,the size of epsilon in the direction of the normal
		MPoint offset = cartCageVertices[i] + cageVertexNormal;
		
		//*******fix it****** check if the new edge of the cage intersects the old one, if true - choose (-normal)
	/*	if (i != 0) {
			MPoint cagePrev(compCageVertices[prev].real(), compCageVertices[prev].imag());

			if (doIntersect(cartCageVertices[prev], cartCageVertices[i], cagePrev, offset)) {

				std::cerr << "inside doIntersect in vertex " << i << std::endl;
				std::cerr << "cageVertices[prev] = " << cartCageVertices[prev] << std::endl;
				std::cerr << "cageVertices[i] = " << cartCageVertices[i] << std::endl;
				std::cerr << "cagePrev = " << cagePrev << std::endl;
				std::cerr << "offset = " << offset << std::endl;
				std::cerr << std::endl;

				offset = cartCageVertices[i] - cageVertexNormal;
			}
		}*/
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
	//find the coordinates that will give us the new vertices so that the real cage point will look like interpolation
	//we will do it using cvx
	MatlabGMMDataExchange::SetEngineDenseMatrix("interpolationCoordinates", mCauchyCoordsForInterpolation);//send the matrix to matlab
	//load the matlab script
//	int res = MatlabInterface::GetEngine().LoadAndRunScript("%DGP_CODE_DIR%/matlab scripts/interpolatedCauchy.m");
	int res = MatlabInterface::GetEngine().LoadAndRunScript("C:/Users/Ben-PC/Documents/MySWprojects/ProjectFinal/DGP-HW1/DGP_CODE_DIR/matlab scripts/interpolatedCauchy.m");
	if (res != 0) {//error if failed to load file
		std::cerr << "ERROR: Matlab script 'inverse.m' failed with error code " << res << std::endl;
	}
	MatlabGMMDataExchange::GetEngineDenseMatrix("interpolationCoordinates", mCauchyCoordsForInterpolation);//get the incersed matrix from matlab
	
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

