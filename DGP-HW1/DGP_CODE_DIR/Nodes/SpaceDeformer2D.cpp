#include "stdafx.h"

#include "SpaceDeformer2D.h"
#include "Utils/STL_Macros.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"

#define EPS 0.00001

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

	gmm::mult(mCoordinates, mCageVertices, mInternalPoints);//multiply coucy coordinates with the new cage vertices and insert it to the internal points

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


MStatus SpaceDeformer2D::doSetup(MItGeometry& iter, MFnMesh& cageMeshFn)
{
	MStatus stat;

	int m = iter.count(&stat); //num internal points
	int n = mCageVertices.nrows(); //num vertices in the cage

	gmm::clear(mCoordinates);
	gmm::resize(mCoordinates, m, n);

	gmm::clear(mInternalPoints);
	gmm::resize(mInternalPoints, m, 1);

	//offset the cage by epsilon in the direction of the normal such that the triangle mesh (the image) is strictly inside the new cage
	MPointArray cageVertices; //cartesian coordinates
	cageMeshFn.getPoints(cageVertices);
	Complex* cage = new Complex[n];//complex coordinates

	for (int i = 0; i < n; i++)
	{
		int next = i + 1;
		int prev = i - 1;
		if (i == n - 1)
			next = 0;
		if (i == 0)
			prev = n - 1;

		//Creating the close edges vectors
		MVector cageVertexNormal;
		cageVertexNormal = cageVertices[i] - cageVertices[prev] + cageVertices[i] - cageVertices[next];
		cageVertexNormal.normalize();
		cageVertexNormal *= EPS; //epsilon to prevent points to be exactly on cage
		MPoint offset = cageVertices[i] + cageVertexNormal;
		Complex c(offset.x, offset.y);
		cage[i] = c;
	}

	//////////////////////////////////////////////////
	//calculate Cj(z) the Cauchy-Green complex barycentric coordinates
	for (iter.reset(); !iter.isDone(); iter.next())
	{
		int i = iter.index();
		MPoint pt = iter.position();

		Complex z(pt[0], pt[1]); //internal point

		mInternalPoints(i, 0) = z;

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

			Complex Zj = cage[j];
			Complex Zjnext = cage[next];
			Complex Zjprev = cage[prev];

			//the cauchy-green complex barycentric coordinates equesion

			Complex multiplicand1(0, (1 / (-2 * M_PI)));
			Complex multiplicand2 = ((Zjnext - z) / (Zjnext - Zj))*log((Zjnext - z) / (Zj - z)) -
							((Zjprev - z) / (Zj - Zjprev))*log((Zj - z) / (Zjprev - z));
			K = multiplicand1*multiplicand2;

			///////////////////////////////
			mCoordinates(i, j) = K;
		}
	}
	return MS::kSuccess;
}

