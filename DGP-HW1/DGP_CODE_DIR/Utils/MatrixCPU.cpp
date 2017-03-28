#include "Utils\MatrixCPU.h"
#include <cassert>

/////////////////////////////////////////////////////////////////////////////////////////
/// GLOBAL FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////

bool convertRealFloatToRealDouble(MatrixCPU<double>& realDoubleA, const MatrixCPU<float>& realFloatA)
{
	if(realDoubleA.mNumRows != realFloatA.mNumRows || realDoubleA.mNumColumns != realFloatA.mNumColumns)
	{
		realDoubleA.resize(realFloatA.mNumRows, realFloatA.mNumColumns);
	}

	int numElements = realFloatA.mNumRows*realFloatA.mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		realDoubleA.mData[i] = realFloatA.mData[i];
	}
	return true;
}

bool convertRealDoubleToRealFloat(MatrixCPU<float>& realFloatA, const MatrixCPU<double>& realDoubleA)
{
	if(realDoubleA.mNumRows != realFloatA.mNumRows || realDoubleA.mNumColumns != realFloatA.mNumColumns)
	{
		realFloatA.resize(realDoubleA.mNumRows, realDoubleA.mNumColumns);
	}

	int numElements = realFloatA.mNumRows*realFloatA.mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		realFloatA.mData[i] = (float)realDoubleA.mData[i];
	}
	return true;
}

bool convertComplexDoubleToComplexFloat(MatrixCPU<ComplexFloat>& complexFloatA, const MatrixCPU<ComplexDouble>& complexDoubleA)
{
	if(complexDoubleA.mNumRows != complexFloatA.mNumRows || complexDoubleA.mNumColumns != complexFloatA.mNumColumns)
	{
		complexFloatA.resize(complexDoubleA.mNumRows, complexDoubleA.mNumColumns);
	}

	int numElements = complexFloatA.mNumRows*complexFloatA.mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		complexFloatA.mData[i]._Val[0] = (float)complexDoubleA.mData[i]._Val[0];
		complexFloatA.mData[i]._Val[1] = (float)complexDoubleA.mData[i]._Val[1];
	}
	return true;
}

bool convertComplexFloatToComplexDouble(MatrixCPU<ComplexDouble>& complexDoubleA, const MatrixCPU<ComplexFloat>& complexFloatA)
{
	if(complexDoubleA.mNumRows != complexFloatA.mNumRows || complexDoubleA.mNumColumns != complexFloatA.mNumColumns)
	{
		complexDoubleA.resize(complexFloatA.mNumRows, complexFloatA.mNumColumns);
	}

	int numElements = complexFloatA.mNumRows*complexFloatA.mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		complexDoubleA.mData[i]._Val[0] = complexFloatA.mData[i]._Val[0];
		complexDoubleA.mData[i]._Val[1] = complexFloatA.mData[i]._Val[1];
	}
	return true;
}



