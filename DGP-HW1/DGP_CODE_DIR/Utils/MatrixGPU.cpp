#include "Utils\MatrixGPU.h"

/////////////////////////////////////////////////////////////////////////////////////////
/// GLOBAL FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////

bool convertRealToComplex(MatrixGPU<ComplexFloat>& complexA, const MatrixGPU<float>& realA, double* gpuTime)
{
	if(complexA.mNumRows != realA.mNumRows || complexA.mNumColumns != realA.mNumColumns || complexA.mNumPitchedRows != realA.mNumPitchedRows)
	{
		assert(0);
		return false;
	}

	CUDATimer timer;
	if(gpuTime)
	{
		*gpuTime = 0.0;
		timer.startTimer();
	}

	cuConvertRealToComplex(complexA.mNumPitchedRows*complexA.mNumColumns, complexA.mData, realA.mData);

	if(gpuTime)
	{
		*gpuTime = timer.stopTimer();
	}

	return true;
}
