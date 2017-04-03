#include "CudaMatrixGPU.h"
#include "CudaMatrixGPUKernel.cu"

#include <cassert>

#define MAX_GRID_DIMENSION (65535)


//parallelization is used only on the rows.
//this function assumes that m is large, otherwise, you get poor utilization of the gpu.
//in order to optimize this function for smaller matrix heights, a different approach should be taken and parallelism should be exploited along the columns as well.
//this can be done by using reduction
bool cuMultMatrixByVectorComplex(int m, int n, int numPitchedRows, const std::complex<float>* d_A, const std::complex<float>* d_X, std::complex<float>* d_b)
{
	if(m <= 0 || n <= 0)
	{
		return false;
	}

	cudaMemset(d_b, 0, sizeof(float2)*m);

	const unsigned int blockSize = 256;

	int numBlocks = m / blockSize + 1;
	int numCycles = n / blockSize + 1;

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	multMatrixByVectorComplexKernel<blockSize><<<grid, threads>>>(m, n, numPitchedRows, numCycles, (const float2*)d_A, (const float2*)d_X, (float2*)d_b);

	return true;
}

bool cuMultMatrixByVectorReal(int m, int n, int numPitchedRows, const float* d_A, const float* d_X, float* d_b)
{
	if(m <= 0 || n <= 0)
	{
		return false;
	}

	cudaMemset(d_b, 0, sizeof(float)*m);

	const unsigned int blockSize = 256;

	int numBlocks = m / blockSize + 1;
	int numCycles = n / blockSize + 1;

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	multMatrixByVectorRealKernel<blockSize><<<grid, threads>>>(m, n, numPitchedRows, numCycles, d_A, d_X, d_b);

	return true;
}

bool cuAddVectorsComplex(int numElements, const std::complex<float>* d_A, const std::complex<float>* d_B, std::complex<float>* d_C)
{
	if(numElements <= 0)
	{
		return false;
	}
	
	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	addVectorsComplexKernel<blockSize><<<grid, threads>>>(numElements, (const float2*)d_A, (const float2*)d_B, (float2*)d_C);

	return true;
}

bool cuAddVectorsReal(int numElements, const float* d_A, const float* d_B, float* d_C)
{
	if(numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	addVectorsRealKernel<blockSize><<<grid, threads>>>(numElements, d_A, d_B, d_C);

	return true;
}

bool cuSubVectorsComplex(int numElements, const std::complex<float>* d_A, const std::complex<float>* d_B, std::complex<float>* d_C)
{
	if(numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	subVectorsComplexKernel<blockSize><<<grid, threads>>>(numElements, (const float2*)d_A, (const float2*)d_B, (float2*)d_C);

	return true;
}

bool cuSubVectorsComplexDouble(int numElements, const std::complex<double>* d_A, const std::complex<double>* d_B, std::complex<double>* d_C)
{
	if (numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	subVectorsComplexDoubleKernel<blockSize> << <grid, threads >> >(numElements, (const double2*)d_A, (const double2*)d_B, (double2*)d_C);

	return true;
}
bool cuSubVectorsReal(int numElements, const float* d_A, const float* d_B, float* d_C)
{
	if(numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;

	const int maxNumElements = blockSize*MAX_GRID_DIMENSION; //the maximum grid size in each dimension is 65535 so if we want to use 1D grid we should take this into account

	int numCycles = I_DIV_UP(numElements, maxNumElements);
	dim3 threads(blockSize, 1);

	for(int i = 0; i < numCycles - 1; i++)
	{
		int numBlocksInCurrentCycle = MAX_GRID_DIMENSION;
		dim3 grid(numBlocksInCurrentCycle, 1);
		subVectorsRealKernel<blockSize><<<grid, threads>>>(maxNumElements, d_A, d_B, d_C);
		d_A +=  (i + 1)*maxNumElements;
		d_B +=  (i + 1)*maxNumElements;
		d_C +=  (i + 1)*maxNumElements;
	}
	int numElementsLeft = numElements - (numCycles - 1)*maxNumElements;
	int numBlocksLeft = I_DIV_UP(numElementsLeft, blockSize);
	dim3 grid(numBlocksLeft, 1);
	subVectorsRealKernel<blockSize><<<grid, threads>>>(numElementsLeft, d_A, d_B, d_C);

	return true;
}

bool cuScaleVectorComplex(int numElements, std::complex<float>* d_A, std::complex<float> alpha)
{
	if(numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	float2 a;
	a.x = alpha.real();
	a.y = alpha.imag();

	scaleVectorComplexKernel<blockSize><<<grid, threads>>>(numElements, (float2*)d_A, a);

	return true;
}

bool cuScaleVectorComplexDouble(int numElements, std::complex<double>* d_A, std::complex<double> alpha)
{
	if (numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	double2 a;
	a.x = alpha.real();
	a.y = alpha.imag();

	scaleVectorComplexDoubleKernel<blockSize> << <grid, threads >> >(numElements, (double2*)d_A, a);

	return true;
}
bool cuExponentVectorComplex(int numElements, std::complex<float>* d_A)
{
	if(numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);

	exponentVectorComplexKernel<blockSize><<<grid, threads>>>(numElements, (float2*)d_A);

	return true;
}

bool cuExponentVectorComplexDouble(int numElements, std::complex<double>* d_A)
{
	if (numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);

	exponentVectorComplexDoubleKernel<blockSize> << <grid, threads >> >(numElements, (double2*)d_A);

	return true;
}

bool cuConvertRealToComplex(int numElements, std::complex<float>* d_A, const float* d_B)
{
	if(numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);

	convertRealToComplexKernel<blockSize><<<grid, threads>>>(numElements, (float2*)d_A, d_B);

	return true;
}

bool cuConjugate(int numElements, std::complex<float>* d_A)
{
	if(numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);

	conjugateComplexKernel<blockSize><<<grid, threads>>>(numElements, (float2*)d_A);

	return true;
}

bool cuAddComplexScalar(int numElements, std::complex<float>* d_A, std::complex<float> alpha)
{
	if(numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;

	const int maxNumElements = blockSize*MAX_GRID_DIMENSION; //the maximum grid size in each dimension is 65535 so if we want to use 1D grid we should take this into account

	int numCycles = I_DIV_UP(numElements, maxNumElements);
	dim3 threads(blockSize, 1);

	float2 a;
	a.x = alpha.real();
	a.y = alpha.imag();

	for(int i = 0; i < numCycles - 1; i++)
	{
		int numBlocksInCurrentCycle = MAX_GRID_DIMENSION;
		dim3 grid(numBlocksInCurrentCycle, 1);
		addComplexScalarKernel<blockSize><<<grid, threads>>>(maxNumElements, (float2*)d_A, a);
		d_A +=  (i + 1)*maxNumElements;
	}
	int numElementsLeft = numElements - (numCycles - 1)*maxNumElements;
	int numBlocksLeft = I_DIV_UP(numElementsLeft, blockSize);
	dim3 grid(numBlocksLeft, 1);
	addComplexScalarKernel<blockSize><<<grid, threads>>>(numElementsLeft, (float2*)d_A, a);

	return true;
}

bool cuAddRealScalar(int numElements, float* d_A, float alpha)
{
	if(numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);

	addRealScalarKernel<blockSize><<<grid, threads>>>(numElements, d_A, alpha);

	return true;
}

