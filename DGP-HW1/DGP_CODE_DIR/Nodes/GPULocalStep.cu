#include "GPULocalStep.h"
#include "GPULocalStepKernel.cu"

#include <cassert>

#define MAX_GRID_DIMENSION (65535)

bool cuProjectPointsToPolygonNoK_split(int numElements, std::complex<double>* log_fz, std::complex<double>* nu_f, const double log_SigmaA, const double sigmaB, const double k, const double xIntersection, const double epsilon, const double m){
	if (numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	projectPointsToPolygonNokKernel << <grid, threads >> >(numElements, (double2*)log_fz, (double2*)nu_f, log_SigmaA, sigmaB, k, xIntersection, epsilon, m);

	return true;

}

bool cuProjectPointsToPolygonNoK_unsplit(int numElements, std::complex<double>* x_vec, const double log_SigmaA, const double sigmaB, const double k, const double xIntersection, const double epsilon, const double m){
	if (numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	projectPointsToPolygonNokKernel_unsplit << <grid, threads >> >(numElements, (double2*)x_vec, log_SigmaA, sigmaB, k, xIntersection, epsilon, m);

	return true;

}

bool cuProjectPointsToPolygonWithK_split(int numElements, std::complex<double>* log_fz, std::complex<double>* nu_f, const double log_SigmaA, const double sigmaB, const double k, const double epsilon, const double m){
	if (numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	projectPointsToPolygonWithkKernel << <grid, threads >> >(numElements, (double2*)log_fz, (double2*)nu_f, log_SigmaA, sigmaB, k, epsilon, m);

	return true;
}

bool cuProjectPointsToPolygonWithK_unsplit(int numElements, std::complex<double>* x_vec, const double log_SigmaA, const double sigmaB, const double k, const double epsilon, const double m){
	if (numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	projectPointsToPolygonWithkKernel_unsplit << <grid, threads >> >(numElements, (double2*)x_vec, log_SigmaA, sigmaB, k, epsilon, m);

	return true;
}

bool cuProjectPointToPolygonMinSeg_split(int numElements, std::complex<double>* log_fz, std::complex<double>* nu_f, const double* mXvaluesOfIntersections, const double* mYvaluesOfIntersections, const int NumOfsegments, const int epsilon){

	if (numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	projectPointToPolygonMinSegKernel <<<grid, threads >>>(numElements, (double2*)log_fz, (double2*)nu_f, mXvaluesOfIntersections, mYvaluesOfIntersections, NumOfsegments, epsilon);

	return true;

}

bool cuProjectPointToPolygonMinSeg_unsplit(int numElements, std::complex<double>* x_vec, const double* mXvaluesOfIntersections, const double* mYvaluesOfIntersections, const int NumOfsegments, const int epsilon){

	if (numElements <= 0)
	{
		return false;
	}

	const unsigned int blockSize = 256;
	int numBlocks = I_DIV_UP(numElements, blockSize);
	assert(numBlocks > 0 && numBlocks <= MAX_GRID_DIMENSION);

	dim3 threads(blockSize, 1);
	dim3 grid(numBlocks, 1);
	projectPointToPolygonMinSegKernel_unsplit << <grid, threads >> >(numElements, (double2*)x_vec, mXvaluesOfIntersections, mYvaluesOfIntersections, NumOfsegments, epsilon);

	return true;

}
