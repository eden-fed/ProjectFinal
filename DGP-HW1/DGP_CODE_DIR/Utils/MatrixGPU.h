#pragma once

////////////////////////////////////////////
/////////        MatrixGPU        //////////
///////// Copyright (C) Ofir Weber /////////
////////////////////////////////////////////


#include <cublas.h>
//#include <cutil.h>
#include <cuda_runtime_api.h>
#include <builtin_types.h>
//#include <cutil_inline_runtime.h>

#include <Utils\CudaMatrixGPU.h>
#include <Utils\Utilities.h>
#include <Utils\Timer.h>

#include <complex>
#include <vector>
#include <string>
#include <cassert>

#include "Utils\FastMatrix.h"



bool convertRealToComplex(MatrixGPU<ComplexFloat>& complexA, const MatrixGPU<float>& realA, double* gpuTime = NULL);

template <class T>
class MatrixGPU
{
	friend class MatrixCPU<T>;
	friend bool convertRealToComplex(MatrixGPU<ComplexFloat>& A, const MatrixGPU<float>& B, double* gpuTime);

public:

	MatrixGPU();
	MatrixGPU(int m, int n);
	MatrixGPU(const MatrixGPU <T>& mat);
	MatrixGPU(const MatrixCPU <T>& mat);
	MatrixGPU(const T* dataCPU, int m, int n);
	const MatrixGPU<T>& operator=(const MatrixGPU<T>& mat);
	const MatrixGPU<T>& operator=(const MatrixCPU<T>& mat);
	~MatrixGPU();
	bool concatenateColVec(const MatrixGPU<T>& mat1, const MatrixGPU<T>& mat2);//***
	bool splitInMiddleColVec(MatrixGPU<T>& mat1, MatrixGPU<T>& mat2);//***
	inline const T get_slow(int i, int j) const;
	void set_slow(int i, int j, const T& value);
	inline bool isValid() const;
	bool resize(int m, int n);
	void freeData();
	int nCols() const;
	int nRows() const;
	int nPitchedRows() const;
	const T* getData();
	bool mult(const MatrixGPU<T>& A, const MatrixGPU<T>& B, char opA = 'N', char opB = 'N', double* gpuTime = NULL);
	ComplexDouble dotColVec(const MatrixGPU<T>& v, double* gpuTime = NULL);//***
	double norm(double* gpuTime = NULL);//***
	bool add(const MatrixGPU<T>& A, const MatrixGPU<T>& B, double* gpuTime = NULL);
	bool add(T alpha);
	bool sub(const MatrixGPU<T>& A, const MatrixGPU<T>& B, double* gpuTime = NULL);
	bool sub(T alpha);
	bool scale(T alpha, double* gpuTime = NULL);
	bool exponent(double* gpuTime = NULL);
	bool conjugate();
	bool exportToMatlabFile(const std::string& matrixName = "") const;


protected:

	bool cuda_add(const MatrixGPU<float>& A, const MatrixGPU<float>& B);
	bool cuda_add(const MatrixGPU<ComplexFloat>& A, const MatrixGPU<ComplexFloat>& B);
	bool cuda_add(MatrixGPU<ComplexFloat>& A, ComplexFloat alpha) const;
	bool cuda_add(MatrixGPU<float>& A, float alpha) const;
	bool cuda_sub(const MatrixGPU<float>& A, const MatrixGPU<float>& B);
	bool cuda_sub(const MatrixGPU<ComplexFloat>& A, const MatrixGPU<ComplexFloat>& B);
	bool cuda_sub(const MatrixGPU<ComplexDouble>& A, const MatrixGPU<ComplexDouble>& B);
	bool cuda_conjugate(MatrixGPU<ComplexFloat>& A) const;
	bool cuda_conjugate(MatrixGPU<float>& A) const;
	int roundToCoalesced(int numElements);
	void cublas_gemm(const MatrixGPU<ComplexFloat>& A, const MatrixGPU<ComplexFloat>& B, int m, int n, int k, char opA, char opB);
	void cublas_gemm(const MatrixGPU<ComplexDouble>& A, const MatrixGPU<ComplexDouble>& B, int m, int n, int k, char opA, char opB);
	void cublas_gemm(const MatrixGPU<float>& A, const MatrixGPU<float>& B, int m, int n, int k, char opA, char opB);
	void cublas_gemm(const MatrixGPU<double>& A, const MatrixGPU<double>& B, int m, int n, int k, char opA, char opB);
	void cublas_gemv(const MatrixGPU<float>& A, const MatrixGPU<float>& v, int m, int n, char opA);
	void cublas_gemv(const MatrixGPU<double>& A, const MatrixGPU<double>& v, int m, int n, char opA);
	void cublas_gemv(const MatrixGPU<ComplexFloat>& A, const MatrixGPU<ComplexFloat>& v, int m, int n, char opA);
	void cublas_gemv(const MatrixGPU<ComplexDouble>& A, const MatrixGPU<ComplexDouble>& v, int m, int n, char opA);
	ComplexDouble cublas_dot(const MatrixGPU<ComplexDouble>& v, int n);
	double cublas_norm(int n);

protected:

	//data is stored in column major manner (Fortran style) in a coalesced manner
	T* mData;
	int mNumRows;
	int mNumColumns;
	int mNumPitchedRows;
};




template <class T>
MatrixGPU<T>::MatrixGPU()
{
	mData = NULL;
	mNumRows = 0;
	mNumColumns = 0;
	mNumPitchedRows = 0;
}

template <class T>
MatrixGPU<T>::MatrixGPU(int m, int n)
{
	mData = NULL;
	mNumRows = 0;
	mNumColumns = 0;
	mNumPitchedRows = 0;

	if(m <= 0 || n <= 0)
	{
		assert(0);
		return;
	}

	int pitchRows = roundToCoalesced(m);

	cudaError_t stat = cudaMalloc((void**)&mData, pitchRows*sizeof(T)*n);

	if(mData == NULL || stat != cudaSuccess)
	{
		assert(0);
		return;
	}
	else
	{
		mNumRows = m;
		mNumColumns = n;
		mNumPitchedRows = pitchRows;
	}

	stat = cudaMemset((void*)mData, 0, pitchRows*sizeof(T)*n);

	if(stat != cudaSuccess)
	{
		assert(0);
		return;
	}
}

//copy c'tor GPU to GPU
template <class T>
MatrixGPU<T>::MatrixGPU(const MatrixGPU<T>& mat)
{
	mData = NULL;
	mNumRows = 0;
	mNumColumns = 0;
	mNumPitchedRows = 0;

	if(mat.mNumRows == 0 && mat.mNumColumns == 0 && mat.mNumPitchedRows == 0) //an empty matrix
	{
		return;
	}

	if(mat.mNumRows <= 0 || mat.mNumColumns <= 0 || mat.mNumPitchedRows <= 0)
	{
		assert(0);
		return;
	}

	cudaError_t stat = cudaMalloc((void**)&mData, mat.mNumPitchedRows*mat.mNumColumns*sizeof(T));

	if(mData == NULL || stat != cudaSuccess)
	{
		assert(0);
		return;
	}
	mNumRows = mat.mNumRows;
	mNumColumns = mat.mNumColumns;
	mNumPitchedRows = mat.mNumPitchedRows;

	checkCudaErrors(cudaMemcpy(mData, mat.mData, mat.mNumPitchedRows*mat.mNumColumns*sizeof(T), cudaMemcpyDeviceToDevice));
}

//copy c'tor CPU to GPU
template <class T>
MatrixGPU<T>::MatrixGPU(const MatrixCPU<T>& mat)
{
	mData = NULL;
	mNumRows = 0;
	mNumColumns = 0;
	mNumPitchedRows = 0;

	if(mat.mNumRows == 0 && mat.mNumColumns == 0) //an empty matrix
	{
		return;
	}

	if(mat.mNumRows <= 0 || mat.mNumColumns <= 0)
	{
		assert(0);
		return;
	}

	int pitchRows = roundToCoalesced(mat.mNumRows);

	cudaError_t stat = cudaMalloc((void**)&mData, pitchRows*sizeof(T)*mat.mNumColumns);

	if(mData == NULL || stat != cudaSuccess)
	{
		assert(0);
		return;
	}
	mNumRows = mat.mNumRows;
	mNumColumns = mat.mNumColumns;
	mNumPitchedRows = pitchRows;

	checkCudaErrors(safeCudaMemcpy2D(mData, pitchRows*sizeof(T), mat.mData, mat.mNumRows*sizeof(T), mat.mNumRows*sizeof(T), mat.mNumColumns, cudaMemcpyHostToDevice));
}

//c'tor that receive raw CPU data
template <class T>
MatrixGPU<T>::MatrixGPU(const T* dataCPU, int m, int n)
{
	mData = NULL;
	mNumRows = 0;
	mNumColumns = 0;
	mNumPitchedRows = 0;

	if(m <= 0 || n <= 0)
	{
		assert(0);
		return;
	}

	int pitchRows = roundToCoalesced(m);

	cudaError_t stat = cudaMalloc((void**)&mData, pitchRows*sizeof(T)*n);

	if(mData == NULL || stat != cudaSuccess)
	{
		assert(0);
		return;
	}
	mNumRows = m;
	mNumColumns = n;
	mNumPitchedRows = pitchRows;

	checkCudaErrors(safeCudaMemcpy2D(mData, pitchRows*sizeof(T), dataCPU, mat.mNumRows*sizeof(T), mat.mNumRows*sizeof(T), mat.mNumColumns, cudaMemcpyHostToDevice));
}

//assignment operator GPU to GPU
template <class T>
const MatrixGPU<T>& MatrixGPU<T>::operator=(const MatrixGPU<T>& mat)
{
	if(mat.mNumRows <= 0 || mat.mNumColumns <= 0 || mat.mNumPitchedRows <= 0)
	{
		assert(0);
		return *this;
	}
	if(mat.mNumPitchedRows*mat.mNumColumns != mNumPitchedRows*mNumColumns)
	{
		if(mData)
		{
			cudaError_t stat = cudaFree(mData);
			assert(stat == cudaSuccess);
			mData = NULL;
		}
		mNumRows = 0;
		mNumColumns = 0;
		mNumPitchedRows = 0;

		cudaError_t stat = cudaMalloc((void**)&mData, mat.mNumPitchedRows*mat.mNumColumns*sizeof(T));

		if(mData == NULL || stat != cudaSuccess)
		{
			assert(0);
			return *this;
		}
	}
	mNumRows = mat.mNumRows;
	mNumColumns = mat.mNumColumns;
	mNumPitchedRows = mat.mNumPitchedRows;

	checkCudaErrors(cudaMemcpy(mData, mat.mData, mat.mNumPitchedRows*mat.mNumColumns*sizeof(T), cudaMemcpyDeviceToDevice));

	return *this;
}

//assignment operator CPU to GPU
template <class T>
const MatrixGPU<T>& MatrixGPU<T>::operator=(const MatrixCPU<T>& mat)
{
	if(mat.mNumRows <= 0 || mat.mNumColumns <= 0)
	{
		assert(0);
		return *this;
	}
	int pitchRows = roundToCoalesced(mat.mNumRows);
	if(pitchRows*mat.mNumColumns != mNumPitchedRows*mNumColumns)
	{
		if(mData)
		{
			cudaError_t stat = cudaFree(mData);
			assert(stat == cudaSuccess);
			mData = NULL;
		}
		mNumRows = 0;
		mNumColumns = 0;

		cudaError_t stat = cudaMalloc((void**)&mData, pitchRows*sizeof(T)*mat.mNumColumns);

		if(mData == NULL || stat != cudaSuccess)
		{
			assert(0);
			return *this;
		}
	}
	mNumRows = mat.mNumRows;
	mNumColumns = mat.mNumColumns;
	mNumPitchedRows = pitchRows;

	size_t sourcePitchInBytes = mat.mNumRows*sizeof(T);
	size_t destinationPitchInBytes = pitchRows*sizeof(T);
	size_t dataInBytes = mat.mNumRows*sizeof(T);

	checkCudaErrors(safeCudaMemcpy2D(mData, destinationPitchInBytes, mat.mData, sourcePitchInBytes, dataInBytes, mat.mNumColumns, cudaMemcpyHostToDevice));

	return *this;
}

/*concatenate column vectors vertically*/
template <class T>
bool MatrixGPU<T>::concatenateColVec(const MatrixGPU<T>& mat1, const MatrixGPU<T>& mat2){
	if (mat1.mNumRows == 0 && mat1.mNumColumns == 0 && mat1.mNumPitchedRows == 0) //an empty matrix
	{
		*this = mat2;
		return true;
	}
	if (mat2.mNumRows == 0 && mat2.mNumColumns == 0 && mat2.mNumPitchedRows == 0) //an empty matrix
	{
		*this = mat1;
		return true;
	}
	if (mat1.mNumRows <= 0 || mat1.mNumColumns != 1 || mat1.mNumPitchedRows <= 0 || mat2.mNumRows <= 0 || mat2.mNumColumns != 1 || mat2.mNumPitchedRows <= 0)
	{
		assert(0);
		return false;
	}

	int pitchRows = roundToCoalesced(mat1.mNumRows + mat2.mNumRows);

	if (pitchRows != mNumPitchedRows*mNumColumns)
	{
		if (mData)
		{
			cudaError_t stat = cudaFree(mData);
			assert(stat == cudaSuccess);
			mData = NULL;
		}
		mNumRows = 0;
		mNumColumns = 0;
		mNumPitchedRows = 0;

		cudaError_t stat = cudaMalloc((void**)&mData, (pitchRows)*mat1.mNumColumns*sizeof(T));

		if (mData == NULL || stat != cudaSuccess)
		{
			assert(0);
			return false;
		}
	}
	mNumRows = mat1.mNumRows + mat2.mNumRows;
	mNumColumns = 1;
	mNumPitchedRows = pitchRows;

	checkCudaErrors(cudaMemcpy(mData, mat1.mData, mat1.mNumPitchedRows*sizeof(T), cudaMemcpyDeviceToDevice));
	checkCudaErrors(cudaMemcpy(mData + mat1.mNumRows, mat2.mData, mat2.mNumRows*sizeof(T), cudaMemcpyDeviceToDevice));

	return true;

}

/*split vertically column vector to 2 vectors in the same size*/
template <class T>
bool MatrixGPU<T>::splitInMiddleColVec(MatrixGPU<T>& mat1, MatrixGPU<T>& mat2){

	if (mNumRows <= 0 || mNumColumns != 1 || mNumPitchedRows <= 0)
	{
		assert(0);
		return false;
	}
	if (mNumRows % 2 != 0 || mNumPitchedRows % 2 != 0)
	{
		assert(0);
		return false;
	}

	int pitchRows = roundToCoalesced(mNumRows / 2);
	if ((mat1.mNumPitchedRows)*mat1.mNumColumns != pitchRows)
	{

		mat1.freeData();
		cudaError_t stat = cudaMalloc((void**)&mat1.mData, pitchRows*mNumColumns*sizeof(T));

		if (mat1.mData == NULL || stat != cudaSuccess)
		{
			assert(0);
			return false;
		}
	}
	mat1.mNumRows = mNumRows/2;
	mat1.mNumColumns = 1;
	mat1.mNumPitchedRows = pitchRows;

	if ((mat2.mNumPitchedRows)*mat2.mNumColumns != pitchRows)
	{

		mat2.freeData();
		cudaError_t stat = cudaMalloc((void**)&mat2.mData, pitchRows*mNumColumns*sizeof(T));

		if (mat2.mData == NULL || stat != cudaSuccess)
		{
			assert(0);
			return false;
		}
	}
	mat2.mNumRows = mNumRows / 2;
	mat2.mNumColumns = 1;
	mat2.mNumPitchedRows = pitchRows;

	checkCudaErrors(cudaMemcpy(mat1.mData, mData, mat1.mNumRows*sizeof(T), cudaMemcpyDeviceToDevice));
	checkCudaErrors(cudaMemcpy(mat2.mData, mData + mat1.mNumRows, mat2.mNumRows*sizeof(T), cudaMemcpyDeviceToDevice));

	return true;
}

template <class T>
MatrixGPU<T>::~MatrixGPU()
{
	if(mData)
	{
		cudaError_t stat = cudaFree(mData);
		assert(stat == cudaSuccess);
		mData = NULL;
	}
}

template <class T>
void MatrixGPU<T>::freeData()
{
	if(mData)
	{
		cudaError_t stat = cudaFree(mData);
		assert(stat == cudaSuccess);
		mData = NULL;
	}
	mNumRows = 0;
	mNumColumns = 0;
	mNumPitchedRows = 0;
}

template <class T>
bool MatrixGPU<T>::resize(int m, int n)
{
	if(m <= 0 || n <= 0)
	{
		assert(0);
		return false;
	}

	if(mData)
	{
		cudaError_t stat = cudaFree(mData);
		assert(stat == cudaSuccess);
		mData = NULL;
		mNumRows = 0;
		mNumColumns = 0;
		mNumPitchedRows = 0;
	}

	int pitchRows = roundToCoalesced(m);

	cudaError_t stat = cudaMalloc((void**)&mData, pitchRows*sizeof(T)*n);

	if(mData == NULL || stat != cudaSuccess)
	{
		assert(0);
		return false;
	}
	mNumRows = m;
	mNumColumns = n;
	mNumPitchedRows = pitchRows;

	stat = cudaMemset((void*)mData, 0, pitchRows*sizeof(T)*n);

	if(stat != cudaSuccess)
	{
		assert(0);
		return false;
	}

	return true;
}


//i - row index, j - column index
//this function is extremely slow and shouldn't be use regularly. It is good for debugging or changing a single element but shouldn't used inside a looptemplate <class T>
template <class T>
inline const T MatrixGPU<T>::get_slow(int i, int j) const
{
	T element;
	checkCudaErrors(cudaMemcpy(&element, mData + j*mNumPitchedRows + i, sizeof(T), cudaMemcpyDeviceToHost));
	return element;
}


//i - row index, j - column index
//this function is extremely slow and shouldn't be use regularly. It is good for debugging or changing a single element but shouldn't used inside a loop
template <class T>
void MatrixGPU<T>::set_slow(int i, int j, const T& value)
{
	assert(i >= 0 && i < mNumRows && j >= 0 && j < mNumColumns);
	checkCudaErrors(cudaMemcpy(mData + j*mNumPitchedRows + i, &value, sizeof(T), cudaMemcpyHostToDevice));
}


template <class T>
inline bool MatrixGPU<T>::isValid() const
{
	return (mData != NULL && mNumRows > 0 && mNumColumns > 0 && mNumPitchedRows > 0);
}

template <class T>
int MatrixGPU<T>::nCols() const
{
	return mNumColumns;
}

template <class T>
int MatrixGPU<T>::nRows() const
{
	return mNumRows;
}

template <class T>
int MatrixGPU<T>::nPitchedRows() const
{
	return mNumPitchedRows;
}

template <class T>
const T* MatrixGPU<T>::getData()
{
	return mData;
}


template <class T>
bool MatrixGPU<T>::mult(const MatrixGPU<T>& A, const MatrixGPU<T>& B, char opA, char opB, double* gpuTime)
{
	int m = 0;
	int n = 0;
	int k = 0;

	if(opA == 'N' && opB == 'N')
	{
		m = A.mNumRows;
		k = A.mNumColumns;
		n = B.mNumColumns;
		if((mNumRows != m) || (mNumColumns != n) || (B.mNumRows != k))
		{
			assert(0);
			return false;
		}
	}
	if(opA != 'N' && opB == 'N')
	{
		m = A.mNumColumns;
		k = A.mNumRows;
		n = B.mNumColumns;
		if((mNumRows != m) || (mNumColumns != n) || (B.mNumRows != k))
		{
			assert(0);
			return false;
		}
	}
	if(opA == 'N' && opB != 'N')
	{
		m = A.mNumRows;
		k = A.mNumColumns;
		n = B.mNumRows;
		if((mNumRows != m) || (mNumColumns != n) || (B.mNumColumns != k))
		{
			assert(0);
			return false;
		}
	}
	if(opA != 'N' && opB != 'N')
	{
		m = A.mNumColumns;
		k = A.mNumRows;
		n = B.mNumRows;
		if((mNumRows != m) || (mNumColumns != n) || (B.mNumColumns != k))
		{
			assert(0);
			return false;
		}
	}
	if(m <= 0 || n <= 0 || k <= 0)
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

	if(opB == 'N' && n == 1) //if B is a column vector
	{
		cublas_gemv(A, B, m, k, opA);
	}
	else
	{
		cublas_gemm(A, B, m, n, k, opA, opB);
	}

	if(gpuTime)
	{
		*gpuTime = timer.stopTimer();
	}

	return true;
}

//result=u dot v, where u is this class
template <class T>
ComplexDouble MatrixGPU<T>::dotColVec(const MatrixGPU<T>& v, double* gpuTime)
{
	ComplexDouble result=Complex(0.0);

	int n = v.mNumRows;
	if ((mNumRows != n) || mNumColumns != 1 ||v.mNumColumns != 1){
		assert(0);
		return result;
	}
	if (n <= 0){
		assert(0);
		return result;
	}

	CUDATimer timer;
	if (gpuTime)
	{
		*gpuTime = 0.0;
		timer.startTimer();
	}

	result=cublas_dot(v, n);

	if (gpuTime)
	{
		*gpuTime = timer.stopTimer();
	}

	return result;
}

template <class T>
double MatrixGPU<T>::norm(double* gpuTime = NULL)
{
	double result = 0;

	int n = mNumRows;
	if (mNumColumns != 1 || n <= 0){
		assert(0);
		return result;
	}

	CUDATimer timer;
	if (gpuTime)
	{
		*gpuTime = 0.0;
		timer.startTimer();
	}

	result = cublas_norm(n);

	if (gpuTime)
	{
		*gpuTime = timer.stopTimer();
	}

	return result;
}

template <class T>
double  MatrixGPU<T>::cublas_norm(int n)
{
	double result = cublasDznrm2(n, (const cuDoubleComplex*)mData, 1);
	cublasStatus err = cublasGetError();

	return result;
}

template <class T>
ComplexDouble MatrixGPU<T>::cublas_dot(const MatrixGPU<ComplexDouble>& v, int n)
{
	//cuDoubleComplex result=cublasZdotu(n, (const cuDoubleComplex*)mData, 1, (const cuDoubleComplex*)v.mData, 1);
	cuDoubleComplex result = cublasZdotc(n, (const cuDoubleComplex*)mData, 1, (const cuDoubleComplex*)v.mData, 1);

	cublasStatus err = cublasGetError();

	ComplexDouble returnResult(result.x, result.y);
	return returnResult;
}


template <class T>
void MatrixGPU<T>::cublas_gemv(const MatrixGPU<float>& A, const MatrixGPU<float>& v, int m, int n, char opA)
{
	if(v.nCols() != 1)
	{
		assert(0);
	}
	float one = 1.0f;
	float zero = 0.0f;

	cublasSgemv(opA, m, n, one, A.mData, A.mNumPitchedRows, v.mData, 1, zero, mData, 1);
	cublasStatus err = cublasGetError();
}


template <class T>
void MatrixGPU<T>::cublas_gemv(const MatrixGPU<double>& A, const MatrixGPU<double>& v, int m, int n, char opA)
{
	if(v.nCols() != 1)
	{
		assert(0);
	}
	double one = 1.0;
	double zero = 0.0;

	cublasDgemv(opA, m, n, one, A.mData, A.mNumPitchedRows, v.mData, 1, zero, mData, 1);
	cublasStatus err = cublasGetError();
}


template <class T>
void MatrixGPU<T>::cublas_gemv(const MatrixGPU<ComplexFloat>& A, const MatrixGPU<ComplexFloat>& v, int m, int n, char opA)
{
	if(v.nCols() != 1)
	{
		assert(0);
	}
	cuFloatComplex one, zero;
	one.x = 1.0f;
	one.y = 0.0f;
	zero.x = 0.0f;
	zero.y = 0.0f;

	cublasCgemv(opA, m, n, one, (const cuFloatComplex*)A.mData, A.mNumPitchedRows, (const cuFloatComplex*)v.mData, 1, zero, (cuFloatComplex*)mData, 1);
	cublasStatus err = cublasGetError();
}


template <class T>
void MatrixGPU<T>::cublas_gemv(const MatrixGPU<ComplexDouble>& A, const MatrixGPU<ComplexDouble>& v, int m, int n, char opA)
{
	if(v.nCols() != 1)
	{
		assert(0);
	}
	cuDoubleComplex one, zero;
	one.x = 1.0;
	one.y = 0.0;
	zero.x = 0.0;
	zero.y = 0.0;

	cublasZgemv(opA, m, n, one, (const cuDoubleComplex*)A.mData, A.mNumPitchedRows, (const cuDoubleComplex*)v.mData, 1, zero, (cuDoubleComplex*)mData, 1);
	cublasStatus err = cublasGetError();
}


template <class T>
void MatrixGPU<T>::cublas_gemm(const MatrixGPU<ComplexFloat>& A, const MatrixGPU<ComplexFloat>& B, int m, int n, int k, char opA, char opB)
{
	cuFloatComplex one, zero;
	one.x = 1.0f;
	one.y = 0.0f;
	zero.x = 0.0f;
	zero.y = 0.0f;

	cublasCgemm(opA, opB, m, n, k, one, (const cuFloatComplex*)A.mData, A.mNumPitchedRows, (const cuFloatComplex*)B.mData, B.mNumPitchedRows, zero, (cuFloatComplex*)mData, mNumPitchedRows);
	cublasStatus err = cublasGetError();
}


template <class T>
void MatrixGPU<T>::cublas_gemm(const MatrixGPU<ComplexDouble>& A, const MatrixGPU<ComplexDouble>& B, int m, int n, int k, char opA, char opB)
{
	cuDoubleComplex one, zero;
	one.x = 1.0;
	one.y = 0.0;
	zero.x = 0.0;
	zero.y = 0.0;

	cublasZgemm(opA, opB, m, n, k, one, (const cuDoubleComplex*)A.mData, A.mNumPitchedRows, (const cuDoubleComplex*)B.mData, B.mNumPitchedRows, zero, (cuDoubleComplex*)mData, mNumPitchedRows);
	cublasStatus err = cublasGetError();
}


template <class T>
void MatrixGPU<T>::cublas_gemm(const MatrixGPU<float>& A, const MatrixGPU<float>& B, int m, int n, int k, char opA, char opB)
{
	float one = 1.0f;
	float zero = 0.0f;

	cublasSgemm(opA, opB, m, n, k, one, A.mData, A.mNumPitchedRows, B.mData, B.mNumPitchedRows, zero, mData, mNumPitchedRows);
	cublasStatus err = cublasGetError();
}


template <class T>
void MatrixGPU<T>::cublas_gemm(const MatrixGPU<double>& A, const MatrixGPU<double>& B, int m, int n, int k, char opA, char opB)
{
	double one = 1.0f;
	double zero = 0.0f;

	cublasDgemm(opA, opB, m, n, k, one, A.mData, A.mNumPitchedRows, B.mData, B.mNumPitchedRows, zero, mData, mNumPitchedRows);
	cublasStatus err = cublasGetError();
}


template <class T>
bool MatrixGPU<T>::add(const MatrixGPU<T>& A, const MatrixGPU<T>& B, double* gpuTime)
{
	if((mNumRows <= 0) || (mNumColumns <= 0) || (mNumRows != A.mNumRows) || (mNumRows != B.mNumRows) || (mNumColumns != A.mNumColumns) || (mNumColumns != B.mNumColumns) || (A.mNumPitchedRows != B.mNumPitchedRows) || (A.mNumPitchedRows != mNumPitchedRows))
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

	bool res = cuda_add(A, B);

	if(gpuTime)
	{
		*gpuTime = timer.stopTimer();
	}

	return true;
}

template <class T>
bool MatrixGPU<T>::cuda_add(const MatrixGPU<float>& A, const MatrixGPU<float>& B)
{
	return cuAddVectorsReal(mNumPitchedRows*mNumColumns, A.mData, B.mData, mData);
}

template <class T>
bool MatrixGPU<T>::cuda_add(const MatrixGPU<ComplexFloat>& A, const MatrixGPU<ComplexFloat>& B)
{
	return cuAddVectorsComplex(mNumPitchedRows*mNumColumns, A.mData, B.mData, mData);
}

template <class T>
bool MatrixGPU<T>::sub(const MatrixGPU<T>& A, const MatrixGPU<T>& B, double* gpuTime)
{
	if((mNumRows <= 0) || (mNumColumns <= 0) || (mNumRows != A.mNumRows) || (mNumRows != B.mNumRows) || (mNumColumns != A.mNumColumns) || (mNumColumns != B.mNumColumns) || (A.mNumPitchedRows != B.mNumPitchedRows) || (A.mNumPitchedRows != mNumPitchedRows))
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

	bool res = cuda_sub(A, B);

	if(gpuTime)
	{
		*gpuTime = timer.stopTimer();
	}

	return true;
}

template <class T>
bool MatrixGPU<T>::cuda_sub(const MatrixGPU<float>& A, const MatrixGPU<float>& B)
{
	return cuSubVectorsReal(mNumPitchedRows*mNumColumns, A.mData, B.mData, mData);
}

template <class T>
bool MatrixGPU<T>::cuda_sub(const MatrixGPU<ComplexFloat>& A, const MatrixGPU<ComplexFloat>& B)
{
	return cuSubVectorsComplex(mNumPitchedRows*mNumColumns, A.mData, B.mData, mData);
}

template <class T>
bool MatrixGPU<T>::cuda_sub(const MatrixGPU<ComplexDouble>& A, const MatrixGPU<ComplexDouble>& B)
{
	return cuSubVectorsComplexDouble(mNumPitchedRows*mNumColumns, A.mData, B.mData, mData);
}
template <class T>
bool MatrixGPU<T>::scale(T alpha, double* gpuTime)
{
	if(mNumColumns <= 0 || mNumRows <= 0)
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

	bool res = cuScaleVectorComplexDouble(mNumPitchedRows*mNumColumns, mData, alpha);

	if(gpuTime)
	{
		*gpuTime = timer.stopTimer();
	}

	return true;
}

template <class T>
bool MatrixGPU<T>::exponent(double* gpuTime)
{
	if(mNumColumns <= 0 || mNumRows <= 0)
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

	bool res = cuExponentVectorComplexDouble(mNumPitchedRows*mNumColumns, mData);

	if(gpuTime)
	{
		*gpuTime = timer.stopTimer();
	}

	return true;
}

template <class T>
int MatrixGPU<T>::roundToCoalesced(int numElements)
{
	if((numElements % 32) == 0)
	{
		return numElements;
	}
	else
	{
		return ((numElements / 32) + 1)*32;
	}
}

template <class T>
bool MatrixGPU<T>::conjugate()
{
	return cuda_conjugate(*this);
}

template <class T>
bool MatrixGPU<T>::cuda_conjugate(MatrixGPU<ComplexFloat>& A) const
{
	return cuConjugate(A.mNumPitchedRows*A.mNumColumns, A.mData);
}

template <class T>
bool MatrixGPU<T>::cuda_conjugate(MatrixGPU<float>& A) const
{
	return true;
}

template <class T>
bool MatrixGPU<T>::add(T alpha)
{
	if(mNumRows <= 0 || mNumColumns <= 0 || mNumPitchedRows <= 0)
	{
		assert(0);
		return false;
	}

	else
	{
		return cuda_add(*this, alpha);
	}
}

template <class T>
bool MatrixGPU<T>::cuda_add(MatrixGPU<ComplexFloat>& A, ComplexFloat alpha) const
{
	return cuAddComplexScalar(A.mNumPitchedRows*A.mNumColumns, A.mData, alpha);
}

template <class T>
bool MatrixGPU<T>::cuda_add(MatrixGPU<float>& A, float alpha) const
{
	return cuAddRealScalar(A.mNumPitchedRows*A.mNumColumns, A.mData, alpha);
}

template <class T>
bool MatrixGPU<T>::sub(T alpha)
{
	return(add(-alpha));
}

template <class T>
bool MatrixGPU<T>::exportToMatlabFile(const std::string& matrixName) const
{
	MatrixCPU<T> CPU_Matrix = *this;
	return CPU_Matrix.exportToMatlabFile(matrixName);
}