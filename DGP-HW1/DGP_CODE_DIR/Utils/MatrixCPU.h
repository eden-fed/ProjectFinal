#pragma once


////////////////////////////////////////////
/////////        MatrixCPU        //////////
///////// Copyright (C) Ofir Weber /////////
////////////////////////////////////////////

#include <windows.h>
#include <Commdlg.h>

#include <Utils\Timer.h>

#include <complex>
#include <vector>
#include <string>
#include <cassert>


#include "Utils\FastMatrix.h"
#include "Utils\Utilities.h"

#include "MatlabInterface.h"

#define MAX_FILENAME_LENGTH 2000

//some non-member global functions

bool convertRealToComplex(MatrixCPU<ComplexFloat>& complexA, const MatrixCPU<float>& realA);
bool convertRealFloatToRealDouble(MatrixCPU<double>& realDoubleA, const MatrixCPU<float>& realFloatA);
bool convertRealDoubleToRealFloat(MatrixCPU<float>& realFloatA, const MatrixCPU<double>& realDoubleA);
bool convertComplexDoubleToComplexFloat(MatrixCPU<ComplexFloat>& complexFloatA, const MatrixCPU<ComplexDouble>& complexDoubleA);
bool convertComplexFloatToComplexDouble(MatrixCPU<ComplexDouble>& complexDoubleA, const MatrixCPU<ComplexFloat>& complexFloatA);

template <class T> void setMatlabScalarVariable(const T& x, const std::string& variableName);
template <class T> void getMatlabScalarVariable(T& x, const std::string& variableName);
template <class T> bool getRealPart(MatrixCPU<T>& realA, const MatrixCPU<std::complex<T>>& complexA);
template <class T> bool getImagPart(MatrixCPU<T>& imagA, const MatrixCPU<std::complex<T>>& complexA);
template <class T> bool setRealImagPart(const MatrixCPU<T>& realA, const MatrixCPU<T>& imagA, MatrixCPU<std::complex<T>>& complexA);
template <class T> bool getRealImagPart(MatrixCPU<T>& realA, MatrixCPU<T>& imagA, const MatrixCPU<std::complex<T>>& complexA);
template <class T> bool convertRealToComplex(MatrixCPU<std::complex<T>>& complexA, const MatrixCPU<T>& realA);


template <class T>
class MatrixCPU
{
	friend class MatrixGPU<T>;
 	friend class MatrixCPU<ComplexFloat>;
	friend class MatrixCPU<ComplexDouble>;

	friend bool convertRealToComplex(MatrixCPU<ComplexFloat>& complexA, const MatrixCPU<float>& realA);
	friend bool convertRealFloatToRealDouble(MatrixCPU<double>& realDoubleA, const MatrixCPU<float>& realFloatA);
	friend bool convertRealDoubleToRealFloat(MatrixCPU<float>& realFloatA, const MatrixCPU<double>& realDoubleA);
	friend bool convertComplexDoubleToComplexFloat(MatrixCPU<ComplexFloat>& complexFloatA, const MatrixCPU<ComplexDouble>& complexDoubleA);
	friend bool convertComplexFloatToComplexDouble(MatrixCPU<ComplexDouble>& complexDoubleA, const MatrixCPU<ComplexFloat>& complexFloatA);
	
 	template <class T> friend bool getRealPart(MatrixCPU<T>& realA, const MatrixCPU<std::complex<T>>& complexA);
 	template <class T> friend bool getImagPart(MatrixCPU<T>& imagA, const MatrixCPU<std::complex<T>>& complexA);
	template <class T> friend bool getRealImagPart(MatrixCPU<T>& realA, MatrixCPU<T>& imagA, const MatrixCPU<std::complex<T>>& complexA);
	template <class T> friend bool convertRealToComplex(MatrixCPU<std::complex<T>>& complexA, const MatrixCPU<T>& realA);
	template <class T> friend bool setRealImagPart(const MatrixCPU<T>& realA, const MatrixCPU<T>& imagA, MatrixCPU<std::complex<T>>& complexA);

public:


	MatrixCPU();
	MatrixCPU(int m, int n, bool clearToZero = true);
	MatrixCPU(const MatrixGPU<T>& mat);
	MatrixCPU(const MatrixCPU<ComplexFloat>& mat);
	MatrixCPU(const MatrixCPU<float>& mat);
	MatrixCPU(const MatrixCPU<ComplexDouble>& mat);
	MatrixCPU(const MatrixCPU<double>& mat);
	void allocateMemory(int numRows, int numColumns);
	void transferData(MatrixCPU<ComplexFloat>& to, const MatrixCPU<ComplexFloat>& from) const;
	void transferData(MatrixCPU<float>& to, const MatrixCPU<float>& from) const;
	void transferData(MatrixCPU<double>& to, const MatrixCPU<double>& from) const;
	void transferData(MatrixCPU<ComplexDouble>& to, const MatrixCPU<ComplexDouble>& from) const;
	void transferData(MatrixCPU<ComplexFloat>& to, const MatrixCPU<ComplexDouble>& from) const;
	void transferData(MatrixCPU<ComplexDouble>& to, const MatrixCPU<double>& from) const;
	void transferData(MatrixCPU<ComplexFloat>& to, const MatrixCPU<double>& from) const;

	const MatrixCPU<T>& operator=(const MatrixCPU<T>& mat);
	const MatrixCPU<T>& operator=(const MatrixGPU<T>& mat);
	~MatrixCPU();
	inline T& operator()(int i, int j);
	inline const T& operator()(int i, int j) const;
	inline bool isValid() const;
	bool resize(int m, int n, bool clearToZero = true);
	bool resizeKeep(int m, int n);
	void freeData();
	int nCols() const;
	int nRows() const;
	bool mult(const MatrixCPU<T>& A, const MatrixCPU<T>& B, char opA = 'N', char opB = 'N', double* cpuTime = NULL);
	bool divide(const MatrixCPU<T>& denominator, double* cpuTime = NULL);
	bool add(const MatrixCPU<T>& A, const MatrixCPU<T>& B, double* cpuTime = NULL);
	bool add(T alpha);
	bool sub(const MatrixCPU<T>& A, const MatrixCPU<T>& B, double* cpuTime = NULL);
	bool sub(T alpha);
	bool scale(T alpha, double* cpuTime = NULL);
	bool exponent();
	T sum() const;
	T maxValue() const;
	double maxAbsValue() const;
	T minValue() const;
	void calcRowSum(std::vector<T>& sum) const;
	double computePseudoInverse(MatrixCPU<T>& Pinv, int* rank = NULL, double tolerance = 0.0) const;
	double computePseudoInverse(MatrixCPU<T>& Pinv, int corank) const;
	bool conjugate();
	bool exportToMatlabFile(const std::string& matrixName = "") const;
	std::string setAsMatlabVariable(const std::string& variableName = "") const;
	void getMatlabVariable(const std::string& variableName);


protected:

	bool conj(MatrixCPU<ComplexFloat>& A);
	bool conj(MatrixCPU<float>& A);
	bool conj(MatrixCPU<ComplexDouble>& A);
	bool conj(MatrixCPU<double>& A);
	float maxValue(const MatrixCPU<float>& A) const;
	ComplexFloat maxValue(const MatrixCPU<ComplexFloat>& A) const;
	float minValue(const MatrixCPU<float>& A) const;
	ComplexFloat minValue(const MatrixCPU<ComplexFloat>& A) const;
	double maxAbsValue(const MatrixCPU<T>& A) const;

	std::string getRandomMatrixName(const std::string& prefix = "") const;

	void setAsMatlabVariable(const MatrixCPU<ComplexDouble>& matrix, const std::string& variableName) const;
	void setAsMatlabVariable(const MatrixCPU<ComplexFloat>& matrix, const std::string& variableName) const;
	void setAsMatlabVariable(const MatrixCPU<double>& matrix, const std::string& variableName) const;
	void setAsMatlabVariable(const MatrixCPU<float>& matrix, const std::string& variableName) const;

	void getMatlabVariable(MatrixCPU<ComplexDouble>& matrix, const std::string& variableName);
	void getMatlabVariable(MatrixCPU<double>& matrix, const std::string& variableName);
	void getMatlabVariable(MatrixCPU<ComplexFloat>& matrix, const std::string& variableName);
	void getMatlabVariable(MatrixCPU<float>& matrix, const std::string& variableName);

	void exportToMatlabFile(const MatrixCPU<float>& A, FILE* fp) const;
	void exportToMatlabFile(const MatrixCPU<ComplexFloat>& A, FILE* fp) const;
	void exportToMatlabFile(const MatrixCPU<double>& A, FILE* fp) const;
	void exportToMatlabFile(const MatrixCPU<ComplexDouble>& A, FILE* fp) const;


protected:

	//data is stored in column major manner (Fortran style)
	T* mData;
	int mNumRows;
	int mNumColumns;
};

template <class T>
MatrixCPU<T>::MatrixCPU()
{
	mData = NULL;
	mNumRows = 0;
	mNumColumns = 0;
}

template <class T>
MatrixCPU<T>::MatrixCPU(int m, int n, bool clearToZero = true)
{
	mData = NULL;
	mNumRows = 0;
	mNumColumns = 0;

	// zl: I'd like to allow zero size matrices, they fit well with the logic, and fit with the empty constructor.
	//if(m <= 0 || n <= 0)
	if (m < 0 || n < 0)
	{
		assert(0);
		return;
	}

	mData = (T*)malloc(sizeof(T)*m*n); //debug - mallochost? paged locked memory
	if(mData == NULL)
	{
		assert(0);
		return;
	}
	else
	{
		mNumRows = m;
		mNumColumns = n;
		if(clearToZero)
		{
			memset(mData, 0, sizeof(T)*m*n);
		}
	}
}

//copy c'tor CPU to CPU
template <class T>
MatrixCPU<T>::MatrixCPU(const MatrixCPU<ComplexFloat>& mat)
{
	allocateMemory(mat.mNumRows, mat.mNumColumns);
	transferData(*this, mat);
}

//copy c'tor CPU to CPU
template <class T>
MatrixCPU<T>::MatrixCPU(const MatrixCPU<float>& mat)
{
	allocateMemory(mat.mNumRows, mat.mNumColumns);
	transferData(*this, mat);
}

//copy c'tor CPU to CPU
template <class T>
MatrixCPU<T>::MatrixCPU(const MatrixCPU<ComplexDouble>& mat)
{
	allocateMemory(mat.nRows(), mat.nCols());
	transferData(*this, mat);
}

//copy c'tor CPU to CPU
template <class T>
MatrixCPU<T>::MatrixCPU(const MatrixCPU<double>& mat)
{
	allocateMemory(mat.mNumRows, mat.mNumColumns);
	transferData(*this, mat);
}

template <class T>
void MatrixCPU<T>::allocateMemory(int numRows, int numColumns)
{
	mData = NULL;
	mNumRows = 0;
	mNumColumns = 0;

	if(numRows == 0 && numColumns == 0) //an empty matrix
	{
		return;
	}

	if(numRows <= 0 || numColumns <= 0)
	{
		assert(0);
		return;
	}

	mData = (T*)malloc(sizeof(T)*numRows*numColumns);
	if(mData == NULL)
	{
		assert(0);
		return;
	}
	mNumRows = numRows;
	mNumColumns = numColumns;
}

template <class T>
void MatrixCPU<T>::transferData(MatrixCPU<ComplexDouble>& to, const MatrixCPU<ComplexDouble>& from) const
{
	memcpy(to.mData, from.mData, from.mNumRows*from.mNumColumns*sizeof(ComplexDouble));
}

template <class T>
void MatrixCPU<T>::transferData(MatrixCPU<ComplexFloat>& to, const MatrixCPU<ComplexFloat>& from) const
{
	memcpy(to.mData, from.mData, from.mNumRows*from.mNumColumns*sizeof(ComplexFloat));
}

template <class T>
void MatrixCPU<T>::transferData(MatrixCPU<float>& to, const MatrixCPU<float>& from) const
{
	memcpy(to.mData, from.mData, from.mNumRows*from.mNumColumns*sizeof(float));
}

template <class T>
void MatrixCPU<T>::transferData(MatrixCPU<double>& to, const MatrixCPU<double>& from) const
{
	memcpy(to.mData, from.mData, from.mNumRows*from.mNumColumns*sizeof(double));
}

template <class T>
void MatrixCPU<T>::transferData(MatrixCPU<ComplexFloat>& to, const MatrixCPU<ComplexDouble>& from) const
{
	int numElements = from.mNumRows*from.mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		to.mData[i]._Val[0] = (float)from.mData[i]._Val[0];
		to.mData[i]._Val[1] = (float)from.mData[i]._Val[1];
	}
}

template <class T>
void MatrixCPU<T>::transferData(MatrixCPU<ComplexDouble>& to, const MatrixCPU<double>& from) const
{
	int numElements = from.mNumRows*from.mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		to.mData[i]._Val[0] = from.mData[i];
		to.mData[i]._Val[1] = 0.0;
	}
}

template <class T>
void MatrixCPU<T>::transferData(MatrixCPU<ComplexFloat>& to, const MatrixCPU<double>& from) const
{
	int numElements = from.mNumRows*from.mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		to.mData[i]._Val[0] = (float)from.mData[i];
		to.mData[i]._Val[1] = 0.0f;
	}
}

//copy c'tor GPU to CPU
template <class T>
MatrixCPU<T>::MatrixCPU(const MatrixGPU<T>& mat)
{
	mData = NULL;
	mNumRows = 0;
	mNumColumns = 0;

	if(mat.mNumRows == 0 && mat.mNumColumns == 0 && mat.mNumPitchedRows == 0) //an empty matrix
	{
		return;
	}

	if(mat.mNumRows <= 0 || mat.mNumColumns <= 0 || mat.mNumPitchedRows <= 0)
	{
		assert(0);
		return;
	}

	mData = (T*)malloc(sizeof(T)*mat.mNumRows*mat.mNumColumns);
	if(mData == NULL)
	{
		assert(0);
		return;
	}
	mNumRows = mat.mNumRows;
	mNumColumns = mat.mNumColumns;

	checkCudaErrors(safeCudaMemcpy2D(mData, mNumRows*sizeof(T), mat.mData, mat.mNumPitchedRows*sizeof(T), mat.mNumRows*sizeof(T), mat.mNumColumns, cudaMemcpyDeviceToHost));
}

//assignment operator CPU to CPU
template <class T>
const MatrixCPU<T>& MatrixCPU<T>::operator=(const MatrixCPU<T>& mat)
{
	if(mat.mNumRows <= 0 || mat.mNumColumns <= 0)
	{
		assert(0);
		return *this;
	}

	if(mat.mNumRows*mat.mNumColumns != mNumRows*mNumColumns)
	{
		if(mData)
		{
			free(mData);
			mData = NULL;
		}
		mNumRows = 0;
		mNumColumns = 0;


		mData = (T*)malloc(sizeof(T)*mat.mNumRows*mat.mNumColumns);
		if(mData == NULL)
		{
			assert(0);
			return *this;
		}
	}

	mNumRows = mat.mNumRows;
	mNumColumns = mat.mNumColumns;

	memcpy(mData, mat.mData, mat.mNumRows*mat.mNumColumns*sizeof(T));
	return *this;
}

//assignment operator GPU to CPU
template <class T>
const MatrixCPU<T>& MatrixCPU<T>::operator=(const MatrixGPU<T>& mat)
{
	if(mat.mNumRows <= 0 || mat.mNumColumns <= 0 || mat.mNumPitchedRows <= 0)
	{
		assert(0);
		return *this;
	}

	if(mat.mNumRows*mat.mNumColumns != mNumRows*mNumColumns)
	{
		if(mData)
		{
			free(mData);
			mData = NULL;
		}
		mNumRows = 0;
		mNumColumns = 0;

		mData = (T*)malloc(sizeof(T)*mat.mNumRows*mat.mNumColumns);
		if(mData == NULL)
		{
			assert(0);
			return *this;
		}
	}

	mNumRows = mat.mNumRows;
	mNumColumns = mat.mNumColumns;

	checkCudaErrors(safeCudaMemcpy2D(mData, mNumRows*sizeof(T), mat.mData, mat.mNumPitchedRows*sizeof(T), mat.mNumRows*sizeof(T), mat.mNumColumns, cudaMemcpyDeviceToHost));
	return *this;
}

template <class T>
MatrixCPU<T>::~MatrixCPU()
{
	if(mData)
	{
		free(mData);
		mData = NULL;
	}
}

template <class T>
void MatrixCPU<T>::freeData()
{
	if(mData)
	{
		free(mData);
		mData = NULL;
	}
	mNumRows = 0;
	mNumColumns = 0;
}

template <class T>
bool MatrixCPU<T>::resize(int m, int n, bool clearToZero = true)
{
	// zl: I'd like to allow zero size matrices, they fit well with the logic.
	//if(m <= 0 || n <= 0)
	if (m < 0 || n < 0)
	{
		assert(0);
		return false;
	}

	if(mData)
	{
		free(mData);
		mData = NULL;
		mNumRows = 0;
		mNumColumns = 0;
	}

	mData = (T*)malloc(sizeof(T)*m*n);
	if(mData != NULL)
	{
		mNumRows = m;
		mNumColumns = n;
		if(clearToZero)
		{
			memset(mData, 0, sizeof(T)*m*n);
		}
		return true;
	}
	else
	{
		assert(0);
		return false;
	}
}

template <class T>
bool MatrixCPU<T>::resizeKeep(int m, int n)
{
	assert(mData != NULL);

	if(m <= 0 || n <= 0 || mNumRows <= 0 || mNumColumns <= 0)
	{
		assert(0);
		return false;
	}
	if(m == mNumRows && n == mNumColumns)
	{
		return false;
	}

	T* tempDataPointer = NULL;
	tempDataPointer = (T*)malloc(sizeof(T)*m*n);
	if(tempDataPointer == NULL)
	{
		assert(0);
		return false;
	}
	
	for(int j = 0; j < n; j++)
	{
		for(int i = 0; i < m; i++)
		{
			if(i < mNumRows && j < mNumColumns)
			{
				tempDataPointer[j*m + i] = (*this)(i, j);
			}
			else
			{
				tempDataPointer[j*m + i] = 0;
			}
		}
	}
	mNumRows = m;
	mNumColumns = n;
	free(mData);
	mData = tempDataPointer;

	return true;
}

//i - row index, j - column index
template <class T>
inline T& MatrixCPU<T>::operator()(int i, int j)
{
	assert(i >= 0 && j >= 0 && i < nRows() && j < nCols());
	return mData[j*mNumRows + i];
}

//i - row index, j - column index
template <class T>
inline const T& MatrixCPU<T>::operator()(int i, int j) const
{
	assert(i >= 0 && j >= 0 && i < nRows() && j < nCols());
	return mData[j*mNumRows + i];
}

template <class T>
inline bool MatrixCPU<T>::isValid() const
{
	return (mData != NULL && mNumRows > 0 && mNumColumns > 0);
}

template <class T>
int MatrixCPU<T>::nCols() const
{
	return mNumColumns;
}

template <class T>
int MatrixCPU<T>::nRows() const
{
	return mNumRows;
}

//element-wise division of the current matrix by the denominator matrix.
//if denominator is a vector rather than a matrix, each row of current matrix is divided by the same corresponding element in denominator
template <class T>
bool MatrixCPU<T>::divide(const MatrixCPU<T>& denominator, double* cpuTime)
{
	CPUTimer timer;
	if(cpuTime)
	{
		*cpuTime = 0.0;
		timer.tic();
	}

	int m = denominator.mNumRows;
	int n = denominator.mNumColumns;
	if((mNumRows == m) && (mNumColumns == n))
	{
		for(int j = 0; j < mNumColumns; j++)
		{
			for(int i = 0; i < mNumRows; i++)
			{
				mData[j*mNumRows + i] /= denominator(i, j);
			}
		}
	}
	else if((mNumRows == m) && (n == 1))
	{
		for(int j = 0; j < mNumColumns; j++)
		{
			for(int i = 0; i < mNumRows; i++)
			{
				mData[j*mNumRows + i] /= denominator(i, 0);
			}
		}
	}
	else
	{
		assert(0);
		return false;
	}

	if(cpuTime)
	{
		*cpuTime = timer.toc();
	}

	return true;
}

//C = A*B, where C is this class
template <class T>
bool MatrixCPU<T>::mult(const MatrixCPU<T>& A, const MatrixCPU<T>& B, char opA, char opB, double* cpuTime)
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
	
	CPUTimer timer;
	if(cpuTime)
	{
		*cpuTime = 0.0;
		timer.tic();
	}
	
	//data is being transferred to the GPU
	MatrixGPU<T> A_gpu = A;
	MatrixGPU<T> B_gpu = B;
	MatrixGPU<T> C_gpu(mNumRows, mNumColumns);

	//multiplication is performed on the GPU
	C_gpu.mult(A_gpu, B_gpu, opA, opB);

	//result is being copied back from GPU to CPU
	*this = C_gpu;


	if(cpuTime)
	{
		*cpuTime = timer.toc();
	}
	
	return true;
}


//C = A + B, where C is this class
template <class T>
bool MatrixCPU<T>::add(const MatrixCPU<T>& A, const MatrixCPU<T>& B, double* cpuTime)
{
	if((mNumRows <= 0) || (mNumColumns <= 0) || (mNumRows != A.mNumRows) || (mNumRows != B.mNumRows) || (mNumColumns != A.mNumColumns) || (mNumColumns != B.mNumColumns))
	{
		assert(0);
		return false;
	}
	CPUTimer timer;
	if(cpuTime)
	{
		*cpuTime = 0.0;
		timer.tic();
	}

	int numElements = mNumRows*mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		mData[i] = A.mData[i] + B.mData[i];
	}

	if(cpuTime)
	{
		*cpuTime = timer.toc();
	}

	return true;
}


//C = A - B, where C is this class
template <class T>
bool MatrixCPU<T>::sub(const MatrixCPU<T>& A, const MatrixCPU<T>& B, double* cpuTime)
{
	if((mNumRows <= 0) || (mNumColumns <= 0) || (mNumRows != A.mNumRows) || (mNumRows != B.mNumRows) || (mNumColumns != A.mNumColumns) || (mNumColumns != B.mNumColumns))
	{
		assert(0);
		return false;
	}

	CPUTimer timer;
	if(cpuTime)
	{
		*cpuTime = 0.0;
		timer.tic();
	}

	int numElements = mNumRows*mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		mData[i] = A.mData[i] - B.mData[i];
	}

	if(cpuTime)
	{
		*cpuTime = timer.toc();
	}

	return true;
}

template <class T>
bool MatrixCPU<T>::scale(T alpha, double* cpuTime)
{
	if((mNumRows <= 0) || (mNumColumns <= 0))
	{
		assert(0);
		return false;
	}

	CPUTimer timer;
	if(cpuTime)
	{
		*cpuTime = 0.0;
		timer.tic();
	}

	int numElements = mNumRows*mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		mData[i] = alpha*mData[i];
	}

	if(cpuTime)
	{
		*cpuTime = timer.toc();
	}
	return true;
}

template <class T>
bool MatrixCPU<T>::exponent()
{
	if(mNumColumns <= 0 || mNumRows <= 0)
	{
		assert(0);
		return false;
	}


	int numElements = mNumRows*mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		mData[i] = exp(mData[i]);
	}

	return true;
}

template <class T>
T MatrixCPU<T>::sum() const
{
	int numElements = mNumRows*mNumColumns;

	T sum = 0;
	for(int i = 0; i < numElements; i++)
	{
		sum += mData[i];
	}
	return sum;
}


template <class T>
T MatrixCPU<T>::maxValue() const
{
	if(mNumRows <= 0 || mNumColumns <= 0)
	{
		assert(0);
	}

	return maxValue(*this);
}

template <class T>
T MatrixCPU<T>::minValue() const
{
	if(mNumRows <= 0 || mNumColumns <= 0)
	{
		assert(0);
	}

	return minValue(*this);
}

template <class T>
float MatrixCPU<T>::maxValue(const MatrixCPU<float>& A) const
{
	float maxValue = A.mData[0];

	int numElements = A.mNumRows*A.mNumColumns;

	for(int i = 1; i < numElements; i++)
	{
		maxValue = max(maxValue, mData[i]);
	}
	return maxValue;
}

//for complex data types, we treat the real and imag parts as separate numbers and find both the maximum of the real numbers and the maximum of the imaginary numbers
template <class T>
ComplexFloat MatrixCPU<T>::maxValue(const MatrixCPU<ComplexFloat>& A) const
{

	float maxReal = A.mData[0]._Val[0];
	float maxImag = A.mData[0]._Val[1];

	int numElements = A.mNumRows*A.mNumColumns;

	for(int i = 1; i < numElements; i++)
	{
		ComplexFloat value = mData[i];
		float real = value._Val[0];
		float imag = value._Val[1];
		maxReal	= max(maxReal, real);
		maxImag	= max(maxImag, imag);
	}
	return ComplexFloat(maxReal, maxImag);
}

template <class T>
float MatrixCPU<T>::minValue(const MatrixCPU<float>& A) const
{
	float minValue = A.mData[0];

	int numElements = A.mNumRows*A.mNumColumns;

	for(int i = 1; i < numElements; i++)
	{
		minValue = min(minValue, mData[i]);
	}
	return minValue;
}

//for complex data types, we treat the real and imag parts as separate numbers and find both the minimum of the real numbers and the minimum of the imaginary numbers
template <class T>
ComplexFloat MatrixCPU<T>::minValue(const MatrixCPU<ComplexFloat>& A) const
{

	float minReal = A.mData[0]._Val[0];
	float minImag = A.mData[0]._Val[1];

	int numElements = A.mNumRows*A.mNumColumns;

	for(int i = 1; i < numElements; i++)
	{
		ComplexFloat value = mData[i];
		float real = value._Val[0];
		float imag = value._Val[1];
		minReal	= min(minReal, real);
		minImag	= min(minImag, imag);
	}
	return ComplexFloat(minReal, minImag);
}


template <class T>
double MatrixCPU<T>::maxAbsValue() const
{
	if(mNumRows <= 0 || mNumColumns <= 0)
	{
		assert(0);
	}

	return maxAbsValue(*this);
}


template <class T>
double MatrixCPU<T>::maxAbsValue(const MatrixCPU<T>& A) const
{
	double maxAbsValue = abs(mData[0]);

	int numElements = A.nRows()*A.nCols();

	for(int i = 1; i < numElements; i++)
	{
		maxAbsValue = max(maxAbsValue, abs(mData[i]));
	}
	return maxAbsValue;
}


template <class T>
double MatrixCPU<T>::computePseudoInverse(MatrixCPU<T>& Pinv, int* rank = NULL, double tolerance = 0.0) const
{
	//bool doTimer = true;
	bool doTimer = false;

	CPUTimer timer;

	if(doTimer)
	{
		timer.tic();
	}

	int m = mNumRows;
	int n = mNumColumns;

	if(m <= 0 || n <= 0)
	{
		assert(0);
		return 0.0;
	}

	std::string matrixName = setAsMatlabVariable();
	std::string pinvMatrixName = getRandomMatrixName("Pinv_");
	std::string rankPinvName = getRandomMatrixName("rankPinv_");
	std::string condPinvName = getRandomMatrixName("condPinv_");

	//[X, rankPinv, condPinv] = pinv_extra(A, tol)
	std::string command = "[" + pinvMatrixName + ", " + rankPinvName + ", " + condPinvName + "] = pinv_extra(" + matrixName;

	if(tolerance != 0.0)
	{
		std::ostringstream oss;
		oss << tolerance;
		command += ", " + oss.str();
	}

	command += ");";

	MatlabInterface::GetEngine().EvalToCout(command.c_str());

	Pinv.getMatlabVariable(pinvMatrixName);

	MatrixCPU<double> cond; //this is a double scalar (not a matrix)
	cond.getMatlabVariable(condPinvName);
	assert(cond.nCols() == 1 && cond.nRows() == 1);

	if(rank != NULL)
	{
		MatrixCPU<double> rank_; //this is an integer scalar (not a matrix). the reason we use double rather than int is because Matlab stores everything in double (by default)
		rank_.getMatlabVariable(rankPinvName);
		assert(rank_.nCols() == 1 && rank_.nRows() == 1);

		*rank = (int)rank_(0, 0);

		assert(*rank == rank_(0, 0));
		assert(*rank > 0);
	}


	//bool doClear = false;
	bool doClear = true;

	if(doClear)
	{
		//clear the temporary matrices from Matlab
		command = "clear " + matrixName + " " + pinvMatrixName + " " + rankPinvName + " " + condPinvName + ";";
		MatlabInterface::GetEngine().EvalToCout(command.c_str());
	}

	if(doTimer)
	{
		double cpuTime = timer.toc();
		cout << "computePseudoInverse time: " << cpuTime << " ms" << endl;
	}

	return cond(0, 0);
}


template <class T>
double MatrixCPU<T>::computePseudoInverse(MatrixCPU<T>& Pinv, int corank) const
{
	//bool doTimer = true;
	bool doTimer = false;

	CPUTimer timer;

	if(doTimer)
	{
		timer.tic();
	}

	int m = mNumRows;
	int n = mNumColumns;

	if(m <= 0 || n <= 0)
	{
		assert(0);
		return 0.0;
	}

	std::string matrixName = setAsMatlabVariable();
	std::string pinvMatrixName = getRandomMatrixName("Pinv_");
	std::string condPinvName = getRandomMatrixName("condPinv_");

	//[X, condPinv] = pinv_corank(A, corank)
	std::string command = "[" + pinvMatrixName + ", " + condPinvName + "] = pinv_corank(" + matrixName;

	std::ostringstream oss;
	oss << corank;
	command += ", " + oss.str() + ");";

	MatlabInterface::GetEngine().EvalToCout(command.c_str());

	Pinv.getMatlabVariable(pinvMatrixName);

	MatrixCPU<double> cond; //this is a double scalar (not a matrix)
	cond.getMatlabVariable(condPinvName);
	assert(cond.nCols() == 1 && cond.nRows() == 1);

	//bool doClear = false;
	bool doClear = true;

	if(doClear)
	{
		//clear the temporary matrices from Matlab
		command = "clear " + matrixName + " " + pinvMatrixName + " " + condPinvName + ";";
		MatlabInterface::GetEngine().EvalToCout(command.c_str());
	}


	if(doTimer)
	{
		double cpuTime = timer.toc();
		cout << "computePseudoInverse (with specified corank) time: " << cpuTime << " ms" << endl;
	}

	return cond(0, 0);
}


template<class T>
void MatrixCPU<T>::calcRowSum(std::vector<T>& sum) const
{
	sum.clear();
	sum.resize(mNumRows);
	for(int j = 0; j < mNumColumns; j++)
	{
		for(int i = 0; i < mNumRows; i++)
		{
			sum[i] += mData[j*mNumRows + i];
		}
	}
}

template <class T>
bool MatrixCPU<T>::conjugate()
{
	if(mNumRows <= 0 || mNumColumns <= 0)
	{
		assert(0);
		return false;
	}

	else
	{
		return conj(*this);
	}
}

template <class T>
bool MatrixCPU<T>::conj(MatrixCPU<ComplexFloat>& A)
{
	int numElements = 2*A.mNumRows*A.mNumColumns;

	float* floatData = (float*)A.mData;

	for(int index = 1; index < numElements; index+=2)
	{
		floatData[index] = -floatData[index];
	}

	return true;
}

template <class T>
bool MatrixCPU<T>::conj(MatrixCPU<float>& A)
{
	return true;
}

template <class T>
bool MatrixCPU<T>::conj(MatrixCPU<ComplexDouble>& A)
{
	int numElements = 2*A.mNumRows*A.mNumColumns;

	double* doubleData = (double*)A.mData;

	for(int index = 1; index < numElements; index+=2)
	{
		doubleData[index] = -doubleData[index];
	}

	return true;
}

template <class T>
bool MatrixCPU<T>::conj(MatrixCPU<double>& A)
{
	return true;
}

template <class T>
bool MatrixCPU<T>::add(T alpha)
{
	if(mNumRows <= 0 || mNumColumns <= 0)
	{
		assert(0);
		return false;
	}

	int numElements = mNumRows*mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		mData[i] = mData[i] + alpha;
	}
	return true;
}

template <class T>
bool MatrixCPU<T>::sub(T alpha)
{
	return add(-alpha);
}


template <class T>
std::string MatrixCPU<T>::getRandomMatrixName(const std::string& prefix = "") const
{
	//uintptr_t id = (uintptr_t)this;
	std::ostringstream oss;
	int id1 = getRandomNumber();
	int id2 = getRandomNumber();
	oss << id1 << id2;
	std::string name("M_");
	name += prefix;
	name += oss.str();

	return name;
}


template <class T>
std::string MatrixCPU<T>::setAsMatlabVariable(const std::string& variableName = "") const
{
	std::string matlabVariableName;

	if(variableName == "")
	{
		matlabVariableName = getRandomMatrixName();
	}
	else
	{
		matlabVariableName = variableName;
	}
	setAsMatlabVariable(*this, matlabVariableName);
	
	return matlabVariableName;
}


template <class T>
void MatrixCPU<T>::setAsMatlabVariable(const MatrixCPU<ComplexDouble>& matrix, const std::string& variableName) const
{
	MatlabInterface::GetEngine(false).SetEngineComplexMatrix<double>(variableName.data(), matrix.mNumRows, matrix.mNumColumns, matrix.mData, true);
}

template <class T>
void MatrixCPU<T>::setAsMatlabVariable(const MatrixCPU<ComplexFloat>& matrix, const std::string& variableName) const
{
	MatlabInterface::GetEngine(false).SetEngineComplexMatrix<float>(variableName.data(), matrix.mNumRows, matrix.mNumColumns, matrix.mData, true);
}

template <class T>
void MatrixCPU<T>::setAsMatlabVariable(const MatrixCPU<double>& matrix, const std::string& variableName) const
{
	MatlabInterface::GetEngine(false).SetEngineRealMatrix<double>(variableName.data(), matrix.mNumRows, matrix.mNumColumns, matrix.mData, true);
}

template <class T>
void MatrixCPU<T>::setAsMatlabVariable(const MatrixCPU<float>& matrix, const std::string& variableName) const
{
	MatlabInterface::GetEngine(false).SetEngineRealMatrix<float>(variableName.data(), matrix.mNumRows, matrix.mNumColumns, matrix.mData, true);
}

template <class T>
void MatrixCPU<T>::getMatlabVariable(const std::string& variableName)
{
	unsigned int m = 0;
	unsigned int n = 0;
	MatlabInterface::GetEngine(false).GetMatrixDimensions(variableName.data(), m, n);
	if(m != 0 && n != 0)
	{
		resize(m, n, false);
		getMatlabVariable(*this, variableName);
	}
}

template <class T>
void MatrixCPU<T>::getMatlabVariable(MatrixCPU<ComplexDouble>& matrix, const std::string& variableName)
{
	MatlabInterface::GetEngine(false).GetEngineComplexMatrix<double>(variableName.data(), matrix.mNumRows, matrix.mNumColumns, matrix.mData, true);
}

template <class T>
void MatrixCPU<T>::getMatlabVariable(MatrixCPU<double>& matrix, const std::string& variableName)
{
	MatlabInterface::GetEngine(false).GetEngineRealMatrix<double>(variableName.data(), matrix.mNumRows, matrix.mNumColumns, matrix.mData, true);
}

template <class T>
void MatrixCPU<T>::getMatlabVariable(MatrixCPU<ComplexFloat>& matrix, const std::string& variableName)
{
	MatlabInterface::GetEngine(false).GetEngineComplexMatrix<float>(variableName.data(), matrix.mNumRows, matrix.mNumColumns, matrix.mData, true);
}

template <class T>
void MatrixCPU<T>::getMatlabVariable(MatrixCPU<float>& matrix, const std::string& variableName)
{
	MatlabInterface::GetEngine(false).GetEngineRealMatrix<float>(variableName.data(), matrix.mNumRows, matrix.mNumColumns, matrix.mData, true);
}

template <class T>
void MatrixCPU<T>::exportToMatlabFile(const MatrixCPU<float>& A, FILE* fp) const
{
	for(int i = 0; i < A.mNumRows; i++)
	{
		for(int j = 0; j < A.mNumColumns; j++)
		{
			float f = A.mData[i + j*A.mNumRows];
			fprintf(fp, "%.30f", f);
			if(j < A.mNumColumns - 1)
			{
				fprintf(fp, " ");
			}
			else
			{
				fprintf(fp, "\n");
			}
		}
	}
}

template <class T>
void MatrixCPU<T>::exportToMatlabFile(const MatrixCPU<ComplexFloat>& A, FILE* fp) const
{
	for(int i = 0; i < A.mNumRows; i++)
	{
		for(int j = 0; j < A.mNumColumns; j++)
		{
			ComplexFloat C = A.mData[i + j*A.mNumRows];
			fprintf(fp, "%.30f%+.30fi", C._Val[0], C._Val[1]);
			if(j < A.mNumColumns - 1)
			{
				fprintf(fp, " ");
			}
			else
			{
				fprintf(fp, "\n");
			}
		}
	}
}

template <class T>
void MatrixCPU<T>::exportToMatlabFile(const MatrixCPU<double>& A, FILE* fp) const
{
	for(int i = 0; i < A.mNumRows; i++)
	{
		for(int j = 0; j < A.mNumColumns; j++)
		{
			double d = A.mData[i + j*A.mNumRows];
			fprintf(fp, "%.30lf", d);
			if(j < A.mNumColumns - 1)
			{
				fprintf(fp, " ");
			}
			else
			{
				fprintf(fp, "\n");
			}
		}
	}
}

template <class T>
void MatrixCPU<T>::exportToMatlabFile(const MatrixCPU<ComplexDouble>& A, FILE* fp) const
{
	for(int i = 0; i < A.mNumRows; i++)
	{
		for(int j = 0; j < A.mNumColumns; j++)
		{
			ComplexDouble C = A.mData[i + j*A.mNumRows];
			fprintf(fp, "%.30lf%+.30lfi", C._Val[0], C._Val[1]);
			if(j < A.mNumColumns - 1)
			{
				fprintf(fp, " ");
			}
			else
			{
				fprintf(fp, "\n");
			}
		}
	}
}

//in order to load the file into Matlab use:  A = dlmread('C:\Users\Ofir\Desktop\FileName.bat');
template <class T>
bool MatrixCPU<T>::exportToMatlabFile(const std::string& matrixName) const
{
	char staticFileName[MAX_FILENAME_LENGTH] = {0};
	matrixName.copy(staticFileName, matrixName.size());

	OPENFILENAMEA ofn =
	{
		sizeof(OPENFILENAMEA),
		NULL,
		NULL, 
		"Matlab Matrix in Ascii format (*.bat)\0*.bat\0All Files (*.*)\0*.*\0\0",
		NULL, 0,
		1,
		staticFileName, sizeof(staticFileName),
		NULL, 0,
		NULL,
		"Export matrix to Matlab as Ascii file ...",
		0,
		0,
		0,
		"bat",
		0, 
		NULL, 
		NULL
	};

	if(!GetSaveFileNameA(&ofn))
	{
		return false;
	}

	FILE* fp = NULL;

	if((fp = fopen(staticFileName, "w")) == NULL)
	{
		return false;
	}
	
	exportToMatlabFile(*this, fp);
	
	fclose(fp);
	
	return true;
}


/////////////////  Global template functions  ///////////////////////////

template <class T>
bool setRealImagPart(const MatrixCPU<T>& realA, const MatrixCPU<T>& imagA, MatrixCPU<std::complex<T>>& complexA)
{
	if(complexA.mNumRows != realA.mNumRows || complexA.mNumColumns != realA.mNumColumns || complexA.mNumRows != imagA.mNumRows || complexA.mNumColumns != imagA.mNumColumns)
	{
		assert(0);
		return false;
	}

	int numElements = complexA.mNumRows*complexA.mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		std::complex<T> C(realA.mData[i], imagA.mData[i]);
		complexA.mData[i] = C;
	}

	return true;
}

template <class T>
bool getRealImagPart(MatrixCPU<T>& realA, MatrixCPU<T>& imagA, const MatrixCPU<std::complex<T>>& complexA)
{
	if(complexA.mNumRows != realA.mNumRows || complexA.mNumColumns != realA.mNumColumns || complexA.mNumRows != imagA.mNumRows || complexA.mNumColumns != imagA.mNumColumns)
	{
		assert(0);
		return false;
	}

	int numElements = complexA.mNumRows*complexA.mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		std::complex<T> C = complexA.mData[i];
		realA.mData[i] = C._Val[0];
		imagA.mData[i] = C._Val[1];
	}

	return true;
}

template <class T>
bool getRealPart(MatrixCPU<T>& realA, const MatrixCPU<std::complex<T>>& complexA)
{
	if(complexA.mNumRows != realA.mNumRows || complexA.mNumColumns != realA.mNumColumns)
	{
		assert(0);
		return false;
	}

	int numElements = complexA.mNumRows*complexA.mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		realA.mData[i] = ((T*)(complexA.mData))[2*i];
	}

	return true;
}

template <class T>
bool getImagPart(MatrixCPU<T>& imagA, const MatrixCPU<std::complex<T>>& complexA)
{
	if(complexA.mNumRows != imagA.mNumRows || complexA.mNumColumns != imagA.mNumColumns)
	{
		assert(0);
		return false;
	}

	int numElements = complexA.mNumRows*complexA.mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		imagA.mData[i] = ((T*)(complexA.mData))[2*i + 1];
	}

	return true;
}

template <class T>
bool convertRealToComplex(MatrixCPU<std::complex<T>>& complexA, const MatrixCPU<T>& realA)
{
	if(complexA.mNumRows != realA.mNumRows || complexA.mNumColumns != realA.mNumColumns)
	{
		assert(0);
		return false;
	}

	int numElements = realA.mNumRows*realA.mNumColumns;

	for(int i = 0; i < numElements; i++)
	{
		std::complex<T> C;
		C._Val[0] = realA.mData[i];
		C._Val[1] = 0.0f;

		complexA.mData[i] = C;
	}
	return true;
}


template <class T>
void setMatlabScalarVariable(const T& x, const std::string& variableName)
{
	MatrixCPU<T> matrix(1, 1);
	matrix(0, 0) = x;
	matrix.setAsMatlabVariable(variableName);
}


template <class T>
void getMatlabScalarVariable(T& x, const std::string& variableName)
{
	MatrixCPU<T> matrix(1, 1);
	matrix.getMatlabVariable(variableName);
	x = matrix(0, 0);
}

