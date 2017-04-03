#pragma once

#include <complex>

#define I_DIV_UP(a, b) ((((a) % (b)) != 0) ? ((a) / (b) + 1) : ((a) / (b)))

extern "C"
{
	//this is obsolete. it used to be faster than cublas due to inefficient implementation of cublas but now cublas is optimal and there is no reason to maintain this function
	bool cuMultMatrixByVectorComplex(int m, int n, int numPitchedRows, const std::complex<float>* d_A, const std::complex<float>* d_X, std::complex<float>* d_b);
	//this is obsolete. it used to be faster than cublas due to inefficient implementation of cublas but now cublas is optimal and there is no reason to maintain this function
	bool cuMultMatrixByVectorReal(int m, int n, int numPitchedRows, const float* d_A, const float* d_X, float* d_b);
	bool cuAddVectorsReal(int numElements, const float* d_A, const float* d_B, float* d_C);
	bool cuAddVectorsComplex(int numElements, const std::complex<float>* d_A, const std::complex<float>* d_B, std::complex<float>* d_C);
	bool cuSubVectorsComplex(int numElements, const std::complex<float>* d_A, const std::complex<float>* d_B, std::complex<float>* d_C);
	bool cuSubVectorsComplexDouble(int numElements, const std::complex<double>* d_A, const std::complex<double>* d_B, std::complex<double>* d_C);//***
	bool cuSubVectorsReal(int numElements, const float* d_A, const float* d_B, float* d_C);
	bool cuScaleVectorComplex(int numElements, std::complex<float>* d_A, std::complex<float> alpha);
	bool cuScaleVectorComplexDouble(int numElements, std::complex<double>* d_A, std::complex<double> alpha);//***
	bool cuExponentVectorComplex(int numElements, std::complex<float>* d_A);
	bool cuExponentVectorComplexDouble(int numElements, std::complex<double>* d_A);//***
	bool cuConvertRealToComplex(int numElements, std::complex<float>* d_A, const float* d_B);
	bool cuConjugate(int numElements, std::complex<float>* d_A);
	bool cuAddRealScalar(int numElements, float* d_A, float alpha);
	bool cuAddComplexScalar(int numElements, std::complex<float>* d_A, std::complex<float> alpha);
};

