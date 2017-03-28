#pragma once

#include "Utils/STL_Macros.h"


////////////////////////////////////////////
/////////        FastMatrix        /////////
///////// Copyright (C) Ofir Weber /////////
////////////////////////////////////////////


//forward declaration
template <class T> class MatrixGPU;
template <class T> class MatrixCPU;


typedef MatrixCPU<std::complex<float>> ComplexFloatCPUMatrix;
typedef MatrixGPU<std::complex<float>> ComplexFloatGPUMatrix;
typedef MatrixCPU<std::complex<double>> ComplexDoubleCPUMatrix;
typedef MatrixGPU<std::complex<double>> ComplexDoubleGPUMatrix;
typedef MatrixCPU<float> FloatCPUMatrix;
typedef MatrixGPU<float> FloatGPUMatrix;
typedef MatrixCPU<double> DoubleCPUMatrix;
typedef MatrixGPU<double> DoubleGPUMatrix;
