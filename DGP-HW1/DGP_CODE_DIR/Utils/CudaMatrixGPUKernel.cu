#ifndef _CUDA_MATRIX_GPU_KERNEL_
#define _CUDA_MATRIX_GPU_KERNEL_

//this function is limited by the amount of shared memory and should only be called when blockSize >= n
template<int blockSize>
__device__ void multMatrixByVectorComplexLimitedKernel(int m, int n, int numPitchedRows, const float2* d_A, const float2* d_X, float2* d_b)
{
	const int rowIndex = blockSize*blockIdx.x + threadIdx.x;

	__shared__ float2 X[blockSize];

	if(threadIdx.x <= n)
	{
		X[threadIdx.x] = d_X[threadIdx.x];
	}
	__syncthreads();

	if(rowIndex < m)
	{
		float2 sum = d_b[rowIndex];

		for(int i = 0; i < n; i++)
		{
			float2 complex1 = d_A[i*numPitchedRows + rowIndex];
			float2 complex2 = X[i];

			sum.x += complex1.x*complex2.x - complex1.y*complex2.y;
			sum.y += complex1.y*complex2.x + complex1.x*complex2.y;
		}

		d_b[rowIndex] = sum;
	}
}

//this function is limited by the amount of shared memory and should only be called when blockSize >= n
template<int blockSize>
__device__ void multMatrixByVectorRealLimitedKernel(int m, int n, int numPitchedRows, const float* d_A, const float* d_X, float* d_b)
{
	const int rowIndex = blockSize*blockIdx.x + threadIdx.x;

	__shared__ float X[blockSize];

	if(threadIdx.x <= n)
	{
		X[threadIdx.x] = d_X[threadIdx.x];
	}
	__syncthreads();

	if(rowIndex < m)
	{
		float sum = d_b[rowIndex];

		for(int i = 0; i < n; i++)
		{
			sum += X[i]*d_A[i*numPitchedRows + rowIndex];
		}

		d_b[rowIndex] = sum;
	}
}

template<int blockSize>
__global__ void multMatrixByVectorComplexKernel(int m, int n, int numPitchedRows, int numCycles, const float2* d_A, const float2* d_X, float2* d_b)
{
	int i = 0;
	for(i = 0; i < numCycles - 1; i++)
	{
		multMatrixByVectorComplexLimitedKernel<blockSize>(m, blockSize, numPitchedRows, d_A + i*numPitchedRows*blockSize, d_X + i*blockSize, d_b);
		__syncthreads();
	}
	int numColumnsLeft = n - (numCycles - 1)*blockSize;
	multMatrixByVectorComplexLimitedKernel<blockSize>(m, numColumnsLeft, numPitchedRows, d_A + i*numPitchedRows*blockSize, d_X + i*blockSize, d_b);
}

template<int blockSize>
__global__ void multMatrixByVectorRealKernel(int m, int n, int numPitchedRows, int numCycles, const float* d_A, const float* d_X, float* d_b)
{
	int i = 0;
	for(i = 0; i < numCycles - 1; i++)
	{
		multMatrixByVectorRealLimitedKernel<blockSize>(m, blockSize, numPitchedRows, d_A + i*numPitchedRows*blockSize, d_X + i*blockSize, d_b);
		__syncthreads();
	}
	int numColumnsLeft = n - (numCycles - 1)*blockSize;
	multMatrixByVectorRealLimitedKernel<blockSize>(m, numColumnsLeft, numPitchedRows, d_A + i*numPitchedRows*blockSize, d_X + i*blockSize, d_b);
}

template<int blockSize>
__global__ void addVectorsRealKernel(int numElements, const float* d_A, const float* d_B, float* d_C)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if(index < numElements)
	{
		d_C[index] = d_A[index] + d_B[index];
	}
}

template<int blockSize>
__global__ void addVectorsComplexKernel(int numElements, const float2* d_A, const float2* d_B, float2* d_C)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if(index < numElements)
	{
		float2 complexA = d_A[index];
		float2 complexB = d_B[index];
		float2 complexC;
		complexC.x = complexA.x + complexB.x;
		complexC.y = complexA.y + complexB.y;

		d_C[index] = complexC;
	}
}

template<int blockSize>
__global__ void subVectorsComplexKernel(int numElements, const float2* d_A, const float2* d_B, float2* d_C)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if(index < numElements)
	{
		float2 complexA = d_A[index];
		float2 complexB = d_B[index];
		float2 complexC;
		complexC.x = complexA.x - complexB.x;
		complexC.y = complexA.y - complexB.y;

		d_C[index] = complexC;
	}
}

template<int blockSize>
__global__ void subVectorsComplexDoubleKernel(int numElements, const double2* d_A, const double2* d_B, double2* d_C)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if (index < numElements)
	{
		double2 complexA = d_A[index];
		double2 complexB = d_B[index];
		double2 complexC;
		complexC.x = complexA.x - complexB.x;
		complexC.y = complexA.y - complexB.y;

		d_C[index] = complexC;
	}
}

template<int blockSize>
__global__ void subVectorsRealKernel(int numElements, const float* d_A, const float* d_B, float* d_C)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if(index < numElements)
	{
		d_C[index] = d_A[index] - d_B[index];
	}
}

template<int blockSize>
__global__ void scaleVectorComplexKernel(int numElements, float2* d_A, float2 alpha)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if(index < numElements)
	{
		float2 C = d_A[index];
		float temp = C.x;
		C.x = alpha.x*C.x - alpha.y*C.y;
		C.y = alpha.y*temp + alpha.x*C.y;

		d_A[index] = C;
	}
}

template<int blockSize>
__global__ void scaleVectorComplexDoubleKernel(int numElements, double2* d_A, double2 alpha)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if (index < numElements)
	{
		double2 C = d_A[index];
		double temp = C.x;
		C.x = alpha.x*C.x - alpha.y*C.y;
		C.y = alpha.y*temp + alpha.x*C.y;

		d_A[index] = C;
	}
}

template<int blockSize>
__global__ void exponentVectorComplexKernel(int numElements, float2* d_A)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if(index < numElements)
	{
		float2 C = d_A[index];
		float sinY = 0.0f;
		float cosY = 0.0f;
		sincosf(C.y, &sinY, &cosY);
		float expX = expf(C.x);
		C.x = expX*cosY;
		C.y = expX*sinY;

		d_A[index] = C;
	}
}


template<int blockSize>
__global__ void exponentVectorComplexDoubleKernel(int numElements, double2* d_A)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if (index < numElements)
	{
		double2 C = d_A[index];
		double sinY = 0.0f;
		double cosY = 0.0f;
		sincos(C.y, &sinY, &cosY);
		double expX = exp(C.x);
		C.x = expX*cosY;
		C.y = expX*sinY;

		d_A[index] = C;
	}
}

template<int blockSize>
__global__ void conjugateComplexKernel(int numElements, float2* d_A)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if(index < numElements)
	{
		float2 C = d_A[index];
		C.y = -C.y;

		d_A[index] = C;
	}
}

template<int blockSize>
__global__ void convertRealToComplexKernel(int numElements, float2* d_A, const float* d_B)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if(index < numElements)
	{
		float2 C;
		C.x = d_B[index];
		C.y = 0.0f;

		d_A[index] = C;
	}
}

template<int blockSize>
__global__ void addComplexScalarKernel(int numElements, float2* d_A, float2 alpha)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if(index < numElements)
	{
		float2 C = d_A[index];
		C.x += alpha.x;
		C.y += alpha.y;

		d_A[index] = C;
	}
}

template<int blockSize>
__global__ void addRealScalarKernel(int numElements, float* d_A, float alpha)
{
	const int index = blockSize*blockIdx.x + threadIdx.x;

	if(index < numElements)
	{
		d_A[index] += alpha;
	}
}


#endif //#ifndef _CUDA_MATRIX_GPU_KERNEL_
