#pragma once

#include <cuda_runtime_api.h>
#include <iostream>

class CPUTimer
{
public:
	CPUTimer();

public:
	void tic();
	double toc();
    double print(std::ostream &os = std::cout);

protected:
	double mPerformanceFrequency;
	double mStartTime;
};

class CUDATimer
{
public:
	CUDATimer();
	~CUDATimer();

public:
	void startTimer();
	float stopTimer();

protected:
	cudaEvent_t startEvent;
	cudaEvent_t stopEvent;
};
