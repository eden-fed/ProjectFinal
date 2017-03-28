#include "Timer.h"

#include <Windows.h>
#include <time.h>
#include <iostream>

using namespace std;

void get_time_str(double sec, char *st)
{
	int isec = int(sec);
	int milli = int((sec-isec)*1000) % 1000;

	time_t tsec = isec;
	tm now = *gmtime(&tsec); 
	char dest[100]={0};
	strftime(dest, sizeof(dest)-1, "%H:%M:%S", &now);
	sprintf(st, "%s.%03d", dest, milli);
}

CPUTimer::CPUTimer()
{
	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);
	mPerformanceFrequency = (double)(freq.QuadPart / 1000.0);
	tic();
}

void CPUTimer::tic()
{
	LARGE_INTEGER timer;
	QueryPerformanceCounter(&timer);
	mStartTime = (double)(timer.QuadPart);
}

double CPUTimer::toc()
{
	LARGE_INTEGER timer;
	QueryPerformanceCounter(&timer);
	double elapsedTime = ((double)(timer.QuadPart) - mStartTime) / mPerformanceFrequency;
	return elapsedTime;
}

double CPUTimer::print(std::ostream &os)
{
	double sec = toc()/1000.0;
	os << "CPUTimer reads ";
	char st[100];
	get_time_str(sec, st);
	os << st;	
	os.flush();
	os << " sec." << endl;
	return sec;
}

CUDATimer::CUDATimer()
{
	cudaEventCreate(&startEvent);
	cudaEventCreate(&stopEvent);
}

CUDATimer::~CUDATimer()
{
	cudaEventDestroy(startEvent);
	cudaEventDestroy(stopEvent);
}

void CUDATimer::startTimer()
{
	cudaEventCreate(&startEvent);
	cudaEventCreate(&stopEvent);
	cudaEventRecord(startEvent, 0);
}

float CUDATimer::stopTimer()
{
	float elapsedTime = 0.0f;

	cudaEventRecord(stopEvent, 0);
	cudaEventSynchronize(stopEvent);
	cudaEventElapsedTime(&elapsedTime, startEvent, stopEvent);

	return elapsedTime;
}

