#include "math.h"

#ifndef _GPU_LOCAL_STEP_KERNEL_
#define _GPU_LOCAL_STEP_KERNEL_

#define MAX_NUMBER_OF_SEGMENTS 3

/*__shared__ bool allPointsInPolygon = true;
__global__ void CheckStopCondition(const int numElements, double2* log_fz, double2* nu_f, const double log_SigmaA, const double sigmaB, const double k, const double epsilon, const double m){
//use reduction?
}*/
__global__ void projectPointsToPolygonNokKernel(const int numElements, double2* log_fz, double2* nu_f, const double log_SigmaA, const double sigmaB, const double k, const double xIntersection, const double epsilon, const double m)
{
	const int index = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (index < numElements)
	{
		double2 log_fz_cur = log_fz[index];
		double2 nu_f_cur = nu_f[index];
		double x = sqrt(pow(nu_f_cur.x, 2) + pow(nu_f_cur.y, 2));
		double y = log_fz_cur.x;
		const double prevX = x;
		const double Yintersection = -log((1 - xIntersection) / sigmaB);

		if (!((x <= k + epsilon) && (x <= log_SigmaA - y + epsilon) && (log(sigmaB) + m*x <= y + epsilon))){
			if (y >= x + log_SigmaA){
				x = 0;
				y = log_SigmaA;
			}
			else if (y > x + Yintersection - xIntersection){
				x = (x - y + log_SigmaA) / 2;
				y = (y - prevX + log_SigmaA) / 2;
			}
			else if (y >= -x / m + Yintersection + xIntersection / m){
				x = xIntersection;
				y = Yintersection;
			}
			else if (y > -x / m + log(sigmaB)){
				x = (x + m*y - m*log(sigmaB)) / (pow(m, 2) + 1);
				y = m*((prevX + m*y - m*log(sigmaB)) / (pow(m, 2) + 1)) + log(sigmaB);
			}
			else{
				x = 0;
				y = log(sigmaB);
			}
			nu_f_cur.x = (x / prevX)*(nu_f_cur.x);
			nu_f_cur.y = (x / prevX)*(nu_f_cur.y);
			nu_f[index] = nu_f_cur;

			log_fz_cur.x = y;
			log_fz_cur.y = log_fz_cur.y;
			log_fz[index] = log_fz_cur;
		}
		//else - update if point is in polygon somehow
	}

}

__global__ void projectPointsToPolygonNokKernel_unsplit(const int numElements, double2* x_vec, const double log_SigmaA, const double sigmaB, const double k, const double xIntersection, const double epsilon, const double m)
{
	const int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < numElements)
	{
		double2 log_fz_cur = x_vec[index];
		double2 nu_f_cur = x_vec[index + numElements];
		double x = sqrt(pow(nu_f_cur.x, 2) + pow(nu_f_cur.y, 2));
		double y = log_fz_cur.x;
		const double prevX = x;
		const double Yintersection = -log((1 - xIntersection) / sigmaB);

		if (!((x <= k + epsilon) && (x <= log_SigmaA - y + epsilon) && (log(sigmaB) + m*x <= y + epsilon))){
			if (y >= x + log_SigmaA){
				x = 0;
				y = log_SigmaA;
			}
			else if (y > x + Yintersection - xIntersection){
				x = (x - y + log_SigmaA) / 2;
				y = (y - prevX + log_SigmaA) / 2;
			}
			else if (y >= -x / m + Yintersection + xIntersection / m){
				x = xIntersection;
				y = Yintersection;
			}
			else if (y > -x / m + log(sigmaB)){
				x = (x + m*y - m*log(sigmaB)) / (pow(m, 2) + 1);
				y = m*((prevX + m*y - m*log(sigmaB)) / (pow(m, 2) + 1)) + log(sigmaB);
			}
			else{
				x = 0;
				y = log(sigmaB);
			}
			nu_f_cur.x = (x / prevX)*(nu_f_cur.x);
			nu_f_cur.y = (x / prevX)*(nu_f_cur.y);
			x_vec[index + numElements] = nu_f_cur;

			log_fz_cur.x = y;
			log_fz_cur.y = log_fz_cur.y;
			x_vec[index] = log_fz_cur;
		}
		//else - update if point is in polygon somehow
	}

}

__global__ void projectPointsToPolygonWithkKernel(const int numElements, double2* log_fz, double2* nu_f, const double log_SigmaA, const double sigmaB, const double k, const double epsilon, const double m)
{
	const int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < numElements)
	{
		double2 log_fz_cur = log_fz[index];
		double2 nu_f_cur = nu_f[index];

		double x = sqrt(pow(nu_f_cur.x, 2) + pow(nu_f_cur.y, 2));
		double y = log_fz_cur.x;
		const double prevX = x;

		if (!((x <= k + epsilon) && (x <= log_SigmaA - y + epsilon) && (log(sigmaB) + m*x <= y + epsilon))){
			if (y >= x + log_SigmaA){
				x = 0;
				y = log_SigmaA;
			}
			else if (y > x + log_SigmaA - 2 * k){
				x = (x - y + log_SigmaA) / 2;
				y = (y - prevX + log_SigmaA) / 2;
			}
			else if (y >= log_SigmaA - k){
				x = k;
				y = log_SigmaA;
			}
			else if (y > log(sigmaB / (1 - k))){
				x = k;
			}
			else if (y >= -x / m + log(sigmaB / (1 - k)) + k / m){
				x = k;
				y = log(sigmaB / (1 - k));
			}
			else if (y > -x / m + log(sigmaB)){
				double prevX = x;
				x = (x + m*y - m*log(sigmaB)) / (pow(m, 2) + 1);
				y = m*((prevX + m*y - m*log(sigmaB)) / (pow(m, 2) + 1)) + log(sigmaB);
			}
			else{
				x = 0;
				y = log(sigmaB);
			}
			nu_f_cur.x = (x / prevX)*(nu_f_cur.x);
			nu_f_cur.y = (x / prevX)*(nu_f_cur.y);
			nu_f[index] = nu_f_cur;

			log_fz_cur.x = y;
			log_fz_cur.y = log_fz_cur.y;
			log_fz[index] = log_fz_cur;
		}
		//else - update if point is in polygon somehow
	}
}

__global__ void projectPointsToPolygonWithkKernel_unsplit(const int numElements, double2* x_vec, const double log_SigmaA, const double sigmaB, const double k, const double epsilon, const double m)
{
	const int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < numElements)
	{
		double2 log_fz_cur = x_vec[index];
		double2 nu_f_cur = x_vec[index + numElements];

		double x = sqrt(pow(nu_f_cur.x, 2) + pow(nu_f_cur.y, 2));
		double y = log_fz_cur.x;
		const double prevX = x;

		if (!((x <= k + epsilon) && (x <= log_SigmaA - y + epsilon) && (log(sigmaB) + m*x <= y + epsilon))){
			if (y >= x + log_SigmaA){
				x = 0;
				y = log_SigmaA;
			}
			else if (y > x + log_SigmaA - 2 * k){
				x = (x - y + log_SigmaA) / 2;
				y = (y - prevX + log_SigmaA) / 2;
			}
			else if (y >= log_SigmaA - k){
				x = k;
				y = log_SigmaA;
			}
			else if (y > log(sigmaB / (1 - k))){
				x = k;
			}
			else if (y >= -x / m + log(sigmaB / (1 - k)) + k / m){
				x = k;
				y = log(sigmaB / (1 - k));
			}
			else if (y > -x / m + log(sigmaB)){
				double prevX = x;
				x = (x + m*y - m*log(sigmaB)) / (pow(m, 2) + 1);
				y = m*((prevX + m*y - m*log(sigmaB)) / (pow(m, 2) + 1)) + log(sigmaB);
			}
			else{
				x = 0;
				y = log(sigmaB);
			}
			nu_f_cur.x = (x / prevX)*(nu_f_cur.x);
			nu_f_cur.y = (x / prevX)*(nu_f_cur.y);
			x_vec[index + numElements] = nu_f_cur;

			log_fz_cur.x = y;
			log_fz_cur.y = log_fz_cur.y;
			x_vec[index] = log_fz_cur;
		}
		//else - update if point is in polygon somehow
	}
}

__global__ void projectPointToPolygonMinSegKernel(const int numElements, double2* log_fz, double2* nu_f, const double* mXvaluesOfIntersections, const double* mYvaluesOfIntersections, const int NumOfsegments, const int epsilon){

	const int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < numElements)
	{
		double2 log_fz_cur = log_fz[index];
		double2 nu_f_cur = nu_f[index];

		double x = sqrt(pow(nu_f_cur.x, 2) + pow(nu_f_cur.y, 2));
		double y = log_fz_cur.x;
		const double prevX = x;
		bool isInsidePolygom = true;
		int i;

		for (i = 0; i < NumOfsegments; i++){
			double crossSegValue = (mXvaluesOfIntersections[i + 1] - mXvaluesOfIntersections[i])*(y - mYvaluesOfIntersections[i]) - (mYvaluesOfIntersections[i + 1] - mYvaluesOfIntersections[i])*(x - mXvaluesOfIntersections[i]);
			if (crossSegValue > 0 - epsilon){
				isInsidePolygom = false;
			}
		}

		if (!isInsidePolygom){

			double closestPointsXvalues[MAX_NUMBER_OF_SEGMENTS];
			double closestPointsYvalues[MAX_NUMBER_OF_SEGMENTS];
			double minDistances[MAX_NUMBER_OF_SEGMENTS];

			//using cross product
			for (i = 0; i < NumOfsegments; i++){
				double dot = (x - mXvaluesOfIntersections[i])*(mXvaluesOfIntersections[i + 1] - mXvaluesOfIntersections[i]) + (y - mYvaluesOfIntersections[i])*(mYvaluesOfIntersections[i + 1] - mYvaluesOfIntersections[i]);
				double projectionOnLine = dot / (pow((mXvaluesOfIntersections[i + 1] - mXvaluesOfIntersections[i]), 2) + pow((mYvaluesOfIntersections[i + 1] - mYvaluesOfIntersections[i]), 2));
				double t = fmax(0.0, fmin(1.0, projectionOnLine));
				double closestPointXvalue = mXvaluesOfIntersections[i] + t*(mXvaluesOfIntersections[i + 1] - mXvaluesOfIntersections[i]);
				double closestPointYvalue = mYvaluesOfIntersections[i] + t*(mYvaluesOfIntersections[i + 1] - mYvaluesOfIntersections[i]);
				closestPointsXvalues[i] = closestPointXvalue;
				closestPointsYvalues[i] = closestPointYvalue;
				minDistances[i] = pow((x - closestPointXvalue), 2) + pow((y - closestPointYvalue), 2);
			}
			int indexOfMin = 0;
			for (i = 1; i < NumOfsegments; i++)
				indexOfMin = (minDistances[i] < minDistances[indexOfMin] ? i : indexOfMin);

			x = closestPointsXvalues[indexOfMin];
			y = closestPointsYvalues[indexOfMin];

			nu_f_cur.x = (x / prevX)*(nu_f_cur.x);
			nu_f_cur.y = (x / prevX)*(nu_f_cur.y);
			nu_f[index] = nu_f_cur;

			log_fz_cur.x = y;
			log_fz_cur.y = log_fz_cur.y;
			log_fz[index] = log_fz_cur;

		}//else - update if point is in polygon somehow
	}
}

__global__ void projectPointToPolygonMinSegKernel_unsplit(const int numElements, double2* x_vec, const double* mXvaluesOfIntersections, const double* mYvaluesOfIntersections, const int NumOfsegments, const int epsilon){

	const int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < numElements)
	{
		double2 log_fz_cur = x_vec[index];
		double2 nu_f_cur = x_vec[index + numElements];

		double x = sqrt(pow(nu_f_cur.x, 2) + pow(nu_f_cur.y, 2));
		double y = log_fz_cur.x;
		const double prevX = x;
		bool isInsidePolygom = true;
		int i;

		for (i = 0; i < NumOfsegments; i++){
			double crossSegValue = (mXvaluesOfIntersections[i + 1] - mXvaluesOfIntersections[i])*(y - mYvaluesOfIntersections[i]) - (mYvaluesOfIntersections[i + 1] - mYvaluesOfIntersections[i])*(x - mXvaluesOfIntersections[i]);
			if (crossSegValue > 0 - epsilon){
				isInsidePolygom = false;
			}
		}

		if (!isInsidePolygom){

			double closestPointsXvalues[MAX_NUMBER_OF_SEGMENTS];
			double closestPointsYvalues[MAX_NUMBER_OF_SEGMENTS];
			double minDistances[MAX_NUMBER_OF_SEGMENTS];

			//using cross product
			for (i = 0; i < NumOfsegments; i++){
				double dot = (x - mXvaluesOfIntersections[i])*(mXvaluesOfIntersections[i + 1] - mXvaluesOfIntersections[i]) + (y - mYvaluesOfIntersections[i])*(mYvaluesOfIntersections[i + 1] - mYvaluesOfIntersections[i]);
				double projectionOnLine = dot / (pow((mXvaluesOfIntersections[i + 1] - mXvaluesOfIntersections[i]), 2) + pow((mYvaluesOfIntersections[i + 1] - mYvaluesOfIntersections[i]), 2));
				double t = fmax(0.0, fmin(1.0, projectionOnLine));
				double closestPointXvalue = mXvaluesOfIntersections[i] + t*(mXvaluesOfIntersections[i + 1] - mXvaluesOfIntersections[i]);
				double closestPointYvalue = mYvaluesOfIntersections[i] + t*(mYvaluesOfIntersections[i + 1] - mYvaluesOfIntersections[i]);
				closestPointsXvalues[i] = closestPointXvalue;
				closestPointsYvalues[i] = closestPointYvalue;
				minDistances[i] = pow((x - closestPointXvalue), 2) + pow((y - closestPointYvalue), 2);
			}
			int indexOfMin = 0;
			for (i = 1; i < NumOfsegments; i++)
				indexOfMin = (minDistances[i] < minDistances[indexOfMin] ? i : indexOfMin);

			x = closestPointsXvalues[indexOfMin];
			y = closestPointsYvalues[indexOfMin];

			nu_f_cur.x = (x / prevX)*(nu_f_cur.x);
			nu_f_cur.y = (x / prevX)*(nu_f_cur.y);
			x_vec[index + numElements] = nu_f_cur;

			log_fz_cur.x = y;
			log_fz_cur.y = log_fz_cur.y;
			x_vec[index] = log_fz_cur;

		}//else - update if point is in polygon somehow
	}
}

#endif