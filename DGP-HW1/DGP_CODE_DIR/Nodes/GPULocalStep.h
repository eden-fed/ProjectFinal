#pragma once

#include <complex>

#define I_DIV_UP(a, b) ((((a) % (b)) != 0) ? ((a) / (b) + 1) : ((a) / (b)))

extern "C"
{
	bool cuProjectPointsToPolygonNoK(int numElements, std::complex<double>* log_fz, std::complex<double>* nu_f, const double log_SigmaA, const double sigmaB, const double k, const double xIntersection, const double epsilon, const double m);
	bool cuProjectPointsToPolygonNoK_HP(int numElements, std::complex<double>* x_vec, const double log_SigmaA, const double sigmaB, const double k, const double xIntersection, const double epsilon, const double m);
	bool cuProjectPointsToPolygonWithK(int numElements, std::complex<double>* log_fz, std::complex<double>* nu_f, const double log_SigmaA, const double sigmaB, const double k, const double epsilon, const double m);
	bool cuProjectPointsToPolygonWithK_HP(int numElements, std::complex<double>* x_vec, const double log_SigmaA, const double sigmaB, const double k, const double epsilon, const double m);
	bool cuProjectPointToPolygonMinSeg(int numElements, std::complex<double>* log_fz, std::complex<double>* nu_f, const double* mXvaluesOfIntersections, const double* mYvaluesOfIntersections, const int NumOfsegments, const int epsilon);
	bool cuProjectPointToPolygonMinSeg_HP(int numElements, std::complex<double>* x_vec, const double* mXvaluesOfIntersections, const double* mYvaluesOfIntersections, const int NumOfsegments, const int epsilon);

}
