#ifndef __MY_TPS_H
#define __MY_TPS_H
#include "R3Shapes/R3Shapes.h"
#include "BasicNumerical/BasicNumerical.h"
#include <vector>
using namespace std;
class ThinPlateSpline3D
{
private:
	// site positions before deformation
//	vector<R3Point> m_sites_src;
	// site positions after deformation
//	vector<R3Point> m_sites_tar;
	// coefficients
//	vector<double> m_coeffs;
private:
	LinAlgVectorReal m_v_sites_src[3];
	LinAlgVectorReal m_v_sites_tar[3];
	LinAlgVectorReal m_v_coeffs[3];
	// control regularization
	RNScalar m_lambda;
public:
	// constructor
	ThinPlateSpline3D();
	ThinPlateSpline3D(vector<R3Point>& src, vector<R3Point>& tar, RNScalar lambda = 0.0);
	// assign source and target sites
	int AssignSites(vector<R3Point>& src, vector<R3Point>& tar);
	// get coefficients
	int ComputeCoefficients();
	// deformation
	int ConductDeformation(LinAlgVectorReal psrc[], LinAlgVectorReal ptar[]);
	int ConductDeformation(vector<R3Point>& psrc, vector<R3Point>& ptar);
	// det
	RNScalar Det();
};
#endif
