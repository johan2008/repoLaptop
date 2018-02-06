#ifndef __MY_TPS_H
#define __MY_TPS_H
#include "R3Shapes/R3Shapes.h"
#include "BasicNumerical/BasicNumerical.h"
#include <vector>
using namespace std;
class ThinPlateSpline2D
{
private:
	// site positions before deformation
//	vector<R3Point> m_sites_src;
	// site positions after deformation
//	vector<R3Point> m_sites_tar;
	// coefficients
//	vector<double> m_coeffs;
private:
	LinAlgVectorReal m_v_sites_src[2];
	LinAlgVectorReal m_v_sites_tar[2];
	LinAlgVectorReal m_v_coeffs[2];
public:
	// constructor
	ThinPlateSpline2D();
	ThinPlateSpline2D(vector<R2Point>& src, vector<R2Point>& tar);
	// assign source and target sites
	int AssignSites(vector<R2Point>& src, vector<R2Point>& tar);
	// get coefficients
	int ComputeCoefficients();
	// deformation
	int ConductDeformation(LinAlgVectorReal psrc[], LinAlgVectorReal ptar[]);
	int ConductDeformation(vector<R2Point>& psrc, vector<R2Point>& ptar);
	// det
	RNScalar Det();
};
#endif
