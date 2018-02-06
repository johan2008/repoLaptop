#include "tps3.h"
// default constructor
ThinPlateSpline3D::ThinPlateSpline3D():m_lambda(0.0)
{}
// constructor : assigning source and target sites
ThinPlateSpline3D::ThinPlateSpline3D(vector<R3Point>& src, vector<R3Point>& tar, RNScalar lambda):m_lambda(lambda)
{
	int suc = AssignSites(src, tar);
	if (!suc)
	{
		fprintf(stderr, "Fail in assigning source sites and target sites!\n");
	}
}

// assign source and target sites
int ThinPlateSpline3D::AssignSites(vector<R3Point>& src, vector<R3Point>& tar)
{
	if (src.size() != tar.size())
	{
		fprintf(stderr, "The size of source sites and target sites do not match!\n");
		return 0;
	}
	vector<double> in_src[3];
	vector<double> in_tar[3];
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<int(src.size()); j++)
		{
			in_src[i].push_back(src[j][i]);
			in_tar[i].push_back(tar[j][i]);
		}
	}
	// put it in m_v_sites_src and m_v_sites_tar
	for (int i=0; i<3; i++)
	{
		m_v_sites_src[i] = LinAlgVectorReal(in_src[i]);
		m_v_sites_tar[i] = LinAlgVectorReal(in_tar[i]);
	}
	return 1;
}

// compute coefficients
int ThinPlateSpline3D::ComputeCoefficients()
{
	int nsites = m_v_sites_src[0].Rows();
	LinAlgMatrixReal matrix(nsites+4, nsites+4);
	// compute alpha, for controlling regularization
	RNScalar alpha = 0;	
	// filling matrix
	// rbf of tps: phi = r*r*log r
	for (int i=0; i<nsites; i++)
	{
		double x[2], y[2], z[2];
		x[0] = m_v_sites_src[0](i);
		y[0] = m_v_sites_src[1](i);
		z[0] = m_v_sites_src[2](i);
		for (int j=0; j<nsites; j++)
		{
			x[1] = m_v_sites_src[0](j);
			y[1] = m_v_sites_src[1](j);
			z[1] = m_v_sites_src[2](j);
			// compute the euclidean distance between two sites
			double r = (x[0] - x[1]) * (x[0] - x[1])
					 + (y[0] - y[1]) * (y[0] - y[1])
					 + (z[0] - z[1]) * (z[0] - z[1]);
			alpha += sqrt(r);
			if (r == 0)
			{
				matrix(i, j) = 0;
			}
			else
			{
				matrix(i, j) = r * log(r)*0.5;
			}
		}
	}
	alpha = alpha/nsites/nsites;
	alpha = alpha * alpha;
	printf("Disturbing for regularization = %.3f\n", alpha*m_lambda);
	for (int i=0; i<nsites; i++)
	{
		matrix(i, i) = alpha*m_lambda;
	}
	// filling the rest of matrix
	for (int i=0; i<nsites; i++)
	{
		for (int j=0; j<3; j++)
		{
			matrix(nsites+j, i) = m_v_sites_src[j](i);
		}
		matrix(nsites+3, i) = 1.0;
		for (int j=0; j<3; j++)
		{
			matrix(i, nsites+j) = m_v_sites_src[j](i);
		}
		matrix(i, nsites+3) = 1.0;
	}

	// filling target vectors
	LinAlgVectorReal b[3];
	b[0] = LinAlgVectorReal(nsites+4, 0);
	b[1] = LinAlgVectorReal(nsites+4, 0);
	b[2] = LinAlgVectorReal(nsites+4, 0);
	for (int i=0; i<nsites; i++)
	{
		b[0](i) = m_v_sites_tar[0](i);
		b[1](i) = m_v_sites_tar[1](i);
		b[2](i) = m_v_sites_tar[2](i);
	}

	// initialize x
	m_v_coeffs[0] = LinAlgVectorReal(nsites+4, 0);
	m_v_coeffs[1] = LinAlgVectorReal(nsites+4, 0);
	m_v_coeffs[2] = LinAlgVectorReal(nsites+4, 0);

	printf("Before solving ...\n");
//	getchar();
	// solve for x
	for (int i=0; i<3; i++)
	{
		printf("Solving %d ...\n", i);
		LinAlgMatrixReal _mat(matrix);
		_mat.Solve(b[i], m_v_coeffs[i]);
	}
	printf("after solving ...\n");

	return 1;
}

int ThinPlateSpline3D::ConductDeformation(LinAlgVectorReal psrc[], LinAlgVectorReal ptar[])
{
	int npoints = psrc[0].Rows();
	for (int i=0; i<3; i++)
	{
		ptar[i] = LinAlgVectorReal(npoints, 0);
	}
	int nsites = m_v_sites_src[0].Rows();
	for (int i=0; i<npoints; i++)
	{
		double x[2], y[2], z[2];
		x[0] = psrc[0](i);
		y[0] = psrc[1](i);
		z[0] = psrc[2](i);
		for (int j=0; j<nsites; j++)
		{
			x[1] = m_v_sites_src[0](j);
			y[1] = m_v_sites_src[1](j);
			z[1] = m_v_sites_src[2](j);
			double r = (x[0] - x[1]) * (x[0] - x[1])
					 + (y[0] - y[1]) * (y[0] - y[1])
					 + (z[0] - z[1]) * (z[0] - z[1]);
			double di = 0;
			if (r == 0)
			{
				di = 0;
			}
			else
			{
				di = r * log(r)*0.5;
			}
			ptar[0](i) = ptar[0](i) + di * m_v_coeffs[0](j);
			ptar[1](i) = ptar[1](i) + di * m_v_coeffs[1](j);
			ptar[2](i) = ptar[2](i) + di * m_v_coeffs[2](j);
		}
		for (int j=0; j<3; j++)
		{
			for (int k=0; k<3; k++)
			{
				ptar[j](i) += psrc[k](i)*m_v_coeffs[j](nsites+k);
			}
			ptar[j](i) += m_v_coeffs[j](nsites+3);
		}
	}
	return 1;
}

int ThinPlateSpline3D::ConductDeformation(vector<R3Point>& psrc, vector<R3Point>& ptar)
{
	int npoints = psrc.size();
	LinAlgVectorReal vpts[3];
	vpts[0] = LinAlgVectorReal(npoints);
	vpts[1] = LinAlgVectorReal(npoints);
	vpts[2] = LinAlgVectorReal(npoints);
	for (int i=0; i<npoints; i++)
	{
		R3Point p = psrc[i];
		vpts[0](i) = p[0];
		vpts[1](i) = p[1];
		vpts[2](i) = p[2];
	}
	LinAlgVectorReal _vpts[3];
	int suc = ConductDeformation(vpts, _vpts);
	if (!suc)
	{
		fprintf(stderr, "Fail to conduct deformation!\n");
		return 0;
	}
	ptar.clear();
	for (int i=0; i<npoints; i++)
	{
		R3Point p(_vpts[0](i), _vpts[1](i), _vpts[2](i));
		ptar.push_back(p);
	}
	return 1;
}

RNScalar ThinPlateSpline3D::Det()
{
	int nsites = m_v_sites_src[0].Rows();
	RNScalar a[3][3];
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			a[i][j] = m_v_coeffs[i](nsites+j);
		}
	}
	RNScalar det = a[0][0]*(a[1][1]*a[2][2] - a[1][2]*a[2][1])
				 - a[0][1]*(a[1][0]*a[2][2] - a[2][0]*a[1][2])
				 + a[0][2]*(a[1][0]*a[2][1] - a[2][0]*a[1][1]);
	return det;
}
