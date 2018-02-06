#include "tps2.h"
// default constructor
ThinPlateSpline2D::ThinPlateSpline2D()
{}
// constructor : assigning source and target sites
ThinPlateSpline2D::ThinPlateSpline2D(vector<R2Point>& src, vector<R2Point>& tar)
{
	int suc = AssignSites(src, tar);
	if (!suc)
	{
		fprintf(stderr, "Fail in assigning source sites and target sites!\n");
	}
}

// assign source and target sites
int ThinPlateSpline2D::AssignSites(vector<R2Point>& src, vector<R2Point>& tar)
{
	if (src.size() != tar.size())
	{
		fprintf(stderr, "The size of source sites and target sites do not match!\n");
		return 0;
	}
	vector<double> in_src[2];
	vector<double> in_tar[2];
	for (int i=0; i<2; i++)
	{
		for (int j=0; j<int(src.size()); j++)
		{
			in_src[i].push_back(src[j][i]);
			in_tar[i].push_back(tar[j][i]);
		}
	}
	// put it in m_v_sites_src and m_v_sites_tar
	for (int i=0; i<2; i++)
	{
		m_v_sites_src[i] = LinAlgVectorReal(in_src[i]);
		m_v_sites_tar[i] = LinAlgVectorReal(in_tar[i]);
	}
	return 1;
}

// compute coefficients
int ThinPlateSpline2D::ComputeCoefficients()
{
	int nsites = m_v_sites_src[0].Rows();
	LinAlgMatrixReal matrix(nsites+3, nsites+3);
	// filling matrix
	// rbf of tps: phi = r*r*log r
	for (int i=0; i<nsites; i++)
	{
		double x[2], y[2];
		x[0] = m_v_sites_src[0](i);
		y[0] = m_v_sites_src[1](i);
		for (int j=0; j<nsites; j++)
		{
			x[1] = m_v_sites_src[0](j);
			y[1] = m_v_sites_src[1](j);
			// compute the euclidean distance between two sites
			double r = (x[0] - x[1]) * (x[0] - x[1])
					 + (y[0] - y[1]) * (y[0] - y[1]);

			if (r == 0)
			{
				matrix(i, j) = 0;
			}
			else
			{
				matrix(i, j) = r * log(r);
			}
		}
//		matrix(i, i) += RN_EPSILON;
	}
	// filling the rest of matrix
	for (int i=0; i<nsites; i++)
	{
		for (int j=0; j<2; j++)
		{
			matrix(nsites+j, i) = m_v_sites_src[j](i);
		}
		matrix(nsites+2, i) = 1.0;
		for (int j=0; j<2; j++)
		{
			matrix(i, nsites+j) = m_v_sites_src[j](i);
		}
		matrix(i, nsites+2) = 1.0;
	}

	// filling target vectors
	LinAlgVectorReal b[2];
	b[0] = LinAlgVectorReal(nsites+3, 0);
	b[1] = LinAlgVectorReal(nsites+3, 0);
	for (int i=0; i<nsites; i++)
	{
		b[0](i) = m_v_sites_tar[0](i);
		b[1](i) = m_v_sites_tar[1](i);
	}

	// initialize x
	m_v_coeffs[0] = LinAlgVectorReal(nsites+3, 0);
	m_v_coeffs[1] = LinAlgVectorReal(nsites+3, 0);

	printf("Before solving ...\n");
//	getchar();
	// solve for x
	for (int i=0; i<2; i++)
	{
		printf("Solving %d ...\n", i);
		LinAlgMatrixReal _mat(matrix);
   /*     printf("B: [ ");*/
		/*for (int j=0; j<b[i].Rows(); j++)*/
		/*{*/
			/*printf("%3.3f    ", b[i](j));*/
		/*}*/
		/*printf(" ]\n");*/

		/*printf("M: \n");*/
		/*for (int j=0; j<_mat.Rows(); j++)*/
		/*{*/
			/*for (int k=0; k<_mat.Cols(); k++)*/
			/*{*/
				/*printf("%3.3f    ", _mat(j, k));*/
			/*}*/
			/*printf("\n");*/
		/*}*/
		/*printf("\n");*/
//		getchar();
		_mat.Solve(b[i], m_v_coeffs[i]);
	}
	printf("after solving ...\n");

	return 1;
}

int ThinPlateSpline2D::ConductDeformation(LinAlgVectorReal psrc[], LinAlgVectorReal ptar[])
{
	int npoints = psrc[0].Rows();
	for (int i=0; i<2; i++)
	{
		ptar[i] = LinAlgVectorReal(npoints, 0);
	}
	int nsites = m_v_sites_src[0].Rows();
	for (int i=0; i<npoints; i++)
	{
		double x[2], y[2];
		x[0] = psrc[0](i);
		y[0] = psrc[1](i);
		for (int j=0; j<nsites; j++)
		{
			x[1] = m_v_sites_src[0](j);
			y[1] = m_v_sites_src[1](j);
			double r = (x[0] - x[1]) * (x[0] - x[1])
					 + (y[0] - y[1]) * (y[0] - y[1]);
			double di = 0;
			if (r == 0)
			{
				di = 0;
			}
			else
			{
				di = r * log(r);
			}
			ptar[0](i) = ptar[0](i) + di * m_v_coeffs[0](j);
			ptar[1](i) = ptar[1](i) + di * m_v_coeffs[1](j);
		}
		for (int j=0; j<2; j++)
		{
			for (int k=0; k<2; k++)
			{
				ptar[j](i) += psrc[k](i)*m_v_coeffs[j](nsites+k);
			}
			ptar[j](i) += m_v_coeffs[j](nsites+2);
		}
	}
	return 1;
}

int ThinPlateSpline2D::ConductDeformation(vector<R2Point>& psrc, vector<R2Point>& ptar)
{
	int npoints = psrc.size();
	LinAlgVectorReal vpts[2];
	vpts[0] = LinAlgVectorReal(npoints);
	vpts[1] = LinAlgVectorReal(npoints);
	for (int i=0; i<npoints; i++)
	{
		R2Point p = psrc[i];
		vpts[0](i) = p[0];
		vpts[1](i) = p[1];
	}
	LinAlgVectorReal _vpts[2];
	int suc = ConductDeformation(vpts, _vpts);
	if (!suc)
	{
		fprintf(stderr, "Fail to conduct deformation!\n");
		return 0;
	}
	ptar.clear();
	for (int i=0; i<npoints; i++)
	{
		R2Point p(_vpts[0](i), _vpts[1](i));
		ptar.push_back(p);
	}
	return 1;
}

RNScalar ThinPlateSpline2D::Det()
{
	int nsites = m_v_sites_src[0].Rows();
	RNScalar a = m_v_coeffs[0](nsites+0);
	RNScalar b = m_v_coeffs[0](nsites+1);
	RNScalar c = m_v_coeffs[1](nsites+0);
	RNScalar d = m_v_coeffs[1](nsites+1);
	return a*d - b*c;
}
