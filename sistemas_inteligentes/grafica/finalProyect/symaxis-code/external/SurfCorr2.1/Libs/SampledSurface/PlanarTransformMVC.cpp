#include "PlanarTransformMVC.h"

PlanarTransformMVC::PlanarTransformMVC()
{
}

LinAlgComplex PlanarTransformMVC::Transform(LinAlgComplex z)
{
	std::vector<double> mvc;
	CalcMVC(m_Z, z.r, z.i, mvc);
	LinAlgComplex tZ;
	Interpolate(m_W, mvc, tZ.r, tZ.i);	
	return tZ;
}

LinAlgComplex PlanarTransformMVC::TransformInv(LinAlgComplex z)
{
	assert(false);
}

void PlanarTransformMVC::FindTransformation(std::vector<LinAlgComplex> &z, 
						std::vector<LinAlgComplex> & w)
{
	for (int i=0; i<(int)z.size(); i++)
	{
		m_Z.push_back(z[i].r);
		m_Z.push_back(z[i].i);
	}

	for (int i=0; i<(int)w.size(); i++)
	{
		m_W.push_back(w[i].r);
		m_W.push_back(w[i].i);
	}
}

void PlanarTransformMVC::SaveTransformation(std::ofstream & textStream)
{
	std::cout<<"[WARNING] Saving PlanarTransformMVC not implemented!"<<std::endl;
}

void PlanarTransformMVC::LoadTransformation(std::ifstream & textStream)
{
	std::cout<<"[WARNING] Loading PlanarTransformMVC not implemented!"<<std::endl;
}

void PlanarTransformMVC::CalcMVC(std::vector<double> & ctrlPnts, 
									double x, double y, 
									std::vector<double> & MVC)
{
	int N = (int)ctrlPnts.size()/2;
	double x1, x2, y1, y2;	
	std::vector<double> tanHalfAlpha;
	std::vector<double> leni;	
	MVC.clear();
	for (int i=0; i<N; i++)
	{
		int in = (i+1)%N;
		x1 = ctrlPnts[2*i + 0] - x;
		y1 = ctrlPnts[2*i + 1] - y;		
		x2 = ctrlPnts[2*in + 0] - x;
		y2 = ctrlPnts[2*in + 1] - y;
		
		double len1 = sqrt(x1*x1 + y1*y1);
		double len2 = sqrt(x2*x2 + y2*y2);
		
		if (len1==0 || len2==0) 
		{
			MVC.clear();
			for (int j=0; j<N; j++)
				MVC.push_back(0.);
			MVC[((len1==0) ? i : in)] = 1.;
			return;
		}
		
		
		y2/=len2;		x2/=len2;
		y1/=len1;		x1/=len1;
		double alphai = atan2(y2, x2) - atan2(y1, x1);
		tanHalfAlpha.push_back(tan(alphai/2.));
		leni.push_back(len1);
	}
	
	double norm=0;
	for (int i=0; i<N; i++)
	{
		int ip = (i-1+N)%N;
		double wi = (tanHalfAlpha[ip] + tanHalfAlpha[i]) / leni[i];
		MVC.push_back(wi);
		norm += wi;
	}
	
	for (int i=0; i<N; i++)
		MVC[i] /= norm;
}


void PlanarTransformMVC::Interpolate(std::vector<double> & ctrlPnts, 
										std::vector<double> & MVC,
										double & newX, double &newY)
{
	assert(MVC.size()==ctrlPnts.size()/2);
	assert(ctrlPnts.size()%2==0);
	newX = 0;
	newY = 0;
	
	for (int i=0; i<(int)MVC.size(); i++)
	{
		newX += ctrlPnts[2*i+0] * MVC[i];
		newY += ctrlPnts[2*i+1] * MVC[i];
	}
}


