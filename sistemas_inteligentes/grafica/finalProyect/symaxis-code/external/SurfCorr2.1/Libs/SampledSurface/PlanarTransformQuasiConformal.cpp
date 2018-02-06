#include "PlanarTransformQuasiConformal.h"

PlanarTransformQuasiConformal::PlanarTransformQuasiConformal(FPIGeneralization nPntCase)
: L(2,2)
{
	m_iNumConstraints = 0;
	m_NPointCase = nPntCase;
	L.SetValues(0.);
	L(0,0) = 1.;
}

LinAlgComplex PlanarTransformQuasiConformal::Transform(LinAlgComplex z)
{
	LinAlgComplex w;
	if (m_iNumConstraints==3 || m_iNumConstraints==4)
	{
		w = mz.Transform(z);
		w = LinearTransform(w);
		w = mw.Transform(w);
	}
	return w;
}

LinAlgComplex PlanarTransformQuasiConformal::TransformInv(LinAlgComplex z)
{
	assert(false);
}

void PlanarTransformQuasiConformal::FindTransformation(std::vector<LinAlgComplex> &z, 
													   std::vector<LinAlgComplex> & w)
{
	if (z.size()!=w.size())
	{
		std::cout<<"[WARNING] Number of constraints is not equal"<<std::endl;
		return;
	}
	
	m_iNumConstraints = (int)z.size();
	
	if (m_iNumConstraints==3)	// conformal map
	{
		mz.FindTransformation(z, w);
		mw = MobiusTransformation();
		L.SetValues(0);
		L(0,0) = 1.;
	}
	else if (m_iNumConstraints==4)	// quasi-conformal (4-point deformation)
	{
		MapQuadrupletToParallelogram(z, mz);
		MapQuadrupletToParallelogram(w, mw);
		LinAlgVectorComplex tZ(m_iNumConstraints);
		LinAlgVectorComplex tW(m_iNumConstraints);	
		for (int i=0; i<m_iNumConstraints; i++)
		{
			tZ.SetValue(i, mz.Transform(z[i]));
			tW.SetValue(i, mw.Transform(w[i]));
		}
		FillAffineOn4Pnt(tZ, tW);
		mw = mw.Inverse();
		
//		std::cout<<"==================="<<std::endl;
//		std::cout<<"z=["<<z[0].r<<"+"<<z[0].i<<"i "<<z[1].r<<"+"<<z[1].i<<"i ";
//		std::cout<<z[2].r<<"+"<<z[2].i<<"i "<<z[3].r<<"+"<<z[3].i<<"i]"<<std::endl;
//		std::cout<<"w=["<<w[0].r<<"+"<<w[0].i<<"i "<<w[1].r<<"+"<<w[1].i<<"i ";
//		std::cout<<w[2].r<<"+"<<w[2].i<<"i "<<w[3].r<<"+"<<w[3].i<<"i]"<<std::endl;
//		
//		std::cout<<"mw=";	mw.Matrix().PrintMatlab();	std::cout<<std::endl;
//		std::cout<<"mz=";	mz.Matrix().PrintMatlab();	std::cout<<std::endl;
//		std::cout<<"L=";	L.PrintMatlab();	std::cout<<std::endl;
	}
}

void PlanarTransformQuasiConformal::SaveTransformation(std::ofstream & textStream)
{
	std::cout<<"[WARNING] Cannot save Quasi-conformal transform"<<std::endl;
}

void PlanarTransformQuasiConformal::LoadTransformation(std::ifstream & textStream)
{
	std::cout<<"[WARNING] Cannot laod Quasi-conformal transform"<<std::endl;
}


LinAlgComplex PlanarTransformQuasiConformal::LinearTransform(LinAlgComplex z)
{
	LinAlgComplex l1(L(0, 0), L(0, 1));
	LinAlgComplex l2(L(1, 0), L(1, 1));	
	
	LinAlgComplex o = l1*z + l2*z.Conjugate();;
	return l1*z + l2*z.Conjugate();	
}

void PlanarTransformQuasiConformal::MapQuadrupletToParallelogram(std::vector<LinAlgComplex> & z, 
																 MobiusTransformation & m)
{
	assert (z.size() == 4);
	
	LinAlgMatrixComplex A(2, 4);
	A.SetValue(0, 0, z[2]-z[0]);
	A.SetValue(0, 1, 0);
	A.SetValue(0, 2, z[0]-z[1]);
	A.SetValue(0, 3, z[1]-z[2]);
	
	A.SetValue(1, 0, z[3]-z[1]);
	A.SetValue(1, 1, z[1]-z[3]);
	A.SetValue(1, 2, z[0]-z[1]);
	A.SetValue(1, 3, z[1]-z[3]);	
	
	LinAlgMatrixComplex U(2, 2);
	LinAlgVectorComplex S(2);
	LinAlgMatrixComplex Vt(4, 4);
	A.SVDDecomposeDestructive(U, S, Vt);	
	
	LinAlgVectorComplex u(4);
	LinAlgVectorComplex w(4);	
	
	for (int i=0; i<4; i++)
	{
		u.SetValue(i, Vt.GetValue(2, i).Conjugate());
		w.SetValue(i, Vt.GetValue(3, i).Conjugate());
	}
	
	LinAlgComplex a2 = w.GetValue(2)*w.GetValue(1)-w.GetValue(3)*w.GetValue(0);
	LinAlgComplex a1 = w.GetValue(2)*u.GetValue(1)+u.GetValue(2)*w.GetValue(1)-w.GetValue(3)*u.GetValue(0)-u.GetValue(3)*w.GetValue(0);
	LinAlgComplex a0 = u.GetValue(2)*u.GetValue(1) - u.GetValue(3)*u.GetValue(0);
	
	LinAlgComplex r1 = (-a1 + (a1*a1 - 4.*a0*a2).Sqrt() ) / (2. * a2);
	LinAlgComplex r2 = (-a1 - (a1*a1 - 4.*a0*a2).Sqrt() ) / (2. * a2);	
	
	LinAlgVectorComplex x(4);
	for (int i=0; i<4; i++)
		x.SetValue(i, u.GetValue(i) + r2 * w.GetValue(i));
	
	LinAlgComplex tau;
	tau = x.GetValue(2) / x.GetValue(0);
	
	if (tau.i < 0)
	{
		for (int i=0; i<4; i++)
			x.SetValue(i, u.GetValue(i) + r1 * w.GetValue(i));		
		tau = x.GetValue(2) / x.GetValue(0);
	}
	m.a = x.GetValue(0);
	m.b = -m.a*z[0];
	m.c = x.GetValue(1);
	m.d = m.a * (z[1]-z[0]) - m.c*z[1];

}

void PlanarTransformQuasiConformal::FillAffineOn4Pnt(LinAlgVectorComplex & tZ, LinAlgVectorComplex & tW)
{
	int N = tZ.Rows();
	assert(tZ.Rows()==tW.Rows());
	assert (N==4);
		
	LinAlgComplex u(tZ(1)-tZ(0));
	LinAlgComplex v(tZ(2)-tZ(1));	
	
	LinAlgComplex U(tW(1)-tW(0));
	LinAlgComplex V(tW(2)-tW(1));	
	
	LinAlgMatrixComplex A(2,2);
	A.SetValue(0, 0, u);		A.SetValue(0, 1, u.Conjugate());
	A.SetValue(1, 0, v);		A.SetValue(1, 1, v.Conjugate());	
	
	LinAlgVectorComplex b(2);
	b.SetValue(0, U);
	b.SetValue(1, V);	
	
	A.Inverse();
	
	LinAlgVectorComplex l1l2(2);	
	A.Multiply(b, l1l2);
	
	L(0, 0) = l1l2.GetValue(0).r;
	L(0, 1) = l1l2.GetValue(0).i;
	L(1, 0) = l1l2.GetValue(1).r;
	L(1, 1) = l1l2.GetValue(1).i;	
}




