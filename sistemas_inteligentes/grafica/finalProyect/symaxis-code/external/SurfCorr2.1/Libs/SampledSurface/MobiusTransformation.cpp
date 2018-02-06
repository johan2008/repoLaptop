#include "VKString.h"
#include "MobiusTransformation.h"
#include <fstream>

////////////////////// MobiusTransformation /////////////////////

MobiusTransformation::MobiusTransformation()
{
	a = 1;
	b = 0;
	c = 0;
	d = 1;
}

MobiusTransformation::MobiusTransformation(const LinAlgMatrixComplex & M)
{
	assert(M.Rows()==2 && M.Cols()==2);
	a = M.GetValue(0, 0);
	b = M.GetValue(0, 1);
	c = M.GetValue(1, 0);
	d = M.GetValue(1, 1);
}


MobiusTransformation::MobiusTransformation(LinAlgComplex a, LinAlgComplex b, LinAlgComplex c, LinAlgComplex d)
{
	this->a = a;
	this->b = b;
	this->c = c;
	this->d = d;
}
		
void MobiusTransformation::SaveTransformation(std::ofstream & textStream)
{
	textStream<<"MobiusTransformationEncoded:\n";
	int sizeDouble = sizeof(double);
	char * numString = new char[sizeDouble * 8];
	
	SaveDouble(a.r, &(numString[sizeDouble*0]));	
	SaveDouble(a.i, &(numString[sizeDouble*1]));	
	SaveDouble(b.r, &(numString[sizeDouble*2]));	
	SaveDouble(b.i, &(numString[sizeDouble*3]));	
	SaveDouble(c.r, &(numString[sizeDouble*4]));	
	SaveDouble(c.i, &(numString[sizeDouble*5]));	
	SaveDouble(d.r, &(numString[sizeDouble*6]));	
	SaveDouble(d.i, &(numString[sizeDouble*7]));
	textStream.write(numString, sizeDouble*8);
	textStream<<" ";
}

void MobiusTransformation::LoadTransformation(std::ifstream & textStream)
{
	std::string tempStr;
	
	textStream>>tempStr;	
	assert(VKString(tempStr.c_str())=="MobiusTransformationEncoded:");
		//hmm... that's not what C++ specs say... cursor should be AFTER the space.
		//could be a problem on some other compilers (using macOS gcc)
	textStream.seekg(1, std::ios_base::cur);	
	int sizeDouble = sizeof(double);
	char * numString = new char[sizeDouble * 8];
	
	textStream.read(numString, sizeDouble*8);
	a.r = GetDouble(&(numString[sizeDouble*0]));
	a.i = GetDouble(&(numString[sizeDouble*1]));
	b.r = GetDouble(&(numString[sizeDouble*2]));
	b.i = GetDouble(&(numString[sizeDouble*3]));
	c.r = GetDouble(&(numString[sizeDouble*4]));
	c.i = GetDouble(&(numString[sizeDouble*5]));
	d.r = GetDouble(&(numString[sizeDouble*6]));
	d.i = GetDouble(&(numString[sizeDouble*7]));
}

void MobiusTransformation::SaveDouble(double val, char * saveTo)
{
	memcpy(saveTo, &val, sizeof(double));
}

double MobiusTransformation::GetDouble(const char * loadFrom)
{
	double val;
	memcpy(&val, loadFrom, sizeof(double));
	return val;
}

void MobiusTransformation::FindTransformation(std::vector<LinAlgComplex> &z, 
											  std::vector<LinAlgComplex> & w)
{
	if (z.size()==3 && w.size()==3)
		ExactFit(z[0], z[1], z[2], w[0], w[1], w[2]);
	else if (z.size()==w.size())
		LeastSquaresFit(z, w);
	else
		assert(false);
}

void MobiusTransformation::ExactFit(LinAlgComplex z1, LinAlgComplex z2, LinAlgComplex z3, 
				LinAlgComplex y1, LinAlgComplex y2, LinAlgComplex y3)
{
	// Compute 2x2 matrix for canonical vectors
	LinAlgComplex az = z2 - z3;;
	LinAlgComplex bz = z1*z3 - z1*z2;
	LinAlgComplex cz = z2 - z1;
	LinAlgComplex dz = z1*z3-z3*z2;
	
	// Compute 2x2 matrix for vertex vectors
	LinAlgComplex ay = y2 - y3;;
	LinAlgComplex by = y1*y3 - y1*y2;
	LinAlgComplex cy = y2 - y1;
	LinAlgComplex dy = y1*y3 - y3*y2;
	
	// Invert 2x2 matrix for vertex vectors
	LinAlgComplex dety = ay*dy - by*cy;
	assert(dety!=0);
	LinAlgComplex aiy =  dy / dety;
	LinAlgComplex biy = -by / dety;
	LinAlgComplex ciy = -cy / dety;
	LinAlgComplex diy =  ay / dety;
	
	// Multiply 2x2 matrices
	a = aiy*az + biy*cz;
	b = aiy*bz + biy*dz;
	c = ciy*az + diy*cz;
	d = ciy*bz + diy*dz;
	
//	std::cout<<"\nT(z1)=";
//	Transform(z1).PrintComplex();
//	std::cout<<"\nT(z2)=";
//	Transform(z2).PrintComplex();
//	std::cout<<"\nT(z3)=";
//	Transform(z3).PrintComplex();
}

bool MobiusTransformation::operator==(const MobiusTransformation & m) const
{
	return a==m.a && b==m.b && c==m.c && d==m.d;
}

bool MobiusTransformation::operator!=(const MobiusTransformation & m) const
{
	return a!=m.a || b!=m.b || c!=m.c || d!=m.d;
}

MobiusTransformation MobiusTransformation::operator*(const MobiusTransformation & m) const
{
	return MobiusTransformation(a * m.a + b * m.c,		// a
								a * m.b + b * m.d,		// b
								c * m.a + d * m.c,		// c
								c * m.b + d * m.d);		// d
}

MobiusTransformation MobiusTransformation::Inverse() const
{
	return MobiusTransformation(d / (a*d-b*c),
								-b / (a*d-b*c),
								-c / (a*d-b*c),
								a / (a*d-b*c));
}

void MobiusTransformation::Normalize()
{
//	std::cout<<"Normalizing Mobius: "<<std::endl;
//	Print();
	LinAlgComplex norm = a*d - b*c;
//	std::cout<<"norm=["<<norm.r<<", "<<norm.i<<"]"<<std::endl;	
	norm = norm.Sqrt();
//	std::cout<<"sqrt(norm)=["<<norm.r<<", "<<norm.i<<"]"<<std::endl;	
	a = a / norm;
	b = b / norm;
	c = c / norm;
	d = d / norm;	
//	std::cout<<"RESULT="<<std::endl;
//	Print();
//	std::cout<<"-------"<<std::endl;
}

LinAlgMatrixComplex MobiusTransformation::Matrix() const
{
	LinAlgMatrixComplex M(2,2);
	M.SetValue(0,0,a);
	M.SetValue(0,1,b);
	M.SetValue(1,0,c);
	M.SetValue(1,1,d);
	return M;
}

double MobiusTransformation::LeastSquaresFit(std::vector<LinAlgComplex> &z, std::vector<LinAlgComplex> & w)
{
	//assert(false);	// shouldn't need to use it... - not verified with the new math library
	
	assert(z.size()==w.size());
	LinAlgMatrixComplex * A[4];
	LinAlgVectorComplex * x[4];
	LinAlgVectorComplex * b[4];
	double error[4];
	for (int i=0; i<4; i++)
	{
		A[i] = new LinAlgMatrixComplex((int)z.size(), 3);
		x[i] = new LinAlgVectorComplex(3);
		b[i] = new LinAlgVectorComplex((int)z.size());
		error[i]=0;
	}
	
	for (int i=0; i<(int)z.size(); i++)
	{
		// set a=1
		A[0]->SetValue(i, 0, LinAlgComplex(1, 0));
		A[0]->SetValue(i, 1, -z[i] * w[i]);
		A[0]->SetValue(i, 2, -w[i]);
		b[0]->SetValue(i, -z[i]);
 
		// set b=1
		A[1]->SetValue(i, 0, z[i]);
		A[1]->SetValue(i, 1, -z[i] * w[i]);
		A[1]->SetValue(i, 2, -w[i]);
		b[1]->SetValue(i, LinAlgComplex(-1, -1));
		
		// set c=1
		A[2]->SetValue(i, 0, z[i]);
		A[2]->SetValue(i, 1, LinAlgComplex(1, 0));
		A[2]->SetValue(i, 2, -w[i]);
		b[2]->SetValue(i, -z[i]*w[i]);
	
		// set d=1
		A[3]->SetValue(i, 0, z[i]);
		A[3]->SetValue(i, 1, LinAlgComplex(1, 0));
		A[3]->SetValue(i, 2, -z[i]*w[i]);
		b[3]->SetValue(i, w[i]);
	}
	
	int bestID=-1;	
	for (int i=0; i<4; i++)
	{
		//A[i]->LeastSquaresFit(*(b[i]), *(x[i]), &(error[i]));
		A[i]->Solve(*(b[i]), *(x[i]));
//		if (error[i]!=error[i]) //nan
//		{
//			std::cout<<"[WARNING] Least Squares Failed (i="<<i<<")"<<std::endl;
//			error[i] = FLT_MAX;
//		}
//		else 
//		{
			// Find proper error
			FillMobiusTransformationFromMatrix('a'+i, *(x[i]));
			error[i] = (Transform(z[i])-w[i]).Power();
			
			// update best error
			if (bestID==-1 || error[bestID]>error[i])
				bestID = i;
//		}
	}
	
	FillMobiusTransformationFromMatrix('a'+bestID, *(x[bestID]));
	
	for (int i=0; i<4; i++)
	{
		delete A[i];
		delete b[i];
		delete x[i];
	}	
	assert(bestID>=0);
	return error[bestID];
}

LinAlgComplex MobiusTransformation::Transform(LinAlgComplex z)
{
	return (a * z + b) / (c * z + d);
}

LinAlgComplex MobiusTransformation::TransformInv(LinAlgComplex z)
{
	return (d * z - b) / (-c * z + a);
}

////////////////////////////////// SUPPPORTING FUNCTIONS ///////////////////////////////
void MobiusTransformation::FillMobiusTransformationFromMatrix(char setTo1, LinAlgVectorComplex & x)
{
	if (setTo1=='a')	// a=1 is the best
	{
		a = 1;
		b = x.GetValue(0);
		c = x.GetValue(1);
		d = x.GetValue(2);
	}
	else if (setTo1=='b')	// b=1 is the best
	{
		a = x.GetValue(0);
		b = 1;
		c = x.GetValue(1);
		d = x.GetValue(2);
	}
	else if (setTo1=='c')	// c=1 is the best
	{
		a = x.GetValue(0);
		b = x.GetValue(1);
		c = 1;
		d = x.GetValue(2);	
	}
	else if (setTo1=='d')	// d=1 is the best
	{
		a = x.GetValue(0);
		b = x.GetValue(1);
		c = x.GetValue(2);
		d = 1;
	}
	else 
	{
		assert(false);
	}
}


void MobiusTransformation::Print() const
{
	std::cout<<"\ta=";
	a.PrintComplex();
	std::cout<<"\n\tb=";
	b.PrintComplex();
	std::cout<<"\n\tc=";
	c.PrintComplex();
	std::cout<<"\n\td=";
	d.PrintComplex();
}


//////////////////////// MobiusTransformationInterpolator ///////////////////
MobiusTransformation MobiusTransformationInterpolator::InterpExpShiftI(
												std::vector<MobiusTransformation> & xforms, 
									            std::vector<double> & weights)
{
	int bestID=0;
	for (int i=0; i<(int)weights.size(); i++)
		if (weights[bestID] < weights[i])
			bestID = i;
//	return xforms[bestID];
	
//	std::cout<<"best="<<bestID<<std::endl;
	
	// retMob = e ^ (∑log(m_i))
	LinAlgMatrixComplex expSumLog(2, 2);
	LinAlgMatrixComplex sumLog(2, 2);
	LinAlgMatrixComplex log_m_i(2, 2);
	LinAlgMatrixComplex m_i(2, 2);

	LinAlgMatrixComplex m_bestInv(2, 2);
	LinAlgMatrixComplex m_best(2, 2);	
	LinAlgMatrixComplex m_i_temp(2, 2);
	
	xforms[bestID].Normalize();	
	m_best.SetValue(0, 0, xforms[bestID].a);
	m_best.SetValue(0, 1, xforms[bestID].b);
	m_best.SetValue(1, 0, xforms[bestID].c);		
	m_best.SetValue(1, 1, xforms[bestID].d);
	
	MobiusTransformation invMob(xforms[bestID]);
//	invMob.Inverse();
//	invMob.Normalize();
	m_bestInv.SetValue(0, 0, invMob.a);
	m_bestInv.SetValue(0, 1, invMob.b);
	m_bestInv.SetValue(1, 0, invMob.c);		
	m_bestInv.SetValue(1, 1, invMob.d);
	m_bestInv.Inverse();
	
//	std::cout<<"===================\n";	
//	std::cout<<"bestM = ";		m_best.PrintMatlab();	std::cout<<std::endl;
//	std::cout<<"A=";
	for (int i=0; i<(int)xforms.size(); i++)
	{
		xforms[i].Normalize();
		m_i.SetValue(0, 0, xforms[i].a);
		m_i.SetValue(0, 1, xforms[i].b);
		m_i.SetValue(1, 0, xforms[i].c);		
		m_i.SetValue(1, 1, xforms[i].d);
//		std::cout<<weights[i]<<"*logm(bestInv*";	m_i.PrintMatlab();	std::cout<<") +";		
		
//		if (weights[i]!=0)
//		{
//			std::cout<<"m_"<<i<<" = ";		m_i.PrintMatlab();	std::cout<<std::endl;		
//			std::cout<<"w_"<<i<<" = "<<weights[i]<<std::endl;
//		}
		m_bestInv.Multiply(m_i, m_i_temp);
		
		m_i_temp.Log(log_m_i);
		log_m_i.MultiplyByReal(weights[i]);
		sumLog.AddMatrix(log_m_i);
	}
	
//	std::cout<<"\n\nsumLog = ";	sumLog.PrintMatlab();
	
	sumLog.Exp(expSumLog);
//	std::cout<<"\n\nEXP= = ";	expSumLog.PrintMatlab();	
	
	m_best.Multiply(expSumLog, m_i_temp);
//	std::cout<<"\n\nFINAL= = ";	m_i_temp.PrintMatlab();		
	
	return MobiusTransformation(m_i_temp);
}

MobiusTransformation MobiusTransformationInterpolator::InterpExp(std::vector<MobiusTransformation> & xforms, 
																 std::vector<double> & weights)
{
	// Approach: 	Linear Combination of Transformations. Marc Alexa. Proc. Siggraph 2002
	//					http://portal.acm.org/citation.cfm?id=566592
	// Some computation: http://en.wikipedia.org/wiki/Logarithm_of_a_matrix
	// if A is diagonalizable: A' = V^-1 * A * V	(V - eigenvectors of A as columns, A' - diagonal)
		// ln (A) = V * ln(A') * V^-1
		// e^A = V * e^A' * V^−1. 	
	
	// if A is not diagonalizable use approach from the paper
	
//	std::cout<<"Interpolating EXP Mobius"<<std::endl;
	
	// retMob = e ^ (∑log(m_i))
	LinAlgMatrixComplex expSumLog(2, 2);
	LinAlgMatrixComplex sumLog(2, 2);
	LinAlgMatrixComplex log_m_i(2, 2);
	LinAlgMatrixComplex m_i(2, 2);
	
//	std::cout<<"=================== \nA = ";
	for (int i=0; i<(int)xforms.size(); i++)
	{
		xforms[i].Normalize();
		m_i.SetValue(0, 0, xforms[i].a);
		m_i.SetValue(0, 1, xforms[i].b);
		m_i.SetValue(1, 0, xforms[i].c);		
		m_i.SetValue(1, 1, xforms[i].d);
		
//		std::cout<<weights[i]<<"*logm(";	m_i.PrintMatlab();	std::cout<<") +";		
//		std::cout<<"m_"<<i<<" = ";		m_i.PrintMatlab();	std::cout<<std::endl;
		m_i.Log(log_m_i);
		
		log_m_i.MultiplyByReal(weights[i]);
		sumLog.AddMatrix(log_m_i);
	}
	
//	std::cout<<"\n\nsumLog = ";	sumLog.PrintMatlab();
	sumLog.Exp(expSumLog);
//	std::cout<<"\n\nEXP= = ";	expSumLog.PrintMatlab();
//	exit(0);
	return MobiusTransformation(expSumLog);
}

MobiusTransformation MobiusTransformationInterpolator::InterpNaive(std::vector<MobiusTransformation> & xforms, 
																   std::vector<double> & weights)
{
	// naive interpolation - replace
	MobiusTransformation retMob;
	
	retMob.a = 0.;
	retMob.b = 0.;
	retMob.c = 0.;
	retMob.d = 0.;	
	assert(xforms.size()==weights.size());
	
	for (int i=0; i<(int)xforms.size(); i++)
	{
		xforms[i].Normalize();
		
		retMob.a = retMob.a + xforms[i].a * weights[i];
		retMob.b = retMob.b + xforms[i].b * weights[i];
		retMob.c = retMob.c + xforms[i].c * weights[i];
		retMob.d = retMob.d + xforms[i].d * weights[i];
	}
	return retMob;	
}

//MobiusTransformationInterpolator::MobiusTransformationInterpolator(const MobiusTransformation & m1, 
//																   const MobiusTransformation & m2)
//: m_P(2,2), m_Pinv(2,2)
//{
//	assert(false);	// todo
//	LinAlgMatrixComplex Finv = m1.Inverse().Matrix();
//	LinAlgMatrixComplex G = m2.Matrix();
//
//	LinAlgMatrixComplex H(2, 2);
//	G.Multiply(Finv, H);
//
//	LinAlgComplex a = H.GetValue(0, 0);
//	LinAlgComplex b = H.GetValue(0, 1);
//	LinAlgComplex c = H.GetValue(1, 0);
//	LinAlgComplex d = H.GetValue(1, 1);	
//	
//	m_Eig1 = .5*a+.5*d+.5*(a*a-2.*a*d+d*d+4.*b*c).Sqrt();
//	//m_Eig2 = 1./2.*a+1./2.*d+1./2.*(a*a-2*a*d+d*d+4*b*c).sqrt();
//	
//	//assert(eig1>0 && eig2>0);
//	
//	m_P.SetValue(0, 0, 1);
//	m_P.SetValue(0, 1, 1);
//	m_P.SetValue(1, 0, -(.5*a-.5*d-.5*(a*a-2*a*d+d*a+4*b*c).Sqrt())/b);
//	m_P.SetValue(1, 1, -(.5*a-.5*d+.5*(a*a-2*a*d+d*d+4*b*c).Sqrt())/b);
//	
//	m_Pinv.SetValue(0, 0, .5*(a-d+(a*a-2*a*d+d*d+4*b*c).Sqrt())/(a*a-2*a*d+d*d+4*b*c).Sqrt());
//	m_Pinv.SetValue(0, 1, 1./(a*a-2*a*d+d*d+4*b*c).Sqrt()*b);
//	m_Pinv.SetValue(1, 0, -.5*(a-d-(a*a-2*a*d+d*d+4*b*c).Sqrt())/(a*a-2*a*d+d*d+4*b*c).Sqrt());
//	m_Pinv.SetValue(1, 1, -1./(a*a-2*a*d+d*d+4*b*c).Sqrt()*b);
//	
//}
//
//MobiusTransformation MobiusTransformationInterpolator::GetTransformation(double t)
//{
//	assert(false);	// todo
//	LinAlgMatrixComplex D(2, 2);
//	D.SetValue(1, 0, 0);
//	D.SetValue(0, 1, 0);
//	D.SetValue(0, 0, m_Eig1.RaiseToPow(t));
//	D.SetValue(1, 1, m_Eig1.RaiseToPow(-t));
//	LinAlgMatrixComplex temp(2, 2);
//	D.Multiply(m_Pinv, temp);	//D * P^-1
//	m_P.Multiply(temp, D);
//	
//	return MobiusTransformation(D);
//}










