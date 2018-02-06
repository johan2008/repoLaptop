#include "BasicNumerical/BasicNumerical.h"
LinAlgVectorComplex::LinAlgVectorComplex(int m)
{
	assert(m>0);
	m_aComplexData = new LinAlgComplex[m];
	//m_pData = gsl_vector_complex_alloc(m);
	this->m = m;
}

void LinAlgVectorComplex::PrintMatlab()
{
	std::cout<<"[ ";
	for (int i=0; i<m; i++)
	{
		//gsl_complex val = (*this)(i);
		std::string endChar=", ";
		if (i==m-1)
			endChar="";
		if (m_aComplexData[i].i>0)
			std::cout<<m_aComplexData[i].r<<"+"<<m_aComplexData[i].i<<"i"<<endChar.c_str();
		else if (m_aComplexData[i].i<0)
			std::cout<<m_aComplexData[i].r<<m_aComplexData[i].i<<"i"<<endChar.c_str();
		else
			std::cout<<m_aComplexData[i].r<<endChar.c_str();
	}
	std::cout<<"]"<<std::flush;
}

LinAlgVectorComplex::~LinAlgVectorComplex()
{
	delete [] m_aComplexData;
	//	gsl_vector_complex_free (m_pData);
}

LinAlgComplex & LinAlgVectorComplex::operator()(int i)
{
	assert(i<m && i>=0);
	return m_aComplexData[i];
}

LinAlgComplex LinAlgVectorComplex::GetValue(int i) const
{
	return (*(LinAlgVectorComplex*)this)(i);
	//return m_aComplexData[i];
	//return LinAlgComplex((*this)(i));
}

void LinAlgVectorComplex::SetValue(int i, const LinAlgComplex & value)
{
	(*this)(i) = value;
//	m_aComplexData[i] = value;
//	(*this)(i).dat[0] = value.r;
//	(*this)(i).dat[1] = value.i;	
}

//void LinAlgVectorComplex::SetValue(int i, double real, double imag)
//{
//	(*this)(i).dat[0] = real;
//	(*this)(i).dat[1] = imag;	
//}

//gsl_complex & LinAlgVectorComplex::operator()(int i)
//{
//	return *(gsl_vector_complex_ptr (m_pData, i));
//}

int LinAlgVectorComplex::Rows() const
{
	return m;
	//return m_pData->size;
}

