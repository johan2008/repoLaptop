#include "LinAlgVectorReal.h"
#include <iostream>
#include <fstream>
#include <sstream>


LinAlgVectorReal::LinAlgVectorReal(int m, double initVal)
{
	m_bFreeData = true;
	this->m = m;
	m_aRealData = new double[m];
	for (int i=0; i<m; i++)
		m_aRealData[i] = initVal;
	
}

LinAlgVectorReal::LinAlgVectorReal(std::vector<double> & inVector)
{
	m = inVector.size();
	m_bFreeData = true;
	m_aRealData = new double[m];
	for (int i=0; i<(int)inVector.size(); i++)
		(*this)(i) = inVector[i];
}

std::vector<double> LinAlgVectorReal::ToVector()
{
	std::vector<double> retvec;
	for (int i=0; i<m; i++)
		retvec.push_back((*this)(i));
	return retvec;
}

LinAlgVectorReal & LinAlgVectorReal::operator=(LinAlgVectorReal & vec)
{
	for (int i=0; i<m; i++)
		m_aRealData[i] = vec.m_aRealData[i];
	
	return *this;
}

void LinAlgVectorReal::Add(LinAlgVectorReal & v2)
{
	assert(v2.m==m);
	for (int i=0; i<m; i++)
		(*this)(i) += v2(i);
}

void LinAlgVectorReal::SetAll(double val)
{
	for (int i=0; i<m; i++)
		(*this)(i) = val;	
}

void LinAlgVectorReal::PrintMatlab(double absMin) const
{
	std::cout<<"[ ";
	for (int i=0; i<m; i++)
	{
		if (vkAbs((*((LinAlgVectorReal *)this))(i)) < absMin)
			std::cout<<0<<" ";
		else
			std::cout<<(*((LinAlgVectorReal *)this))(i)<<" ";			

	}
	std::cout<<"]"<<std::flush;
}

bool LinAlgVectorReal::WriteVectorBin(const char * outputFile)
{
	FILE * file = fopen(outputFile, "wb");
	assert (file!=NULL);

	fwrite(m_aRealData, sizeof(double), m, file);
	fclose(file);
	return true;
}

bool LinAlgVectorReal::ReadVectorBin(const char * inputFile)
{
	FILE * file = fopen(inputFile, "rb");
	if (file==NULL)
		return false;
	fread(m_aRealData, sizeof(double), m, file);
	fclose(file);
	
	return true;
}

bool LinAlgVectorReal::WriteVectorASCII(const char * outputFile)
{
	std::ofstream textStream(outputFile);
	assert(textStream.is_open());
	for (int i=0; i<m; i++)
		textStream<<(*this)(i)<<" ";
	textStream.close();
	return true;
}

bool LinAlgVectorReal::ReadVectorASCII(const char * inputFile)
{
	std::ifstream textStream(inputFile);
	if(!textStream.is_open())
		return false;
	for (int i=0; i<m; i++)
		textStream>>(*this)(i);
	textStream.close();
	return true;
}


LinAlgVectorReal::~LinAlgVectorReal()
{
	if (m_bFreeData)
		delete [] m_aRealData;
}

double & LinAlgVectorReal::operator()(int i)
{
	assert(i<m && i>=0);
	return m_aRealData[i];
}

double LinAlgVectorReal::GetValue(int i) const
{
	return (*(LinAlgVectorReal*)this)(i);
}

void LinAlgVectorReal::SetValue(int i, double value)
{
	(*this)(i) = value;
}


double LinAlgVectorReal::Norm(double p) const
{
	double currSum=0;
	for (int i=0; i<m; i++)
		currSum += pow((*((LinAlgVectorReal *)this))(i), p);
	return currSum;
}

void LinAlgVectorReal::Multiply(double val)
{
	for (int i=0; i<m; i++)
		(*this)(i) *= val;
}

int LinAlgVectorReal::Rows() const
{
	return m;
}




///////// GSL 
//	m_pData = gsl_vector_alloc(m);
//	gsl_vector_set_all (m_pData, initVal);
//	gsl_vector_memcpy(m_pData, vec.m_pData);
	//	gsl_vector_free (m_pData);
//return *(gsl_vector_ptr (m_pData, i));
