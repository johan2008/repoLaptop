#include "LinAlgComplex.h"


///////////////// non-member functions ////////////
LinAlgComplex operator * (double val1, LinAlgComplex val2)
{
	return val2 * val1;
}

LinAlgComplex operator / (double val1, LinAlgComplex val2)
{
	LinAlgComplex val1C(val1);
	return val1C / val2;
}

LinAlgComplex operator + (double val1, LinAlgComplex val2)
{
	return val2 + val1;
}

LinAlgComplex operator - (double val1, LinAlgComplex val2)
{
	return -val2 + val1;
}

///////////////// LinAlgComplex /////////////////
LinAlgComplex::LinAlgComplex(double r, double i)
{
	this->r = r;
	this->i = i;
}

LinAlgComplex LinAlgComplex::Exp()
{
	return exp(r) * LinAlgComplex(cos(i), sin(i));
}

#ifdef WITH_GSL
LinAlgComplex::LinAlgComplex(const gsl_complex & gslComplex)
{
	r = gslComplex.dat[0];
	i = gslComplex.dat[1];	
}

gsl_complex LinAlgComplex::GetGSLComplex()
{
	gsl_complex newComplex;
	newComplex.dat[0] = r;
	newComplex.dat[1] = i;
	return newComplex;
}
#endif

void LinAlgComplex::PrintComplex() const 
{
	if (i>=0)
		std::cout<<r<<"+"<<i<<"i"<<std::flush;
	else
		std::cout<<r<<"-"<<-i<<"i"<<std::flush;
}

double LinAlgComplex::Power() const 
{
	return sqrt(r*r + i*i);
}

double LinAlgComplex::Phase() const 
{
	if (r<0)
		return atan(i / r) + 3.14159265;
	else
		return atan(i / r);
}

LinAlgComplex LinAlgComplex::Conjugate()
{
	return LinAlgComplex(r, -i);
}

LinAlgComplex LinAlgComplex::Sqrt() const
{
	return RaiseToPow(.5);
}

LinAlgComplex LinAlgComplex::RaiseToPow(double val) const
{
	// [re, im] = e^(a + I*b)
	double a = log(Power());
	double t = atan2(i, r);
	return LinAlgComplex(exp(a*val)*cos(t*val), exp(a*val)*sin(t*val) );
}

bool LinAlgComplex::operator == (const LinAlgComplex & c) const
{
	return i==c.i && r==c.r;
}

bool LinAlgComplex::operator == (double real) const
{
	return i == 0 && r==real;
}

bool LinAlgComplex::operator != (const LinAlgComplex & c) const
{
	return i!=c.i || r!=c.r;
}

bool LinAlgComplex::operator != (double real) const
{
	return i!=0 || r!=real;
}

LinAlgComplex & LinAlgComplex::operator = (const LinAlgComplex & c)
{
	r = c.r;
	i = c.i;
	return *this;
}

LinAlgComplex & LinAlgComplex::operator = (double real)
{
	r = real;
	i = 0;
	return *this;
}

LinAlgComplex LinAlgComplex::operator * (const LinAlgComplex & c) const 
{	
	return LinAlgComplex(r * c.r - i * c.i, r*c.i + i * c.r);
}

LinAlgComplex LinAlgComplex::operator * (double real) const 
{
	return LinAlgComplex(r * real, i * real);
}

LinAlgComplex LinAlgComplex::operator / (const LinAlgComplex & c) const 
{
	double denom = c.r*c.r + c.i*c.i;
	if (denom == 0) 
		return LinAlgComplex(0, 0);
	
	return LinAlgComplex((r * c.r + i * c.i)/denom, (i * c.r - r * c.i)/denom) ;
}

LinAlgComplex LinAlgComplex::operator / (double real) const 
{
	return LinAlgComplex(r/real, i/real);
}

LinAlgComplex LinAlgComplex::operator + (const LinAlgComplex & c) const 
{
	return LinAlgComplex(r + c.r, i + c.i);
}

LinAlgComplex LinAlgComplex::operator + (double real) const 
{
	return LinAlgComplex(r + real, i);
}

LinAlgComplex LinAlgComplex::operator - (const LinAlgComplex & c) const 
{
	return LinAlgComplex(r - c.r, i - c.i);
}

LinAlgComplex LinAlgComplex::operator - (double real) const 
{
	return LinAlgComplex(r - real, i);	
}

LinAlgComplex &  LinAlgComplex::operator *= (const LinAlgComplex & c)
{
	r = r * c.r - i * c.i;
	i = r * c.i + i * c.r;
	return (*this);
}

LinAlgComplex & LinAlgComplex::operator *= (double real)
{
	r = r * real;
	i = i * real;
	return (*this);	
}

LinAlgComplex & LinAlgComplex::operator /= (const LinAlgComplex & c)
{
	double denom = c.r*c.r + c.i*c.i;
	if (denom == 0) 
	{
		r = 0;
		i = 0;
	}
	else
	{
		r = (r * c.r + i * c.i)/denom;
		i = (i * c.r - r * c.i)/denom;
	}
	return (*this);
}

LinAlgComplex & LinAlgComplex::operator /= (double real)
{
	r/=real;
	i/=real;
	return (*this);
}

LinAlgComplex & LinAlgComplex::operator += (const LinAlgComplex & c)
{
	r += c.r;
	i += c.i;
	return (*this);		
}

LinAlgComplex & LinAlgComplex::operator += (double real)
{
	r += real;
	return (*this);		
}

LinAlgComplex & LinAlgComplex::operator -= (const LinAlgComplex & c)
{
	assert(false);	// does not seem to return a desired result
	r -= c.r;
	i -= c.i;
	return (*this);		
}

LinAlgComplex & LinAlgComplex::operator -= (double real)
{
	r -= real;
	return (*this);	
}

LinAlgComplex LinAlgComplex::operator - () const
{
	return LinAlgComplex(-r, -i);
}



