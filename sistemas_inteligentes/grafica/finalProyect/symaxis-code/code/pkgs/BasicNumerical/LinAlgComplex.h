#ifndef __LINALG_COMPLEX_H
#define __LINALG_COMPLEX_H
#include <math.h>
class LinAlgComplex
{
	public:
		double r;
		double i;			
	
		LinAlgComplex(double r = 0, double i = 0);
	
#ifdef WITH_GSL
		LinAlgComplex(const gsl_complex & gslComplex);
		gsl_complex GetGSLComplex();
#endif
	

	// printing
		void PrintComplex() const ;
	
	
	// data operations
		double Power() const ;
		double Phase() const ;
	
		LinAlgComplex Sqrt() const;
		LinAlgComplex RaiseToPow(double val) const;
		LinAlgComplex Exp();

		LinAlgComplex Conjugate();
		
		bool operator == (const LinAlgComplex & c) const;
		bool operator == (double real) const;	

		bool operator != (const LinAlgComplex & c) const;
		bool operator != (double real) const;	

		LinAlgComplex & operator = (const LinAlgComplex & c);
		LinAlgComplex & operator = (double real);	
	
		LinAlgComplex operator * (const LinAlgComplex & c) const ;
		LinAlgComplex operator * (double real) const ;

		LinAlgComplex operator / (const LinAlgComplex & c) const ;
		LinAlgComplex operator / (double real) const ;

		LinAlgComplex operator + (const LinAlgComplex & c) const;
		LinAlgComplex operator + (double real) const ;

		LinAlgComplex operator - (const LinAlgComplex & c) const ;
		LinAlgComplex operator - (double real) const ;

		LinAlgComplex & operator *= (const LinAlgComplex & c);
		LinAlgComplex & operator *= (double real);

		LinAlgComplex & operator /= (const LinAlgComplex & c);
		LinAlgComplex & operator /= (double real);
	
		LinAlgComplex & operator += (const LinAlgComplex & c);
		LinAlgComplex & operator += (double real);

		LinAlgComplex & operator -= (const LinAlgComplex & c);
		LinAlgComplex & operator -= (double real);
	
		LinAlgComplex operator - () const;
};

LinAlgComplex operator * (double val1, LinAlgComplex val2);
LinAlgComplex operator / (double val1, LinAlgComplex val2);
LinAlgComplex operator + (double val1, LinAlgComplex val2);
LinAlgComplex operator - (double val1, LinAlgComplex val2);

#endif

