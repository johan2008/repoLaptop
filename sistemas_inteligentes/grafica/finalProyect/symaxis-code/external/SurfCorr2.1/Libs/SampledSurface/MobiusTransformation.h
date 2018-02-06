#include "LinAlgComplex.h"
#include "LinAlgVectorComplex.h"
#include "LinAlgMatrixComplex.h"
#include "PlanarTransform.h"

#ifndef __MOBIUS_TRANSFORMATION_H
#define __MOBIUS_TRANSFORMATION_H

class PlanarTransformQuasiConformal;

class MobiusTransformation : public PlanarTransform
{
	public:
		friend class MobiusTransformationInterpolator;
		friend class PlanarTransformQuasiConformal;
	
		MobiusTransformation();
		MobiusTransformation(LinAlgComplex a, LinAlgComplex b, LinAlgComplex c, LinAlgComplex d);
		MobiusTransformation(const LinAlgMatrixComplex & m);
		virtual ~MobiusTransformation(){}
		
		void SaveTransformation(std::ofstream & textStream);
		void LoadTransformation(std::ifstream & textStream);

		virtual void FindTransformation(std::vector<LinAlgComplex> &z, 
										std::vector<LinAlgComplex> & w);
	
		void ExactFit(LinAlgComplex z1, LinAlgComplex z2, LinAlgComplex z3, 
					  LinAlgComplex y1, LinAlgComplex y2, LinAlgComplex y3);
		double LeastSquaresFit(std::vector<LinAlgComplex> &z, std::vector<LinAlgComplex> & w);
	
		void Print() const;

		virtual LinAlgComplex Transform(LinAlgComplex z);
		virtual LinAlgComplex TransformInv(LinAlgComplex z);	
		void Normalize();
		LinAlgMatrixComplex Matrix() const;		
		MobiusTransformation Inverse() const;

		bool operator==(const MobiusTransformation & m) const;
		bool operator!=(const MobiusTransformation & m) const;

		MobiusTransformation operator*(const MobiusTransformation & m) const; // stacking MTs

	protected:
		LinAlgComplex a, b, c, d;

		void FillMobiusTransformationFromMatrix(char setTo1, LinAlgVectorComplex & x);
	
	private:
		void SaveDouble(double val, char * saveTo);
		double GetDouble(const char * loadFrom);
};

class MobiusTransformationInterpolator
{
	public:
//		MobiusTransformationInterpolator(const MobiusTransformation & m1, 
//										 const MobiusTransformation & m2);
//		MobiusTransformation GetTransformation(double t);
		
		static MobiusTransformation InterpExp(std::vector<MobiusTransformation> & xforms, 
											std::vector<double> & weights);
		static MobiusTransformation InterpExpShiftI(std::vector<MobiusTransformation> & xforms, 
											  std::vector<double> & weights);	
		static MobiusTransformation InterpNaive(std::vector<MobiusTransformation> & xforms, 
												std::vector<double> & weights);
	protected:
//		// diagonalized
//		LinAlgComplex m_Eig1;
//		LinAlgMatrixComplex m_P;
//		LinAlgMatrixComplex m_Pinv;
		
};

#endif
