#ifndef __PLANAR_TRANSFORM_QUASI_CONFORMAL_H
#define __PLANAR_TRANSFORM_QUASI_CONFORMAL_H

#include "PlanarTransform.h"
#include "MobiusTransformation.h"
#include "LinAlgMatrixComplex.h"
#include "LinAlgMatrixReal.h"

class PlanarTransformQuasiConformal : public PlanarTransform
{
	public:
		enum FPIGeneralization
		{
			FPI_NONE,
			FPI_MLS
		};
	
		PlanarTransformQuasiConformal(FPIGeneralization nPntCase = FPI_MLS);
	
		virtual LinAlgComplex Transform(LinAlgComplex z);
		virtual LinAlgComplex TransformInv(LinAlgComplex z);
		
		virtual void FindTransformation(std::vector<LinAlgComplex> &z, 
										std::vector<LinAlgComplex> & w);
	
		virtual void SaveTransformation(std::ofstream & textStream);
		virtual void LoadTransformation(std::ifstream & textStream);

		virtual ~PlanarTransformQuasiConformal(){}
	
	protected:
		LinAlgComplex LinearTransform(LinAlgComplex z);
	
		void FillAffineOn4Pnt(LinAlgVectorComplex & tZ, LinAlgVectorComplex & tW);
		void MapQuadrupletToParallelogram(std::vector<LinAlgComplex> & z, MobiusTransformation & m);
	
		FPIGeneralization m_NPointCase;
	
		int m_iNumConstraints;
		
		MobiusTransformation mz;
		MobiusTransformation mw;
		LinAlgMatrixReal L;
};
#endif

