#ifndef __PLANAR_TRANSFORM_MVC_H
#include "PlanarTransform.h"

class PlanarTransformMVC : public PlanarTransform
{
	public:
		PlanarTransformMVC();
		virtual LinAlgComplex Transform(LinAlgComplex z);
		virtual LinAlgComplex TransformInv(LinAlgComplex z);
		
		virtual void FindTransformation(std::vector<LinAlgComplex> &z, 
										std::vector<LinAlgComplex> & w);
		
		virtual void SaveTransformation(std::ofstream & textStream);
		virtual void LoadTransformation(std::ifstream & textStream);
		
		virtual ~PlanarTransformMVC(){}
		
	protected:
		std::vector<double> m_Z;
		std::vector<double> m_W;
	
		static void CalcMVC(std::vector<double> & ctrlPnts, double origX, double origY, 
							std::vector<double> & MVC);
		static void Interpolate(std::vector<double> & ctrlPnts, std::vector<double> & MVC,
								double & newX, double &newY);
	
};
#endif


