#ifndef __PLANAR_TRANSFORM_LSCM_H
#define __PLANAR_TRANSFORM_LSCM_H

#include "PlanarTransform.h"
class FlatSearchNode;

class PlanarTransformLSCM : public PlanarTransform
{
	public:
		PlanarTransformLSCM(R3Mesh * flatMesh,
							R2Kdtree<FlatSearchNode*> * flatMeshSearch);
		virtual LinAlgComplex Transform(LinAlgComplex z);
		virtual LinAlgComplex TransformInv(LinAlgComplex z);
		
		virtual void FindTransformation(std::vector<LinAlgComplex> &z, 
										std::vector<LinAlgComplex> & w);
		
		virtual void SaveTransformation(std::ofstream & textStream);
		virtual void LoadTransformation(std::ifstream & textStream);
		
		virtual ~PlanarTransformLSCM(){}
	
		virtual void PrepareTransformation(int numConstraints);
		virtual void TransformAllPoints();
	
	protected:
		int m_iFreeVars;
		int m_iConstrainedVars;
	
		void PrepareInverseForNewPositions();
		std::vector<LinAlgComplex> m_NewVertexPositions;
		R2Kdtree<FlatSearchNode*> * m_pInverseFlatMeshSearch;
		FlatSearchNode * m_aAllNodes;
};
#endif

