#ifndef __PLANAR_TRANSFORM_H
#define __PLANAR_TRANSFORM_H

#include "LinAlgComplex.h"
#include "gaps.h"

struct FlatSearchNode;

class PlanarTransform
{
	public:
		PlanarTransform(R3Mesh * flatMesh = NULL, 
						R2Kdtree<FlatSearchNode*> * flatMeshSearch=NULL);	// can only query points on the mesh vertices
	
		/**
		 * Apply transformation and its inverse
		 */
		virtual LinAlgComplex Transform(LinAlgComplex z) = 0;
		virtual LinAlgComplex TransformInv(LinAlgComplex z) = 0;
		
		/**
		 * Find transformation that takes points z to points w
		 */
		virtual void FindTransformation(std::vector<LinAlgComplex> &z, 
										std::vector<LinAlgComplex> & w) = 0;
	
		virtual void SaveTransformation(std::ofstream & textStream) = 0;
		virtual void LoadTransformation(std::ifstream & textStream) = 0;

		virtual ~PlanarTransform(){}
	
		static PlanarTransform * CreateTransform(std::ifstream & textStream, 
												 R3Mesh * flatMesh = NULL);
	
	protected:
		R3Mesh * m_pFlatMesh;
		R2Kdtree<FlatSearchNode*> * m_pFlatMeshSearch;
};
#endif

