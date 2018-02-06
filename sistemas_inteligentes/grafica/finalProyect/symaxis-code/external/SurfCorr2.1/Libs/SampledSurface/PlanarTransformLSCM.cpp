#include "PlanarTransformLSCM.h"
#include "MapFlattening.h"

PlanarTransformLSCM::PlanarTransformLSCM(R3Mesh * flatMesh, 
										 R2Kdtree<FlatSearchNode*> * flatMeshSearch)
: PlanarTransform(flatMesh, flatMeshSearch)
{
	assert(m_pFlatMesh!=NULL);
	assert(m_pFlatMeshSearch!=NULL);
	m_pInverseFlatMeshSearch = NULL;
	m_aAllNodes = NULL;
	m_iFreeVars = m_pFlatMesh->NVertices();
	m_iConstrainedVars = 0;
	for (int i=0; i<m_pFlatMesh->NVertices(); i++)
	{
		R3Point p = m_pFlatMesh->VertexPosition(m_pFlatMesh->Vertex(i));
		m_NewVertexPositions.push_back(LinAlgComplex(p.X(), p.Y()));
	}
	PrepareInverseForNewPositions();
}

LinAlgComplex PlanarTransformLSCM::Transform(LinAlgComplex z)
{
	FlatSearchNode * flatNode = m_pFlatMeshSearch->FindClosest(R2Point(z.r, z.i));
	assert(flatNode!=NULL);
	int id = flatNode->m_iID * 2;
	assert(id < (int)m_NewVertexPositions.size() && id >= 0);
	return m_NewVertexPositions[id];
}

LinAlgComplex PlanarTransformLSCM::TransformInv(LinAlgComplex z)
{
	FlatSearchNode * flatNode = m_pInverseFlatMeshSearch->FindClosest(R2Point(z.r, z.i));
	assert(flatNode!=NULL);
	int id = flatNode->m_iID * 2;
	assert(id < (int)m_NewVertexPositions.size() && id >= 0);
	return m_NewVertexPositions[id];
}

void PlanarTransformLSCM::FindTransformation(std::vector<LinAlgComplex> &z, 
											 std::vector<LinAlgComplex> & w)
{
	assert(z.size()==w.size());
	// re-initialize matrix system (if additional constraints were added)
	if ((int)z.size()*2 != m_iConstrainedVars)
		PrepareTransformation(z.size()*2);
	
	// Transform all vertices & store to: m_pInverseFlatMeshSearch, m_NewVertexPositions
	TransformAllPoints();
	PrepareInverseForNewPositions();
}

void PlanarTransformLSCM::SaveTransformation(std::ofstream & textStream)
{
	std::cout<<"[WARNIGN] PlanarTransformLSCM does not save!"<<std::endl;
}

void PlanarTransformLSCM::LoadTransformation(std::ifstream & textStream)
{
	std::cout<<"[WARNIGN] PlanarTransformLSCM does not load!"<<std::endl;
}

void PlanarTransformLSCM::PrepareTransformation(int numConstraints)
{

}

void PlanarTransformLSCM::TransformAllPoints()
{

}


void PlanarTransformLSCM::PrepareInverseForNewPositions()
{
	if (m_aAllNodes==NULL)
		m_aAllNodes = new FlatSearchNode[m_pFlatMesh->NVertices()];
	
	if (m_pInverseFlatMeshSearch!=NULL)
		delete m_pInverseFlatMeshSearch;
	
	RNArray<FlatSearchNode*> vts;
	for (int i=0; i<(int)m_NewVertexPositions.size(); i++)
	{
		m_aAllNodes[i].m_iID = i;
		m_aAllNodes[i].m_Pnt = R2Point(m_NewVertexPositions[i].r, m_NewVertexPositions[i].i);
		vts.InsertTail(&(m_aAllNodes[i]));
	}
	
	m_pInverseFlatMeshSearch = new R2Kdtree<FlatSearchNode*>(vts, &global_FlatNodePosition, NULL);
}



