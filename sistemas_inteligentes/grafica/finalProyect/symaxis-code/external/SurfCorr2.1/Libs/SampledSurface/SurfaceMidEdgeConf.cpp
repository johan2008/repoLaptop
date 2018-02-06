#include "BasicNumerical.h"
#include "SurfaceMidEdgeConf.h"
#include "AnalysisStats.h"
#include "MeshProcessor.h"

#include "DistanceGeodesic.h"
#include "FeatureAGD.h"

#define SEARCH_FOR_VERTICES
//#define SEARCH_FOR_TRIANGLES

FlatSearchNode::FlatSearchNode()
: m_iID(-1), m_Pnt(0, 0)
{
}

FlatSearchNode::FlatSearchNode(int id, const R2Point & pnt)
{
	m_iID = id;
	m_Pnt = pnt;
}
R2Point global_FlatNodePosition(FlatSearchNode * currNode, void * data)
{
	return currNode->m_Pnt;
}


int FaceNormalIsUp_global(const R3Point&, const R3Vector&, R3Mesh * mesh, R3MeshFace * found, void * reflected)
{
//	std::cout<<"CHeck if up"<<std::endl;
	bool refl = *((bool*)reflected);
	if ((mesh->FaceNormal(found).Z()>0 && !refl)
		|| (mesh->FaceNormal(found).Z()<0 && refl))
		return 1;
	else
		return 0;
}

//R2Point VertexPosition_globalFn(R3MeshVertex * vertex, void * data)
//{
//	R3Mesh * flatMesh = (R3Mesh*)data;
//	R3Point pnt = flatMesh->VertexPosition(vertex);
//	return R2Point(pnt.X(), pnt.Y());
//}

SurfaceMidEdgeConf::SurfaceMidEdgeConf(R3Mesh * mesh, bool complexFlip)
: SampledSurface(mesh)
{
	m_iMaxNumConfCoordsToCache = 2048;
	m_bUpdateTransformedMesh = false;	
	m_bReflection = complexFlip;
}

SurfaceMidEdgeConf::SurfaceMidEdgeConf(R3Mesh * mesh, int faceToRip, int canID1, int canID2, int canID3,
									   R3Mesh ** onlyFlattened, bool invert)
: SampledSurface(mesh)
{
	m_iMaxNumConfCoordsToCache = 2048;	
	m_bUpdateTransformedMesh = false;	
	m_bReflection = invert;
	*onlyFlattened = CreateFlattenedMidEdgeMesh(mesh, faceToRip);
	TransformToCanonicalFrame(mesh, *onlyFlattened, canID1, canID2, canID3);
	Initialize(mesh, *onlyFlattened, invert, true);
}

SurfaceMidEdgeConf::SurfaceMidEdgeConf(R3Mesh * mesh, R3Mesh * flatMesh, bool invert)
: SampledSurface(mesh)
{
	m_iMaxNumConfCoordsToCache = 2048;	
	m_bUpdateTransformedMesh = false;		
	m_bReflection = invert;	
	Initialize(mesh, flatMesh, invert, true);
}

R3Mesh * SurfaceMidEdgeConf::Flatten(int numSmooth, SurfaceFeature * agdFeature)
{
	int faceToRip, v1, v2, v3;
	R3Mesh * onlyFlattened;
	GetBestFlatteningParams(m_pMesh, faceToRip, v1, v2,  v3, agdFeature);
	
	if (numSmooth>0)
	{
		MeshProcessor processor(m_pMesh);
		processor.CopyAndForget();
		processor.Smooth(numSmooth);
		onlyFlattened = CreateFlattenedMidEdgeMesh(processor.m_pMesh, faceToRip);
	}
	else
		onlyFlattened = CreateFlattenedMidEdgeMesh(m_pMesh, faceToRip);
	TransformToCanonicalFrame(m_pMesh, onlyFlattened, v1, v2, v3);
	Initialize(m_pMesh, onlyFlattened, m_bReflection, true);
	return onlyFlattened;
}

MobiusTransformation SurfaceMidEdgeConf::GetCurrentTransform()
{
	return m_CurrentTransform;
}

bool SurfaceMidEdgeConf::Reflected()
{
	return m_bReflection;
}

// anotherlook
LinAlgComplex SurfaceMidEdgeConf::GetConfCoord(const SurfaceSample & s)
{
	if (m_bUpdateTransformedMesh)
		return GetConfCoord(s, m_pMesh, m_pFlatMeshTransformed);
	else
	{
		LinAlgComplex origFlat;
		bool atVertex;
		int nearestVertex = s.NearestVertex(&atVertex);
		std::map<int, LinAlgComplex>::iterator iter = m_iVertexToConfCoord.find(nearestVertex);
		if (atVertex && iter!=m_iVertexToConfCoord.end())
		{
			origFlat = iter->second;
		}
		else
		{
			origFlat = GetConfCoord(s, m_pMesh, m_pFlatMesh);
			if (atVertex && (int)m_iVertexToConfCoord.size() < m_iMaxNumConfCoordsToCache)
				m_iVertexToConfCoord[nearestVertex] = origFlat;
		}

		return m_CurrentTransform.Transform(origFlat);
	}
}

LinAlgComplex SurfaceMidEdgeConf::GetConfCoordOriginal(const SurfaceSample & s)
{
	return GetConfCoord(s, m_pMesh, m_pFlatMesh);
}

// anotherlook
SurfaceSample SurfaceMidEdgeConf::GetSurfaceSample(const LinAlgComplex & pos)
{	
#ifdef SEARCH_FOR_VERTICES
	assert(!m_bUpdateTransformedMesh);
	LinAlgComplex c = m_CurrentTransform.TransformInv(pos);
	return GetSurfaceSampleKdTree(m_CurrentTransform.TransformInv(pos), m_pFlatMesh, 
								  m_pFlatMeshKdVertexSearch, m_pMesh);	
#endif
#ifdef SEARCH_FOR_TRIANGLES
	if (m_bUpdateTransformedMesh)
		return GetSurfaceSampleSearchTree(pos, m_pFlatMeshTransformed, 
										  m_pFlatMeshTransSearch, m_pMesh, m_bReflection);
	else
	{
		LinAlgComplex c = m_CurrentTransform.TransformInv(pos);
		return GetSurfaceSampleSearchTree(m_CurrentTransform.TransformInv(pos), m_pFlatMesh, 
										  m_pFlatMeshTransSearch, m_pMesh, m_bReflection);
	}
#endif
}


// Lookme
SurfaceSample SurfaceMidEdgeConf::GetSurfaceSampleSearchTree(const LinAlgComplex & pos, 
															 R3Mesh * flatMesh, 
															 R3MeshSearchTree * flatMeshSearch, 
															 R3Mesh * mesh, bool reflected)
{
	assert(flatMeshSearch!=NULL);
	
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("ConfMapSearchFlatMesh");	
	R3MeshIntersection closest;
	flatMeshSearch->FindClosest(R3Point(pos.r, pos.i, 0), R3Vector(0,0,1), 
								closest, 0, FLT_MAX, &FaceNormalIsUp_global, &reflected);
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("ConfMapSearchFlatMesh");	
	
	if (closest.type!=R3_MESH_NULL_TYPE)
	{
		return MapMidEdgeToMesh(mesh, flatMesh, 
								SurfaceSample(closest, R3Point(pos.r, pos.i, 0), flatMesh));
	}
	std::cout<<"[ERROR] SurfaceMidEdgeConf: Search tree did not find triangle!"<<std::endl;
	assert(closest.type!=R3_MESH_NULL_TYPE);
	return SurfaceSample();	
}

SurfaceSample SurfaceMidEdgeConf::GetFlatSurfaceSampleKdTree(const LinAlgComplex & pos, 
															 R3Mesh * flatMesh, 
															 R2Kdtree<FlatSearchNode*> * flatMeshSearch)
{
	//	std::cout<<"Search KdTree - 1"<<std::endl;
	assert(flatMeshSearch!=NULL);
	R2Point query(pos.r, pos.i);
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("ConfMapSearchFlatMesh");		
	FlatSearchNode * flatNode = flatMeshSearch->FindClosest(query);
	if (flatNode==NULL)	// vk: XXX todo this is not supposed to be happening.
		return SurfaceSample(0, 1, 0, 0, flatMesh);
	
	assert(flatNode!=NULL);
	assert(flatNode->m_iID < flatMesh->NVertices() && flatNode->m_iID>=0);
	R3MeshVertex * flatVertex = flatMesh->Vertex(flatNode->m_iID);
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("ConfMapSearchFlatMesh");		
	if (flatVertex==NULL)
		std::cout<<"[ERROR] Invalid Query = ["<<pos.r<<", "<<pos.i<<"]"<<std::endl;
	assert(flatVertex!=NULL);
	// for every face adjacent to flatVertex
	int NFaces = flatMesh->VertexValence(flatVertex);
	assert(NFaces!=0);
	R3Point bary;
	bool foundBarycentric = false;
	R3MeshFace * face=NULL;
	for (int i=0; i<NFaces && !foundBarycentric; i++)
	{
		face = flatMesh->FaceOnVertex(flatVertex, 
									  flatMesh->EdgeOnVertex(flatVertex, i));
		if (face!=NULL)
		{
			bary = flatMesh->FaceBarycentric(face, R3Point(query.X(), query.Y(), 0));
			
			if (bary.X()>=0 && bary.Y()>=0 && bary.Z()>=0 
				&& bary.X()<=1 && bary.Y()<=1 && bary.Z()<=1
				&& (bary.X()!=0 || bary.Y()!=0 || bary.Z()!=0))
				foundBarycentric = true;
		}
	}
	if (face==NULL)
	{
		face = flatMesh->FaceOnVertex(flatVertex);
	}
	assert(face!=NULL);	
	if (!foundBarycentric)
	{		
		bary.X(0.);		bary.Y(0.);		bary.Z(0.);
		if (flatMesh->VertexOnFace(face, 0)==flatVertex)
			bary.X(1.);
		else if (flatMesh->VertexOnFace(face, 1)==flatVertex)
			bary.Y(1.);
		else if (flatMesh->VertexOnFace(face, 2)==flatVertex)
			bary.Z(1.);
	}
	else
	{
		double norm = bary.X()+bary.Y()+bary.Z();
		bary.X(bary.X()/norm);
		bary.Y(bary.Y()/norm);
		bary.Z(bary.Z()/norm);		
		
	}
	
	return SurfaceSample(flatMesh->FaceID(face), bary.X(), bary.Y(), bary.Z(), flatMesh);
}


SurfaceSample SurfaceMidEdgeConf::GetSurfaceSampleKdTree(const LinAlgComplex & pos, 
														 R3Mesh * flatMesh, 
														 R2Kdtree<FlatSearchNode*> * flatMeshSearch,
														 R3Mesh * mesh)
{	
	SurfaceSample flatSamp = GetFlatSurfaceSampleKdTree(pos, flatMesh, flatMeshSearch);
	SurfaceSample origSamp = MapMidEdgeToMesh(mesh, flatMesh, flatSamp);
	return origSamp;
}

void SurfaceMidEdgeConf::AddVerticesToFlatMesh(R3Mesh* mesh, R3Mesh * flatMesh)
{
	int probl2 = 0;
	for (int i=0; i<mesh->NVertices(); i++)
	{
		R3Point vertexPos(0,0,0);
		double norm = 0;
		
		R3MeshVertex * vertFlat = flatMesh->CreateVertex(vertexPos);
		for (int j=0; j < mesh->VertexValence(mesh->Vertex(i)); j++)
		{
			R3MeshEdge * e = mesh->EdgeOnVertex(mesh->Vertex(i), j);
			
			assert(mesh->EdgeID(e) < flatMesh->NVertices());
			R3Point diff = flatMesh->VertexPosition(flatMesh->Vertex(mesh->EdgeID(e)));
			vertexPos += diff;
			norm += 1.;
			
			R3MeshEdge * e2 = mesh->EdgeOnVertex(mesh->Vertex(i), e);
			
			if (e2==NULL)
				continue;
			if(mesh->FaceOnVertex(mesh->Vertex(i), e) != mesh->FaceBetweenEdges(e, e2))
				std::cout<<"[WARNING] 1. Problem occured when adding vertices to the flat mesh. "<<std::endl;
			
			R3MeshFace * newFace = NULL;
			if (newFace==NULL)
			{
				newFace = flatMesh->CreateFace(vertFlat, 
												  flatMesh->Vertex(mesh->EdgeID(e2)), 
												  flatMesh->Vertex(mesh->EdgeID(e)));
			}
			
			if (newFace==NULL)
			{
				newFace = flatMesh->CreateFace(vertFlat,
												  flatMesh->Vertex(mesh->EdgeID(e)), 
												  flatMesh->Vertex(mesh->EdgeID(e2)));
			}
			
			if(newFace==NULL)
				probl2++;
		}

		vertexPos.X(vertexPos.X()/norm);
		vertexPos.Y(vertexPos.Y()/norm);
		vertexPos.Z(vertexPos.Z()/norm);
		flatMesh->SetVertexPosition(vertFlat, vertexPos);
	}
	if (probl2>0)
	{
		std::cout<<"[WARNING] 2. Could not add "<<probl2<<" / "<<mesh->NVertices();
		std::cout<<" vertices to the flat mesh"<<std::endl;
	}
}

void SurfaceMidEdgeConf::Initialize(R3Mesh * mesh, R3Mesh * flatMeshTMP, bool complexFlip, bool copy)
{
	m_bReflection = complexFlip;
	m_pFlatMeshTransSearch = NULL;
	m_pFlatMeshTransformed = NULL;
	if (copy)
		m_pFlatMesh = SampledSurface::CreateCopy(flatMeshTMP);
	else 
		m_pFlatMesh = flatMeshTMP;
	
	
	AddVerticesToFlatMesh(mesh, m_pFlatMesh);
	
	if (complexFlip)
	{	
		for (int i=0; i<m_pFlatMesh->NVertices(); i++)
		{
			R3Point vPos = m_pFlatMesh->VertexPosition(m_pFlatMesh->Vertex(i));
			vPos.Y(-vPos.Y());
			m_pFlatMesh->SetVertexPosition(m_pFlatMesh->Vertex(i), vPos);
		}
	}	
	
	// Test mesh
	int incorrectNorm = 0;
	for (int i=0; i<m_pFlatMesh->NFaces(); i++)
	{
		R3MeshFace * face = m_pFlatMesh->Face(i);
		R3Vector norm = m_pFlatMesh->FaceNormal(face);
		if (! (m_bReflection ? (norm.Z()<0) : (norm.Z()>0)))
			incorrectNorm++;
	}
	if (incorrectNorm>0)
	{
		std::cout<<"[WARNING] Flattening contains foldovers. ";
		std::cout<<"NumTriangles="<<incorrectNorm<<" / "<<m_pFlatMesh->NFaces()<<std::endl;
	}
	
	if (m_bUpdateTransformedMesh)
		m_pFlatMeshTransformed = SampledSurface::CreateCopy(m_pFlatMesh);
	UpdateFlatSearchTree();
}

void SurfaceMidEdgeConf::Transform(const SurfaceSample & s1, const SurfaceSample & s2, const SurfaceSample & s3,
				const LinAlgComplex & t1, const LinAlgComplex & t2, const LinAlgComplex & t3)
{
	LinAlgComplex z1 = GetConfCoord(s1, m_pMesh, m_pFlatMesh);
	LinAlgComplex z2 = GetConfCoord(s2, m_pMesh, m_pFlatMesh);
	LinAlgComplex z3 = GetConfCoord(s3, m_pMesh, m_pFlatMesh);	

	m_CurrentTransform.ExactFit(z1, z2, z3, t1, t2, t3);
	if (m_bUpdateTransformedMesh)	
	{
		Transform(m_CurrentTransform, m_pFlatMesh, m_pFlatMeshTransformed);
		UpdateFlatSearchTree();
	}
}

void SurfaceMidEdgeConf::Transform(const MobiusTransformation & m)
{
	m_CurrentTransform = m;
	if (m_bUpdateTransformedMesh)
	{
		assert(m_pFlatMesh->NVertices()==m_pFlatMeshTransformed->NVertices());		
		Transform(m_CurrentTransform, m_pFlatMesh, m_pFlatMeshTransformed);
		UpdateFlatSearchTree();
	}
}

void SurfaceMidEdgeConf::LoadOriginalFlattening()
{
	MobiusTransformation identity;
	Transform(identity);
}

//lookme
LinAlgComplex SurfaceMidEdgeConf::GetConfCoord(const SurfaceSample & s, R3Mesh * mesh, R3Mesh * flat)
{
	assert(!s.Invalid());
	
	SurfaceSample flatSamp = MapMeshToMidEdge(mesh, flat, s);
	assert(!flatSamp.Invalid());
	
	R3Point pos = flatSamp.GetPosition();
	
	assert(!isnan(pos.X()));
	assert(!isnan(pos.Y()));
	return LinAlgComplex(pos.X(), pos.Y());
}

void SurfaceMidEdgeConf::MidEdgeToOriginalMap(int flatVertexID, int & edgeID, 
											  int & vertexID, R3Mesh * mesh)
{
	if (flatVertexID < mesh->NEdges())	//edge - return midpoint
	{
		edgeID = flatVertexID;
		vertexID = -1;
	//	return mesh->EdgeMidpoint(mesh->Edge(flatVertexID));
	}
	else
	{
		edgeID = -1;
		vertexID = flatVertexID - mesh->NEdges();
		assert(mesh->NVertices() > vertexID);
	//	return mesh->VertexPosition(mesh->Vertex(vertexID));
	}
}

int SurfaceMidEdgeConf::OriginalToMidEdgeMap(int v, R3Mesh * mesh, R3Mesh * flat)
{
	assert(v < mesh->NVertices());
	return v + mesh->NEdges();
}

int SurfaceMidEdgeConf::OriginalToMidEdgeMap(int v1, int v2, R3Mesh * mesh, R3Mesh * flat)
{
	assert(v1 < mesh->NVertices() && v2 < mesh->NVertices());
	R3MeshEdge * edge = mesh->EdgeBetweenVertices(mesh->Vertex(v1), mesh->Vertex(v2));
	assert(edge!=NULL);
	int edgeID = mesh->EdgeID(edge);
	return edgeID;
}

/** Outputs p - permutation according to local p */ 
int SurfaceMidEdgeConf::GetTriangleID(int vID[3], int p[3], R3Mesh * mesh)
{
	R3MeshVertex * v[3];
	for (int i=0; i<3; i++)
	{
		v[i] = mesh->Vertex(vID[i]);
		p[i] = -1;
	}
	
	int numFaces = mesh->VertexValence(v[0]);	
	for (int i=0; i<numFaces; i++)
	{
		R3MeshEdge * edge =  mesh->EdgeOnVertex(v[0], i);
		assert(edge!=NULL);
		R3MeshFace * face = mesh->FaceOnVertex(v[0], edge);
		if (face!=NULL && mesh->IsVertexOnFace(v[1], face) && mesh->IsVertexOnFace(v[2], face))
		{
			// find permutation
			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++)
					if (mesh->VertexOnFace(face, i)==v[j])
						p[j] = i;
//			std::cout<<"Permutation: "<<p[0]<<", "<<p[1]<<", "<<p[2]<<std::endl;
			assert(p[0]>=0 && p[1]>=0 && p[2]>=0 && p[0]!=p[1] && p[0]!=p[2] && p[1]!=p[2]);
			return mesh->FaceID(face);			
		}
	}
	
	return -1;
}

SurfaceSample SurfaceMidEdgeConf::MapMidEdgeToMeshNearVertex(R3Mesh * mesh, R3Mesh * flat, 
															 const SurfaceSample & samp, 
															 int localVertexID, 
															 int localToVertexID1, 
															 int localToVertexID2)
{
	double bary[3];
	int newTriangleID=-1;
	int vertexIDs[3];
	int vertPerm[3];
	
	// check if it is in triangle 1: v1, e2, e3
	vertexIDs[0] = OriginalToMidEdgeMap(samp.VertID(localVertexID), mesh, flat);
	
	vertexIDs[1] = OriginalToMidEdgeMap(samp.VertID(localVertexID), 
										samp.VertID(localToVertexID1), mesh, flat);

	vertexIDs[2] = OriginalToMidEdgeMap(samp.VertID(localVertexID), 
										samp.VertID(localToVertexID2), mesh, flat);

	
	newTriangleID = GetTriangleID(vertexIDs, vertPerm, flat);

	if (newTriangleID!=-1)
	{
		bary[vertPerm[1]] = samp.B(localToVertexID1) / .5;
		bary[vertPerm[2]] = samp.B(localToVertexID2) / .5;	
		bary[vertPerm[0]] = samp.B(localVertexID) - .5*(bary[vertPerm[1]] + bary[vertPerm[2]]);
		
		if (bary[vertPerm[0]] >= 0 && bary[vertPerm[0]] <= 1.
			&& bary[vertPerm[1]] >= 0 && bary[vertPerm[1]] <= 1.
			&& bary[vertPerm[2]] >= 0 && bary[vertPerm[2]] <= 1.)	
			return SurfaceSample(newTriangleID, bary[0], bary[1], bary[2], flat);
	}
	return SurfaceSample();
}

SurfaceSample SurfaceMidEdgeConf::MapMeshToMidEdge(R3Mesh * mesh, R3Mesh * flat, 
												   const SurfaceSample & samp)
{
	assert(samp.CheckMeshIsSame(mesh));
	
	SurfaceSample retVal;
	// check if it is in triangle 1: v1, e2, e3 
	retVal = MapMidEdgeToMeshNearVertex(mesh, flat, samp, 0, 1, 2);
	if (!retVal.Invalid())
		return retVal;

	// check if it is in triangle 2: v2, e1, e3 
	retVal = MapMidEdgeToMeshNearVertex(mesh, flat, samp, 1, 0, 2);
	if (!retVal.Invalid())
		return retVal;
	
	// check if it is in triangle 3: v3, e1, e2
	retVal = MapMidEdgeToMeshNearVertex(mesh, flat, samp, 2, 0, 1);
	if (!retVal.Invalid())
		return retVal;
	
	// check if it is in triangle 4: e1, e2, e3
	double bary[3];
	int newTriangleID=-1;
	int vertexIDs[3];
	int vertexPermutation[3];
	
	vertexIDs[0] = OriginalToMidEdgeMap(samp.VertID(1), samp.VertID(2), mesh, flat);
	vertexIDs[1] = OriginalToMidEdgeMap(samp.VertID(0), samp.VertID(2), mesh, flat);
	vertexIDs[2] = OriginalToMidEdgeMap(samp.VertID(0), samp.VertID(1), mesh, flat);
	newTriangleID = GetTriangleID(vertexIDs, vertexPermutation, flat);	
	if (newTriangleID!=-1)
	{
		bary[vertexPermutation[0]] = -samp.B(0) + samp.B(1) + samp.B(2);
		if (bary[vertexPermutation[0]]>=0 && bary[vertexPermutation[0]]<=1)
		{
			bary[vertexPermutation[1]] =  samp.B(0) - samp.B(1) + samp.B(2);
			if (bary[vertexPermutation[1]]>=0 && bary[vertexPermutation[1]]<=1)
			{
				bary[vertexPermutation[2]] =  samp.B(0) + samp.B(1) - samp.B(2);
				if (bary[vertexPermutation[2]]>=0 && bary[vertexPermutation[2]]<=1)
					return SurfaceSample(newTriangleID, bary[0], bary[1], bary[2], flat);
			}
		}	
	}	
	assert(false);
	return SurfaceSample();
}

int SurfaceMidEdgeConf::MidEdgeAdjacentOriginalVertex(int flatVertexFromEdgeID, int notThisVertex,
													  R3Mesh * mesh, R3Mesh * )
{
	assert(flatVertexFromEdgeID < mesh->NEdges());
	assert(notThisVertex < mesh->NVertices());
	R3MeshEdge * edge = mesh->Edge(flatVertexFromEdgeID);
	assert(edge!=NULL);
	R3MeshVertex * vertex = mesh->VertexAcrossEdge(edge, mesh->Vertex(notThisVertex));
	assert(vertex!=NULL);
	return mesh->VertexID(vertex);
}

int SurfaceMidEdgeConf::MidEdgeOpposingOriginalVertex(int edgeID, int otherEdgeID, 
													  R3Mesh * mesh, R3Mesh * )
{
	R3MeshEdge * otherEdge = mesh->Edge(otherEdgeID);
	assert(otherEdge!=NULL);
	R3MeshEdge * edge = mesh->Edge(edgeID);
	assert(edge!=NULL);	
	R3MeshVertex * v1 = mesh->VertexOnEdge(edge, 0);
	R3MeshVertex * v2 = mesh->VertexOnEdge(edge, 1);
	assert(v1!=NULL && v2!=NULL);
	
	R3MeshFace * face = mesh->FaceBetweenEdges(edge, otherEdge);
	for (int i=0; i<3; i++)
	{
		R3MeshVertex * thirdVertex = mesh->VertexOnFace(face, i);
		if (thirdVertex!=v1 && thirdVertex!=v2)
			return mesh->VertexID(thirdVertex);
	}
	assert(false);
	return -1;
}

SurfaceSample SurfaceMidEdgeConf::MapMidEdgeToMesh(R3Mesh * mesh, R3Mesh * flat, 
												   const SurfaceSample & samp)
{
	assert(samp.CheckMeshIsSame(flat));
	int v[3];
	int p[3];
	int e[3];
	double b[3];
	for (int i=0; i<3; i++)
		MidEdgeToOriginalMap(samp.VertID(i), e[i], v[i], mesh);
	
	if (v[0]!=-1 || v[1]!=-1 || v[2]!=-1)	// triangle near vertex
	{
		int locP[3];
		if (v[0]==-1 && v[2]==-1)
		{
			locP[0] = 1;
			v[0] = v[1];
		}
		else if (v[0]==-1 && v[1]==-1)
		{
			locP[0] = 2;
			v[0] = v[2];
		}
		else 
		{
			assert(v[1]==-1 && v[2]==-1);
			locP[0] = 0;
		}
		
		if (e[1]==-1)
		{
			locP[1] = 0;
			locP[2] = 2;
			e[1] = e[0];
		}
		else if (e[2]==-1)
		{
			locP[1] = 1;
			locP[2] = 0;
			e[2] = e[0];
		}
		else
		{
			assert(e[0]==-1);
			locP[1] = 1;
			locP[2] = 2;
		}
		
		v[1] = MidEdgeAdjacentOriginalVertex(e[1], v[0], mesh, flat);
		v[2] = MidEdgeAdjacentOriginalVertex(e[2], v[0], mesh, flat);		
		
		int triangleID = GetTriangleID(v, p, mesh);
		assert(triangleID!=-1);
		b[p[1]] = samp.B(locP[1]) * .5;
		b[p[2]] = samp.B(locP[2]) * .5;		
		b[p[0]] = samp.B(locP[0]) + (samp.B(locP[1]) + samp.B(locP[2])) * .5;
		
		assert(b[0]>=0 && b[0]<=1 && b[1]>=0 && b[1]<=1 && b[2]>=0 && b[2]<=1);
		assert(vkAbs(b[0] + b[1] + b[2] - 1.)<.00001);
		return SurfaceSample(triangleID, b[0], b[1], b[2], mesh);
	}
	else		// middle triangle
	{
		for (int i=0; i<3; i++)
		{
			assert(v[i]==-1);
			assert(e[i]!=-1);
		}
		v[0] = MidEdgeOpposingOriginalVertex(e[0], e[1], mesh, flat);
		v[1] = MidEdgeOpposingOriginalVertex(e[1], e[0], mesh, flat);
		v[2] = MidEdgeOpposingOriginalVertex(e[2], e[1], mesh, flat);
		int triangleID = GetTriangleID(v, p, mesh);
		assert(triangleID!=-1);
		b[p[0]] = 0.5 * (samp.B(1) + samp.B(2));
		b[p[1]] = 0.5 * (samp.B(0) + samp.B(2));
		b[p[2]] = 0.5 * (samp.B(0) + samp.B(1));
		assert(b[0]>=0 && b[0]<=1 && b[1]>=0 && b[1]<=1 && b[2]>=0 && b[2]<=1);
		assert(vkAbs(b[0] + b[1] + b[2] - 1.)<.00001);
		return SurfaceSample(triangleID, b[0], b[1], b[2], mesh);
	}
	assert(false);
	return SurfaceSample();
}

//////////////////////// FLATTENING //////////////////////////////////
void SurfaceMidEdgeConf::GetBestFlatteningParams(R3Mesh * mesh, int & faceToRip, 
												 int & v1, int & v2, int & v3, SurfaceFeature * agd)
{	
	assert(agd!=NULL);
	int minVertex=-1;
	double minValue = -1;
	for (int i=0; i<mesh->NVertices(); i++)	// find minimal vertex (face to rip)
	{
		double currVal = agd->Value(SurfaceSample(i, mesh));
		if (minVertex==-1 || currVal < minValue)
		{
			minValue = currVal;
			minVertex = i;
		}
	}
	assert(minVertex!=-1);
	
	DistanceGeodesic::GeoVertexData * vertexData = DistanceGeodesic::InitializeFunkPrecomputation(mesh, false);	
	double * distToVertex = new double[mesh->NVertices()];
	double * minDistanceToSet = new double[mesh->NVertices()];

	v1 = 0;			v2 = 0;			v3 = 0;
	// find furthest from min(agd)
	DistanceGeodesic::PrecomputeFillFunkhouser(vertexData, mesh, minVertex, distToVertex);
	for (int i=0; i<mesh->NVertices(); i++)
	{
		minDistanceToSet[i] = distToVertex[i];
		if (minDistanceToSet[i]>minDistanceToSet[v1])
			v1 = i;
	}
	// find furthest from MINDIST(min(agd), v1)
	DistanceGeodesic::PrecomputeFillFunkhouser(vertexData, mesh, v1, distToVertex);
	for (int i=0; i<mesh->NVertices(); i++)
	{
		minDistanceToSet[i] = vkMin(distToVertex[i], minDistanceToSet[i]);
		if (minDistanceToSet[i]>minDistanceToSet[v2])
			v2 = i;
	}
	// find furthest from MINDIST(min(agd), v1, v2)
	DistanceGeodesic::PrecomputeFillFunkhouser(vertexData, mesh, v2, distToVertex);
	for (int i=0; i<mesh->NVertices(); i++)
	{
		minDistanceToSet[i] = vkMin(distToVertex[i], minDistanceToSet[i]);
		if (minDistanceToSet[i]>minDistanceToSet[v3])
			v3 = i;
	}
	
	// find face to rip
	faceToRip = minVertex;
	assert ((faceToRip<mesh->NVertices() && faceToRip>=0));
	assert ((mesh->VertexValence(mesh->Vertex(faceToRip))!=0));
	assert ((mesh->FaceOnVertex(mesh->Vertex(faceToRip))!=NULL));

	faceToRip = mesh->FaceID(mesh->FaceOnVertex(mesh->Vertex(faceToRip)));	
	delete [] vertexData;
	delete [] distToVertex;
	delete [] minDistanceToSet;
}

void SurfaceMidEdgeConf::TransformToCanonicalFrame(R3Mesh *original_mesh, R3Mesh *flattened_mesh, 
						int v1ID, int v2ID, int v3ID)
{
	// Pick three vertices (this is not right)
	R3MeshVertex *v1 = original_mesh->Vertex(v1ID);
	R3MeshVertex *v2 = original_mesh->Vertex(v2ID);
	R3MeshVertex *v3 = original_mesh->Vertex(v3ID);
	R3MeshEdge *e1 = original_mesh->EdgeOnVertex(v1);
	R3MeshEdge *e2 = original_mesh->EdgeOnVertex(v2);
	R3MeshEdge *e3 = original_mesh->EdgeOnVertex(v3);
	R3MeshVertex *ve1 = flattened_mesh->Vertex(original_mesh->EdgeID(e1));
	R3MeshVertex *ve2 = flattened_mesh->Vertex(original_mesh->EdgeID(e2));
	R3MeshVertex *ve3 = flattened_mesh->Vertex(original_mesh->EdgeID(e3));
	R3Vector z1 = flattened_mesh->VertexPosition(ve1).Vector();
	R3Vector z2 = flattened_mesh->VertexPosition(ve2).Vector();
	R3Vector z3 = flattened_mesh->VertexPosition(ve3).Vector();

	LinAlgComplex y1(1, 0);
	LinAlgComplex y2(cos(2*3.14159265/3), sin(2*3.14159265/3));	
	LinAlgComplex y3(cos(4*3.14159265/3), sin(4*3.14159265/3));		

	MobiusTransformation m;
	m.ExactFit(LinAlgComplex(z1.X(), z1.Y()), 
			   LinAlgComplex(z2.X(), z2.Y()), 
			   LinAlgComplex(z3.X(), z3.Y()), 
			   y1, y2, y3);
	
	Transform(m, flattened_mesh, flattened_mesh);
}

void SurfaceMidEdgeConf::Transform(MobiusTransformation m, R3Mesh * meshIn, R3Mesh * meshOut)
{
	assert(meshIn->NVertices()==meshOut->NVertices());
	for (int i=0; i<meshIn->NVertices(); i++)
	{
		R3Point pos = meshIn->VertexPosition(meshIn->Vertex(i));
		LinAlgComplex newPos = m.Transform(LinAlgComplex(pos.X(), pos.Y()));
		meshOut->SetVertexPosition(meshOut->Vertex(i), R3Point(newPos.r, newPos.i, 0));
	}
}

double SurfaceMidEdgeConf::Cotangent(const R3Point& position0, const R3Point& position1, const R3Point& position2)
{
	// Return cotan of angle formed by vertex at position1
	R3Vector vec1 = position0 - position1; vec1.Normalize();
	R3Vector vec2 = position2 - position1; vec2.Normalize();
	RNAngle angle = R3InteriorAngle(vec1, vec2);
	if (angle == 0) return FLT_MAX;
	double tan_angle = tan(angle);
	if (tan_angle == 0) return FLT_MAX;
	return 1.0 / tan_angle;
}

LinAlgVectorReal * SurfaceMidEdgeConf::ComputeXCoordinates(R3Mesh *mesh, R3MeshFace *remove_face)
{	
	// Get the indices of the fixed vertices
	int fixed_vertex_id0 = mesh->VertexID(mesh->VertexOnFace(remove_face, 0));
	int fixed_vertex_id1 = mesh->VertexID(mesh->VertexOnFace(remove_face, 1));
	
	// Allocate data
	int n = mesh->NVertices();
	LinAlgMatrixSparseReal A(n, n, LinAlgMatrixSparseReal::SUPPORT_CX_SPARSE);
	LinAlgVectorReal * x = new LinAlgVectorReal(n);
	LinAlgVectorReal b(n);
	
	// Initialize data
	for (int i = 0; i < n; i++) 
	{
		(*x)(i) = 0;
		b(i) = 0;
	}

	std::map<int, bool> nonDelaunay;
	// Fill matrix A
	for (int i1 = 0; i1 < n; i1++) 
	{
		R3MeshVertex *v1 = mesh->Vertex(i1);
		const R3Point& p1 = mesh->VertexPosition(v1);
		for (int j = 0; j < mesh->VertexValence(v1); j++) 
		{
			R3MeshEdge *e = mesh->EdgeOnVertex(v1, j);
			R3MeshVertex *v2 = mesh->VertexAcrossEdge(e, v1);
			const R3Point& p2 = mesh->VertexPosition(v2);
			int i2 = mesh->VertexID(v2);
			double cot = 0;
			for (int k = 0; k < 2; k++) 
			{
				R3MeshFace *f = mesh->FaceOnEdge(e, k);
				if (!f) continue;
				R3MeshVertex *v3 = mesh->VertexAcrossFace(f, e);
				const R3Point& p3 = mesh->VertexPosition(v3);
				cot += Cotangent(p1, p3, p2);
			}
			if (cot < 0)
				nonDelaunay[mesh->EdgeID(e)]=true;
			if (i1!=fixed_vertex_id0 && i1!=fixed_vertex_id1)
			{
				A.AddVal(i1, i1, cot);
				A.AddVal(i1, i2, -cot);
			}
		}
	}

	if (nonDelaunay.size()>0)
		std::cout<<"[WARNING] Flattening with non-delaunay edges: "<<nonDelaunay.size()<<std::endl;
	
	// Fix two vertices
	A.AddVal(fixed_vertex_id0, fixed_vertex_id0, 1);
	A.AddVal(fixed_vertex_id1, fixed_vertex_id1, 1);
	b(fixed_vertex_id0) = -1;
	b(fixed_vertex_id1) = 1;
	
	//std::cout<<"\tMid Edge Flattening: solving linear system."<<std::endl;
	A.Solve(b, *x);
	// Return x coordinates
	return x;
}



LinAlgVectorReal * SurfaceMidEdgeConf::ComputeYCoordinates(R3Mesh *mesh, R3MeshFace *remove_face, 
														   LinAlgVectorReal &x)
{	
	// Allocate y coordinates
	LinAlgVectorReal *y = new LinAlgVectorReal(mesh->NEdges());
	if (!y) { fprintf(stderr, "Unable to allocate data for Y coordinates\n"); return NULL; }
	for (int i = 0; i < mesh->NEdges(); i++) (*y)(i) = 0;
	
	// Compute y coordinates
	R3mesh_mark++;
	R3MeshEdge *fixed_edge = mesh->EdgeOnFace(remove_face, 0);
	mesh->SetEdgeMark(fixed_edge, R3mesh_mark);
	RNQueue<R3MeshEdge *> queue;
	queue.Push(fixed_edge);
	while (!queue.IsEmpty()) {
		R3MeshEdge *edge1 = queue.Pop();
		int edge_id1 = mesh->EdgeID(edge1);
		for (int i = 0; i < 2; i++) {
			R3MeshFace *face = mesh->FaceOnEdge(edge1, i);
			if (!face) continue;
			if (face == remove_face) 
				continue;

			for (int dir = RN_CCW; dir <= RN_CW; dir++) {
				R3MeshEdge *edge0 = mesh->EdgeOnFace(face, edge1, dir);
				if (mesh->EdgeMark(edge0) != R3mesh_mark) {
					mesh->SetEdgeMark(edge0, R3mesh_mark);
					int edge_id0 = mesh->EdgeID(edge0);
					R3MeshVertex *vertex0 = mesh->VertexOnFace(face, edge1, dir);
					R3MeshVertex *vertex1 = mesh->VertexAcrossEdge(edge1, vertex0);
					R3MeshVertex *vertex2 = mesh->VertexAcrossFace(face, edge1);
					const R3Point& position0 = mesh->VertexPosition(vertex0);
					const R3Point& position1 = mesh->VertexPosition(vertex1);
					const R3Point& position2 = mesh->VertexPosition(vertex2);
					int vertex_id0 = mesh->VertexID(vertex0);
					int vertex_id1 = mesh->VertexID(vertex1);
					int vertex_id2 = mesh->VertexID(vertex2);
					RNScalar term1 = (x(vertex_id2) - x(vertex_id0)) * Cotangent(position0, position1, position2);
					RNScalar term2 = (x(vertex_id1) - x(vertex_id0)) * Cotangent(position1, position2, position0);
					RNScalar sign = (dir == RN_CCW) ? -1 : 1;
					RNScalar delta_y = sign * 0.5 * (term1 + term2);
					(*y)(edge_id0) = (*y)(edge_id1) + delta_y;
					queue.Push(edge0);
				}
			}
		}
	}
		
	// Return y coordinates
	return y;
}

R3Mesh * SurfaceMidEdgeConf::CreateFlattenedMidEdgeMesh(R3Mesh *mesh, R3MeshFace *remove_face, 
						LinAlgVectorReal *x, LinAlgVectorReal *y)
{
	// Allocate mesh
	R3Mesh *midedge_mesh = new R3Mesh();
	if (!midedge_mesh) {
		fprintf(stderr, "Unable to allocate mesh\n");
		return NULL;
	}
	
	// Create vertices 
	for (int i = 0; i < mesh->NEdges(); i++) {
		R3MeshEdge *edge = mesh->Edge(i);
		R3MeshVertex *vertex0 = mesh->VertexOnEdge(edge, 0);
		R3MeshVertex *vertex1 = mesh->VertexOnEdge(edge, 1);
		int vertex_id0 = mesh->VertexID(vertex0);
		int vertex_id1 = mesh->VertexID(vertex1);
		RNScalar midedge_x = 0.5 * ((*x)(vertex_id0) + (*x)(vertex_id1));
		R3Point midedge_position(midedge_x, (*y)(i), 0);
		midedge_mesh->CreateVertex(midedge_position);
	}
	
	// Create faces
	for (int i = 0; i < mesh->NFaces(); i++) {
		R3MeshFace *face = mesh->Face(i);
		if (face == remove_face
			|| face==mesh->FaceOnFace(remove_face, 0)
			|| face==mesh->FaceOnFace(remove_face, 1)
			|| face==mesh->FaceOnFace(remove_face, 2)
			)
			continue;
		R3MeshEdge *edge0 = mesh->EdgeOnFace(face, 0);
		R3MeshEdge *edge1 = mesh->EdgeOnFace(face, 1);
		R3MeshEdge *edge2 = mesh->EdgeOnFace(face, 2);
		R3MeshVertex *midedge_vertex0 = midedge_mesh->Vertex(mesh->EdgeID(edge0));
		R3MeshVertex *midedge_vertex1 = midedge_mesh->Vertex(mesh->EdgeID(edge1));
		R3MeshVertex *midedge_vertex2 = midedge_mesh->Vertex(mesh->EdgeID(edge2));
		
		//R3MeshFace * midedge_face = 
		midedge_mesh->CreateFace(midedge_vertex0, midedge_vertex1, midedge_vertex2);
	}
	
	// Return midedge mesh
	return midedge_mesh;
}



R3Mesh * SurfaceMidEdgeConf::CreateFlattenedMidEdgeMesh(R3Mesh *mesh, int removeFaceID)
{	
	// Pick face to exclude
	R3MeshFace *remove_face = mesh->Face(removeFaceID);
	
	// Compute X coordinates
	LinAlgVectorReal * x = ComputeXCoordinates(mesh, remove_face);

	// Compute Y coordinates
	LinAlgVectorReal *y = ComputeYCoordinates(mesh, remove_face, *x);

	// Create mid-edge mesh
	R3Mesh *midedge_mesh = CreateFlattenedMidEdgeMesh(mesh, remove_face, x, y);
	if (!midedge_mesh) return NULL;
		
	// Delete coordinates
	delete x;
	delete y;

	// Return mid-edge mesh
	return midedge_mesh;
}

R2Kdtree<FlatSearchNode*>* SurfaceMidEdgeConf::CreateKdTree(R3Mesh * flattenedMesh, 
													FlatSearchNode ** vertexArray)
{
	*vertexArray = new FlatSearchNode[flattenedMesh->NVertices()];
	RNArray<FlatSearchNode*> vts;
	for (int i=0; i<flattenedMesh->NVertices(); i++)
	{
		R3Point pos = flattenedMesh->VertexPosition(flattenedMesh->Vertex(i));
		(*vertexArray)[i].m_iID = i;
		(*vertexArray)[i].m_Pnt = R2Point(pos.X(), pos.Y());
		vts.InsertTail(&((*vertexArray)[i]));
	}
	
	return new R2Kdtree<FlatSearchNode*>(vts, &global_FlatNodePosition, NULL);
}

void SurfaceMidEdgeConf::UpdateFlatSearchTree()
{
	if (m_bUpdateTransformedMesh)
	{
#ifdef SEARCH_FOR_VERTICES
		assert(false);
#endif
#ifdef SEARCH_FOR_TRIANGLES
		if (m_pFlatMeshTransSearch!=NULL)
			delete m_pFlatMeshTransSearch;
		m_pFlatMeshTransSearch = new R3MeshSearchTree(m_pFlatMeshTransformed);
#endif
	}
	else if (m_pFlatMeshTransSearch==NULL)
	{
#ifdef SEARCH_FOR_VERTICES
		m_pFlatMeshKdVertexSearch = CreateKdTree(m_pFlatMesh, &m_VertexArray);
#endif
#ifdef SEARCH_FOR_TRIANGLES
		m_pFlatMeshTransSearch = new R3MeshSearchTree(m_pFlatMesh);
#endif
	}
}


void SurfaceMidEdgeConf::SetParentWindow(AnalysisWindow * window)
{
	m_pWindow = window;
	m_iCameraID = m_pWindow->AddMesh(GetMesh());
	if (m_bUpdateTransformedMesh)
		assert(m_pWindow->AddMesh(m_pFlatMeshTransformed, m_iCameraID)==m_iCameraID);
	else
		assert(m_pWindow->AddMesh(m_pFlatMesh, m_iCameraID)==m_iCameraID);
}

VKString SurfaceMidEdgeConf::GetSurfaceType()
{
	return "SurfaceMidEdgeConf";
}

////////////// RENDERING /////////////
double SurfaceMidEdgeConf::GetStandardRadius(ParamParser * params, 
											 const VKString & renderingParams,
											 const VKString & surfaceName)
{
	assert(params!=NULL);
	if (params->GetStrValue(renderingParams, "ConfSurf", valid)=="RenderOriginal")
		return 0.005 * GetMesh()->BBox().DiagonalLength();
	else
	{
		if (m_pFlatMeshTransformed==NULL)
			return 0.001 * m_pFlatMesh->BBox().DiagonalLength();
		else
			return 0.001 * m_pFlatMeshTransformed->BBox().DiagonalLength();
	}
}


void SurfaceMidEdgeConf::Draw(ParamParser * params, 
							  const VKString & renderingParams,
							  const VKString & surfaceName)
{
	if (params->GetStrValue(renderingParams, "ConfSurf", valid)=="RenderOriginal")
	{
		m_pWindow->SetCurrentMeshInCollection(0);		
		DrawSurface(params, renderingParams, surfaceName);
	}
	else if (params->GetStrValue(renderingParams, "ConfSurf", valid)=="RenderConformal")
	{
		m_pWindow->SetCurrentMeshInCollection(1);
		DrawConformal(params, renderingParams, surfaceName);
	}
}

void SurfaceMidEdgeConf::LoadConformalCamera()
{
	m_pWindow->LoadCamera(m_iCameraID, 1);	
}

void SurfaceMidEdgeConf::DrawConformalDoNotLoadCamera(ParamParser * params, 
													  const VKString & renderingParams,
													  const VKString & surfaceName)
{
	R3Mesh * flatMesh = m_pFlatMeshTransformed;
	if (m_pWindow->IsFlagSet(renderingParams, "MeshFlags", "Faces", 
							 surfaceName, "MeshFlagsFlips")
		&& !m_bUpdateTransformedMesh)
	{
		static bool warnOnce = true;
		if (warnOnce)
		{
			std::cout<<"[WARNING] Trying to render conformal maps, but impossible due to speedup. Set m_bUpdateTransformedMesh=true"<<std::endl;
			warnOnce = false;
		}
		flatMesh = m_pFlatMesh;
	}
		
	if (m_pWindow->IsFlagSet(renderingParams, "MeshFlags", "Faces", 
							 surfaceName, "MeshFlagsFlips"))
	{
		assert(flatMesh!=NULL);
		glEnable(GL_LINE_STIPPLE); 
		glLineStipple(1, 0xF00F); 
		
		glDisable(GL_LIGHTING);
		glLineWidth(3);
		glBegin(GL_LINES);
		glColor3d(0, 0, 1);
		glVertex3d(0, 0, 0);
		glVertex3d(2, 0, 0);
		glEnd();
		glDisable(GL_LINE_STIPPLE); 
		
		m_pWindow->SetLights(GetMesh());
		static GLfloat material[4];	
		material[3] = .8;
		glEnable(GL_LIGHTING);
		glBegin(GL_TRIANGLES);
		for (int i=0; i<flatMesh->NFaces(); i++)
		{
			if (i < m_pMesh->NFaces())
			{	material[0] = 0;	material[1] = 0;	material[2] = 1;	}
			else
			{	material[0] = .5;	material[1] = .8;	material[2] = .8;	}
				
			
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material); 			
			R3MeshFace * face = flatMesh->Face(i);
			
			R3Point p1 = flatMesh->VertexPosition(flatMesh->VertexOnFace(face, 0));
			R3Point p2 = flatMesh->VertexPosition(flatMesh->VertexOnFace(face, 1));
			R3Point p3 = flatMesh->VertexPosition(flatMesh->VertexOnFace(face, 2));
			
			double zNorm = flatMesh->FaceNormal(face).Z();
			
			if (zNorm>0 && !m_bReflection)
			{
				glVertex2d(p1.X(), p1.Y());
				glVertex2d(p2.X(), p2.Y());
				glVertex2d(p3.X(), p3.Y());			
			}
			else if (zNorm<0 && m_bReflection)
			{
				glVertex2d(p1.X(), p1.Y());
				glVertex2d(p3.X(), p3.Y());
				glVertex2d(p2.X(), p2.Y());
			}
		}
		glEnd();
	}
	
	DrawSamples(params, renderingParams, surfaceName);	
}


void SurfaceMidEdgeConf::DrawConformal(ParamParser * params, 
									   const VKString & renderingParams,
									   const VKString & surfaceName)
{
	LoadConformalCamera();
	DrawConformalDoNotLoadCamera(params, renderingParams, surfaceName);
}

R3Point SurfaceMidEdgeConf::GetDrawablePosition(const SurfaceSample & sample, 
											ParamParser * params,
											const VKString & renderingParams,
											const VKString & surfaceName)
{
	if (params->GetStrValue(renderingParams, "ConfSurf", valid)=="RenderOriginal")
		return m_pWindow->PositionInMyCoords(m_iCameraID, 0, sample.GetPosition());
	else if (params->GetStrValue(renderingParams, "ConfSurf", valid)=="RenderConformal"
			 || params->GetStrValue(renderingParams, "ConfSurf", valid)=="RenderConformalMap1"
			 || params->GetStrValue(renderingParams, "ConfSurf", valid)=="RenderConformalMap2")
	{
		LinAlgComplex val=GetConfCoord(sample);
		//return R3Point(val.r, val.i, 0);
		return m_pWindow->PositionInMyCoords(m_iCameraID, 1, R3Point(val.r, val.i, 0));
	}
	else
		assert(false);
}





