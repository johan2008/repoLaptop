#include "SurfaceSample.h"
#include <iostream>
#include "VkFunctions.h"
#include "MapFlattening.h"

SurfaceSample::SurfaceSample()
{
	m_pMesh=NULL;
	for (int i=0; i<3; i++)
		m_fBarycentricCoords[i] = -1;
	m_iTriangleID = -1;
	m_sType="Invalid";
}

SurfaceSample::SurfaceSample(R3MeshIntersection & closest, R3Point query, R3Mesh * mesh,
							 const VKString & type)
{	
	assert (!InitializedSpecialMesh(mesh, 0, 0, 0));
	assert(closest.type!=R3_MESH_NULL_TYPE);
	m_sType = type;	
	m_pMesh = mesh;
	m_iTriangleID = m_pMesh->FaceID(closest.face);
	assert(m_iTriangleID>=0);
	assert(m_iTriangleID < m_pMesh->NFaces());
	
	// get barycentric coords of triangle
	R3Point baryCoords = mesh->FaceBarycentric(m_pMesh->Face(m_iTriangleID), query);
	
	if (baryCoords.X()<=0)
		baryCoords.X(0);
	if (baryCoords.Y()<=0)
		baryCoords.Y(0);
	if (baryCoords.Z()<=0)
		baryCoords.Z(0);
	double baryNorm = baryCoords.X()+baryCoords.Y()+baryCoords.Z();
	if (baryNorm>=0)
	{
		baryCoords.X(baryCoords.X()/baryNorm);
		baryCoords.Y(baryCoords.Y()/baryNorm);
		baryCoords.Z(baryCoords.Z()/baryNorm);
	}
	else
		baryCoords.X(1.);
	m_fBarycentricCoords[0] = baryCoords.X();
	m_fBarycentricCoords[1] = baryCoords.Y();
	m_fBarycentricCoords[2] = baryCoords.Z();
}

SurfaceSample::SurfaceSample(const SurfaceSample & other)
{
	if (InitializedSpecialMesh(other.m_pMesh, other.m_fBarycentricCoords[0], 
							   other.m_fBarycentricCoords[1], other.m_fBarycentricCoords[2]))
		return;
	m_pMesh = other.m_pMesh;
	for (int i=0; i<3; i++)
		m_fBarycentricCoords[i] = other.m_fBarycentricCoords[i];
	m_iTriangleID = other.m_iTriangleID;
	m_sType = other.m_sType;
}

void SurfaceSample::operator = (const SurfaceSample & other)
{
	if (InitializedSpecialMesh(other.m_pMesh, other.m_fBarycentricCoords[0], 
							   other.m_fBarycentricCoords[1], other.m_fBarycentricCoords[2]))
		return;
	
	m_pMesh = other.m_pMesh;
	for (int i=0; i<3; i++)
		m_fBarycentricCoords[i] = other.m_fBarycentricCoords[i];
	m_iTriangleID = other.m_iTriangleID;
	m_sType = other.m_sType;	
}

SurfaceSample::SurfaceSample(int triangleID, double b1, double b2, double b3, 
							 R3Mesh * mesh,  const VKString & type)
{
	if (InitializedSpecialMesh(mesh, b1, b2, b3))
		return;

	m_sType = type;
	m_fBarycentricCoords[0] = b1;
	m_fBarycentricCoords[1] = b2;
	m_fBarycentricCoords[2] = b3;
	assert(m_fBarycentricCoords[0]>=0 && m_fBarycentricCoords[1]>=0 && m_fBarycentricCoords[2]>=0);
	if (vkAbs(m_fBarycentricCoords[0]+m_fBarycentricCoords[1]+m_fBarycentricCoords[2]-1.)>0.0001)
	{
		std::cout<<"[ERROR] Barycentric coords do not interpolate. ";
		std::cout<<(m_fBarycentricCoords[0]+m_fBarycentricCoords[1]+m_fBarycentricCoords[2])<<std::endl;
		assert(vkAbs(m_fBarycentricCoords[0]+m_fBarycentricCoords[1]+m_fBarycentricCoords[2]-1.)<=0.0001);	
	}
	m_iTriangleID = triangleID;
	m_pMesh = mesh;
}

SurfaceSample::SurfaceSample(int triangleID, R3Point posOnTriangle, 
							 R3Mesh * mesh,  const VKString & type)
{
	assert (!InitializedSpecialMesh(mesh, 0,0,0));

	m_sType = type;
	m_iTriangleID = triangleID;
	m_pMesh = mesh;
	R3Point bary = m_pMesh->FaceBarycentric(m_pMesh->Face(m_iTriangleID), posOnTriangle);
	m_fBarycentricCoords[0] = bary.X();
	m_fBarycentricCoords[1] = bary.Y();
	m_fBarycentricCoords[2] = bary.Z();
	assert(m_fBarycentricCoords[0]>=0 && m_fBarycentricCoords[1]>=0 && m_fBarycentricCoords[2]>=0);
	assert(vkAbs(m_fBarycentricCoords[0]+m_fBarycentricCoords[1]+m_fBarycentricCoords[2]-1.)<=0.0001);	
}

SurfaceSample::SurfaceSample(int vertexID, R3Mesh * mesh, const VKString & type)
{
	assert (!InitializedSpecialMesh(mesh, 0,0,0));
	
	m_sType = type;
	m_pMesh = mesh;
	assert(m_pMesh->VertexValence(m_pMesh->Vertex(vertexID))>0);
	R3MeshFace * face = m_pMesh->FaceOnVertex(m_pMesh->Vertex(vertexID));
	if (face==NULL)
	{
		for (int i=0; i<m_pMesh->VertexValence(m_pMesh->Vertex(vertexID)); i++)
		{
			face = m_pMesh->FaceOnEdge(m_pMesh->EdgeOnVertex(m_pMesh->Vertex(vertexID), i), 0);
			if (face==NULL)
				face = m_pMesh->FaceOnEdge(m_pMesh->EdgeOnVertex(m_pMesh->Vertex(vertexID), i), 1);
			if (face!=NULL)
				break;
		}		
		if (face==NULL)
			std::cout<<"\n[ERROR] Vertex "<<vertexID<<" is not connected to a face."<<std::endl;
		assert(face!=NULL);
		
	}

	m_iTriangleID = m_pMesh->FaceID(face);
	for (int i=0; i<3; i++)
	{
	    if (vertexID==m_pMesh->VertexID(m_pMesh->VertexOnFace(face, i)))	
			m_fBarycentricCoords[i] = 1;
	    else
			m_fBarycentricCoords[i] = 0;
	}
	assert(m_fBarycentricCoords[0]>=0 && m_fBarycentricCoords[1]>=0 && m_fBarycentricCoords[2]>=0);
	assert(vkAbs(m_fBarycentricCoords[0]+m_fBarycentricCoords[1]+m_fBarycentricCoords[2]-1.)<=0.0001);		
}

bool SurfaceSample::Invalid() const
{
	return (m_iTriangleID<0 || m_pMesh==NULL)
			&& (m_pMesh!=Surface2DPlane::m_pInfinitePlanePseudomesh);
}

int SurfaceSample::TriID() const
{
	return m_iTriangleID;
}

double SurfaceSample::B(int coordID) const
{
	return m_fBarycentricCoords[coordID];
}

int SurfaceSample::VertID(int coordID) const
{
	return m_pMesh->VertexID(m_pMesh->VertexOnFace(m_pMesh->Face(m_iTriangleID), coordID));
}

double SurfaceSample::Interpolate(double b1, double b2, double b3) const
{
	return b1 * B(0) + b2 * B(1) + b3 * B(2);
}

int SurfaceSample::NearestVertex(bool * atVertex) const
{
	if (atVertex!=NULL)
		*atVertex = (B(0)==1) || (B(1)==1) || (B(2)==1);

	if (B(0)>=B(1) && B(0)>=B(2))
		return VertID(0);
	else if (B(1)>=B(0) && B(1)>=B(2))
		return VertID(1);
	else
		return VertID(2);
}
	
R3Point SurfaceSample::GetPosition() const
{
	if (m_pMesh==Surface2DPlane::m_pInfinitePlanePseudomesh)
		return R3Point(m_fBarycentricCoords[0], m_fBarycentricCoords[1], 0);
	
	assert(m_pMesh!=NULL);
	assert(m_iTriangleID < m_pMesh->NFaces());
	assert(m_iTriangleID>=0);
	RNMagnitude mag[3];
	for (int i=0; i<3; i++)
		mag[i] = m_fBarycentricCoords[i]; 
	return m_pMesh->FacePoint(m_pMesh->Face(m_iTriangleID), mag);
}

R3Vector SurfaceSample::Normal() const
{
	if (m_pMesh==Surface2DPlane::m_pInfinitePlanePseudomesh)
		return R3Vector(0,0,1);
	
	return m_pMesh->FaceNormal(m_pMesh->Face(m_iTriangleID));
}
 
VKString SurfaceSample::GetSampleType() const
{
	return m_sType;
}

void SurfaceSample::SetSampleType(const VKString & type)
{
	m_sType = type;
}

bool SurfaceSample::operator == (const SurfaceSample & other) const
{
	bool exactlySame = (m_iTriangleID==other.m_iTriangleID);
	exactlySame = exactlySame && (m_fBarycentricCoords[0]==other.m_fBarycentricCoords[0]);
	exactlySame = exactlySame && (m_fBarycentricCoords[1]==other.m_fBarycentricCoords[1]);
	exactlySame = exactlySame && (m_fBarycentricCoords[2]==other.m_fBarycentricCoords[2]);
	return (exactlySame || GetPosition()==other.GetPosition()) && m_pMesh==other.m_pMesh;	
}

bool SurfaceSample::CheckMeshIsSame(const R3Mesh * mesh) const
{
	return mesh==m_pMesh;
}

bool SurfaceSample::InitializedSpecialMesh(R3Mesh * mesh, double b1, double b2, double b3)
{
	if (mesh==Surface2DPlane::m_pInfinitePlanePseudomesh)
	{
		m_fBarycentricCoords[0] = b1;
		m_fBarycentricCoords[1] = b2;		
		m_fBarycentricCoords[2] = 0;		
		m_iTriangleID = -1;
		m_pMesh = mesh;
		m_sType = "FlatSample";
		return true;
	}
	return false;
}




