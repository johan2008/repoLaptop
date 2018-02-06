#include "MapEuclidean.h"

MapEuclidean::MapEuclidean(SampledSurface * M1, SampledSurface * M2)
: SurfaceMap(M1, M2)
{
	m_pSearchTreeOnM1 = new R3MeshSearchTree(M1->GetMesh());
	m_pSearchTreeOnM2 = new R3MeshSearchTree(M2->GetMesh());
}
	
SurfaceSample MapEuclidean::ForwardMap(const SurfaceSample & s)
{
	return MapEuclideanDist(s, m_pSearchTreeOnM2, m_pM2->GetMesh());
}

SurfaceSample MapEuclidean::InverseMap(const SurfaceSample & s)
{
	return MapEuclideanDist(s, m_pSearchTreeOnM1, m_pM1->GetMesh());
}


SurfaceSample MapEuclidean::MapEuclideanDist(const SurfaceSample & s, 
											 R3MeshSearchTree * nearestOnMesh, 
											 R3Mesh * corrMesh)
{
	R3Vector normal = s.Normal();
	
	R3MeshIntersection closest;
	nearestOnMesh->FindClosest(s.GetPosition(), normal, closest, 0, FLT_MAX);
	
	assert(closest.type!=R3_MESH_NULL_TYPE);
	
//	if (error!=NULL)
//		*error = (closest.point-s.GetPosition()).Length();
	
	return SurfaceSample(corrMesh->FaceID(closest.face), closest.point, corrMesh);
}

void MapEuclidean::SaveMap(std::ofstream & textStream)
{
	textStream<<"Format MapEuclidean\n";
}

void MapEuclidean::LoadMap(std::ifstream & textStream)
{
	// nothing here yet: eventually maybe we should support rigid transformations
	std::string tempStr;
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="Format");
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="MapEuclidean");	
}

VKString MapEuclidean::GetSurfaceMapType()
{
	return "MapEuclidean";
}


