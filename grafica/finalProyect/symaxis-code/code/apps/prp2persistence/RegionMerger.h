#ifndef __MY_REGIONMERGER_H
#define __MY_REGIONMERGER_H
#include "R3Shapes/R3Shapes.h"
#include <vector>
#include "BasicDataStructure/BasicHeap.h"
#include "Region.h"
class RegionMerger
{
public:
	RegionMerger(void);
	RegionMerger(R3Mesh* mesh, R3MeshProperty* property, RNScalar lambda = FLT_MAX);
	~RegionMerger(void);
	R3Mesh* GetMesh() const;
public:
	// get m_ValueOnVertices
	RNScalar GetVertexValue(R3MeshVertex* v) const;
	// get m_ValueOnEdges
	RNScalar GetEdgeValue(R3MeshEdge* e) const;
	// get m_RegionOnVertices
	Region* GetVertexRegion(R3MeshVertex* v) const;
	// update m_Persistence at the central vertex
	int UpdateVertexPersistence(Region* r);
	// update the region belongings
	int UpdateVertexRegion(Region* r);
	// decide whether a edge is across two regions
	int IsEdgeOnBorder(R3MeshEdge* e, Region* r0, Region* r1) const;
private:
	R3Mesh* m_Mesh;
	// property value on vertices
	std::vector<RNScalar> m_ValueOnVertices;
	// property value on edges
	std::vector<RNScalar> m_ValueOnEdges;
	// belonging of regions
	std::vector<Region*> m_RegionOnVertices;
	// current persistency on each vertex
	// only update in function UpdatePersistence
	std::vector<RNScalar> m_Persistence;
	// priority heap of regions
	BasicHeap<Region> m_HeapRegions;
	RNScalar m_Lambda;
private:
	// utility
	int ComputeValueFromProperty(R3MeshProperty* property); 
public:
	// solve
	int Run(void);
	// return critical vertices
	int GetCentralVertices(std::vector<R3MeshVertex*>& centrals, std::vector<RNScalar>& persistence) const;
	// return persistence as property
	R3MeshProperty* GetPersistence(void) const;
	// return the assignment of regions
	int GetRegionAssignment(std::vector<int>& assign);
};
#endif
