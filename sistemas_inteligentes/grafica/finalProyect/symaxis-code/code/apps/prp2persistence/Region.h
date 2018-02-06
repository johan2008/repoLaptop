#ifndef __MY_REGION_H
#define __MY_REGION_H
#include "BasicDataStructure/BasicHeap.h"
#include "R3Shapes/R3Shapes.h"
#include <vector>
class RegionMerger;
class Region
{
public:
	Region(void);
	// generate from single vertex
	Region(RegionMerger* r, R3MeshVertex* v);
	// generate from merging two regions
	Region(Region* r0, Region* r1);
private:
	// associated mesh
	RegionMerger* m_RegionMerger;
	// array of vertices (ordered by value, from highest to lowest)
	std::vector<R3MeshVertex*> m_VertexArray;
	// array of edges (ordered by value, from highest to lowest)
	std::vector<R3MeshEdge*> m_EdgeArray;
public:
	// get the persistence value of the region
	RNScalar GetPersistence(void) const;
	// get central vertex
	R3MeshVertex* GetCentralVertex(void) const;
	// get critical edge 
	R3MeshEdge* GetCriticalEdge(void) const;
	// get value of central vertex
	RNScalar GetCentralVertexValue(void) const;
	// get value of critical edge
	RNScalar GetCriticalEdgeValue(void) const;
	// get general vertex
	R3MeshVertex* GetVertex(int idx) const;
	// get general edge
	R3MeshEdge* GetEdge(int idx) const;
public:
	// utility
	int NVertices(void) const;
	int NEdges(void) const;
	// merge ordered arrays of boundary edges
	int MergeEdgeArrays(Region* r0, Region* r1);
	int MergeVertexArrays(Region* r0, Region* r1);
};

inline bool
operator==(const Region& lhs, const Region& rhs)
{
	return (lhs.GetPersistence() == rhs.GetPersistence() && lhs.GetCentralVertexValue() == rhs.GetCentralVertexValue());
}

inline bool
operator!=(const Region& lhs, const Region& rhs)
{
	return !operator==(lhs, rhs);
}

// we would like to pick the region with minimum persistence first
// if two regions have same persistence, we pick the one with smalller central value
inline bool
operator< (const Region& lhs, const Region& rhs)
{
	if (lhs.GetPersistence() == rhs.GetPersistence())
		return lhs.GetCentralVertexValue() < rhs.GetCentralVertexValue();
	return lhs.GetPersistence() < rhs.GetPersistence();
}

inline bool
operator> (const Region& lhs, const Region& rhs)
{
	return operator< (rhs, lhs);
}

inline bool
operator<= (const Region& lhs, const Region& rhs)
{
	return !operator> (lhs, rhs);
}

inline bool
operator>= (const Region& lhs, const Region& rhs)
{
	return !operator< (lhs, rhs);
}

#endif
