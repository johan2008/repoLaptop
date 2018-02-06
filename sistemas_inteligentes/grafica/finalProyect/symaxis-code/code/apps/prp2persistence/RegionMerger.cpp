#include "RegionMerger.h"
RegionMerger::
RegionMerger(void)
	:m_Mesh(NULL), m_Lambda(FLT_MAX)
{}

RegionMerger::
RegionMerger(R3Mesh* mesh, R3MeshProperty* property, RNScalar lambda)
	:m_Mesh(mesh), m_Lambda(lambda)
{
	ComputeValueFromProperty(property);
}

RegionMerger::
~RegionMerger(void)
{
	// delete regions that are still in the heap
	for (int i=0; i<m_HeapRegions.NEntries(); i++)
	{
		Region* r = m_HeapRegions.Kth(i);
		delete r;
	}
}

R3Mesh* RegionMerger::
GetMesh() const
{
	return m_Mesh;
}

RNScalar RegionMerger::
GetVertexValue(R3MeshVertex* v) const
{
	assert(m_Mesh);
	int idx = m_Mesh->VertexID(v);
	assert(idx < int(m_ValueOnVertices.size()));
	return m_ValueOnVertices[idx];
}

RNScalar RegionMerger::
GetEdgeValue(R3MeshEdge* e) const
{
	assert(m_Mesh);
	int idx = m_Mesh->EdgeID(e);
	assert(idx < int(m_ValueOnEdges.size()));
	return m_ValueOnEdges[idx];
}

Region* RegionMerger::
GetVertexRegion(R3MeshVertex* v) const
{
	assert(m_Mesh);
	int idx = m_Mesh->VertexID(v);
	assert(idx < int(m_RegionOnVertices.size()));
	return m_RegionOnVertices[idx];
}

int RegionMerger::
UpdateVertexPersistence(Region* r)
{
	R3MeshVertex* central = r->GetCentralVertex();
	m_Persistence[m_Mesh->VertexID(central)] = r->GetCentralVertexValue();
	return 1;
}

int RegionMerger::
UpdateVertexRegion(Region* r)
{
	for (int i=0; i<r->NVertices(); i++)
	{
		R3MeshVertex* v = r->GetVertex(i);
		int idx = m_Mesh->VertexID(v);
		m_RegionOnVertices[idx] = r;
	}
	return 1;
}

int RegionMerger::
IsEdgeOnBorder(R3MeshEdge* e, Region* r0, Region* r1) const
{
	R3MeshVertex* v[2];
	v[0] = m_Mesh->VertexOnEdge(e, 0);
	v[1] = m_Mesh->VertexOnEdge(e, 1);
	assert(v[0] && v[1]);
	int idx[2];
	idx[0] = m_Mesh->VertexID(v[0]);
	idx[1] = m_Mesh->VertexID(v[1]);
	if (m_RegionOnVertices[idx[0]] == r0 && m_RegionOnVertices[idx[1]] == r1)
	{
		return 1;
	}
	if (m_RegionOnVertices[idx[0]] == r1 && m_RegionOnVertices[idx[1]] == r0)
	{
		return 1;
	}
	return 0;
}

int RegionMerger::
ComputeValueFromProperty(R3MeshProperty * property)
{
	assert(m_Mesh);
	m_ValueOnVertices.resize(m_Mesh->NVertices());
	m_ValueOnEdges.resize(m_Mesh->NEdges());
	m_RegionOnVertices.resize(m_Mesh->NVertices());
	m_Persistence.resize(m_Mesh->NVertices());
	for (int i=0; i<m_Mesh->NVertices(); i++)
	{
		m_ValueOnVertices[i] = property->VertexValue(i);
	}
	for (int i=0; i<m_Mesh->NEdges(); i++)
	{
		R3MeshEdge* edge = m_Mesh->Edge(i);
		R3MeshVertex* v[2];
		v[0] = m_Mesh->VertexOnEdge(edge, 0);
		v[1] = m_Mesh->VertexOnEdge(edge, 1);
		assert(v[0] && v[1]);
		m_ValueOnEdges[i] = (property->VertexValue(v[0]) + property->VertexValue(v[1])) * 0.5;
	}
	for (int i=0; i<m_Mesh->NVertices(); i++)
	{
		m_RegionOnVertices[i] = NULL;
	}
	return 1;
}

int RegionMerger::
Run(void)
{
	// initialize - each vertex is a region
	printf("Initializing the regions ...\n");
	std::vector<Region*> init_region_list;
	for (int i=0; i<m_Mesh->NVertices(); i++)
	{
		R3MeshVertex* v = m_Mesh->Vertex(i);
		Region* region = new Region(this, v);
		// add to init_region_list
		init_region_list.push_back(region);
		// update vertex region
		UpdateVertexRegion(region);
		// update persistence
		UpdateVertexPersistence(region);
	}
	// make heap
	printf("    initialzing the heap ... \n");
	m_HeapRegions.Mount(init_region_list, 1);

	// repeating ...
	// (1) pop the region r with least persistency
	// (2) find out the neighbor r_neighbor
	// (3) remove r and r_neighbor from the heap, create r_merged = r \cup r_neighbor
	// (4) udpate the data in 'regionmerger'
	printf("    repeatedly merging (lambda = %.3f) ...\n", m_Lambda);
	while (true)
	{
		// (1) pop the region r with least persistency
		Region* r = m_HeapRegions.Peek();
		if (r->GetPersistence() > m_Lambda)
		{
			break;
		}
		if (r->GetPersistence() > 0.003)
		{
			printf("......%d/%d (%.3f)\n", m_HeapRegions.NEntries(), m_Mesh->NVertices(), r->GetPersistence());
		}
		// (2) find out the neighbor r_neighbor
		R3MeshEdge* critical_edge = r->GetCriticalEdge();
		if (!critical_edge)
		{
			printf("no connection, have to stop here!\n");
			break;
		}
		R3MeshVertex* v[2];
		v[0] = m_Mesh->VertexOnEdge(critical_edge, 0);
		v[1] = m_Mesh->VertexOnEdge(critical_edge, 1);
	    assert(GetVertexRegion(v[0]) == r ^ GetVertexRegion(v[1]) == r);
		Region* r_neighbor = NULL;
		if (GetVertexRegion(v[0]) == r)
		{
			r_neighbor = GetVertexRegion(v[1]);
		}
		else if (GetVertexRegion(v[1]) == r)
		{
			r_neighbor = GetVertexRegion(v[0]);
		}
		// (3) remove r and r_neighbor from the heap, create r_merged = r \cup r_neighbor
		Region* r_new = new Region(r, r_neighbor);
		m_HeapRegions.Pop();
		m_HeapRegions.Remove(r_neighbor);
		m_HeapRegions.Add(r_new);
		// (4) udpate the data in 'regionmerger'
		UpdateVertexRegion(r_new);
		UpdateVertexPersistence(r_new);
		// (5) clean up
		delete r_neighbor;
		delete r;
	}
	return 1;
}

int RegionMerger::
GetCentralVertices(std::vector<R3MeshVertex*>& centrals, std::vector<RNScalar>& persistence) const
{
	centrals.clear();
	persistence.clear();
	for (int i=0; i<m_HeapRegions.NEntries(); i++)
	{
		Region* r = m_HeapRegions.Kth(i);
		R3MeshVertex* v = r->GetCentralVertex();
		RNScalar val = r->GetCentralVertexValue();
		centrals.push_back(v);
		persistence.push_back(val);
	}
	return 1;
}

R3MeshProperty* RegionMerger::
GetPersistence(void) const
{
	assert(m_Mesh);
	assert(int(m_Persistence.size()) == m_Mesh->NVertices());
	R3MeshProperty* prp = new R3MeshProperty(m_Mesh);
	for (int i=0; i<m_Mesh->NVertices(); i++)
	{
		prp->SetVertexValue(i, m_Persistence[i]);
	}
	return prp;
}

int RegionMerger::
GetRegionAssignment(std::vector<int>& assign)
{
	assert(m_Mesh);
	assert(int(m_RegionOnVertices.size()) == m_Mesh->NVertices());
	assign.clear();
	for (int i=0; i<m_Mesh->NVertices(); i++)
	{
		Region* r = m_RegionOnVertices[i];
		assign.push_back(m_HeapRegions.GetIndex(r));
	}
	return 1;
}
