#include "Region.h"
#include "RegionMerger.h"
Region::
Region(void)
	:m_RegionMerger(NULL)
{}

Region::
Region(RegionMerger* r, R3MeshVertex* v)
	:m_RegionMerger(r)
{
	R3Mesh* mesh = r->GetMesh();
	m_VertexArray.push_back(v);
	for (int i=0; i<mesh->VertexValence(v); i++)
	{
		R3MeshEdge* e = mesh->EdgeOnVertex(v, i);
		m_EdgeArray.push_back(e);
	}
	// sort edges using bubble sort
	for (int i=0; i<NEdges()-1; i++)
	{
		for (int j=0; j<NEdges()-i-1; j++)
		{
			R3MeshEdge* e[2];
			e[0] = m_EdgeArray[j];
			e[1] = m_EdgeArray[j+1];
			if (m_RegionMerger->GetEdgeValue(e[0]) < m_RegionMerger->GetEdgeValue(e[1]))
			{
				R3MeshEdge* temp = e[1];
				m_EdgeArray[j+1] = e[0];
				m_EdgeArray[j] = temp;
			}
		}
	}
}

Region::
Region(Region* r0, Region* r1)
{
	assert(r0->m_RegionMerger == r1->m_RegionMerger);
	m_RegionMerger = r0->m_RegionMerger;
   /* MergeEdgeArrays(r0.m_EdgeArray, r1.m_EdgeArray);*/
	/*MergeVertexArrays(r0.m_VertexArray, r1.m_VertexArray);*/
	MergeEdgeArrays(r0, r1);
	MergeVertexArrays(r0, r1);

}

RNScalar Region::
GetPersistence(void) const
{
	RNScalar val = GetCentralVertexValue() - GetCriticalEdgeValue();
	if (val < 0.0) val = 0.0;
	return val;	
}

R3MeshVertex* Region::
GetCentralVertex(void) const
{
	assert(NVertices() > 0);
	return m_VertexArray[0];
}

R3MeshEdge* Region::
GetCriticalEdge(void) const
{
	if (NEdges() > 0)
		return m_EdgeArray[0];
	else
		return NULL;
}

RNScalar Region::
GetCentralVertexValue(void) const
{
	R3MeshVertex* v = GetCentralVertex();
	assert(m_RegionMerger);
	return m_RegionMerger->GetVertexValue(v);
}

RNScalar Region::
GetCriticalEdgeValue(void) const
{
	R3MeshEdge* e = GetCriticalEdge();
	assert(m_RegionMerger);
	if (!e)
	{
		return 0.0;
	}
	else
	{
		return m_RegionMerger->GetEdgeValue(e);
	}
}

R3MeshVertex* Region::
GetVertex(int idx) const
{
	assert(idx>=0 && idx<NVertices());
	return m_VertexArray[idx];
}

R3MeshEdge* Region::
GetEdge(int idx) const
{
	assert(idx>=0 && idx<NEdges());
	return m_EdgeArray[idx];
}

int Region::
NVertices(void) const
{
	return int(m_VertexArray.size());
}

int Region::
NEdges(void) const
{
	return int(m_EdgeArray.size());
}

int Region::
MergeEdgeArrays(Region* r0, Region* r1)
{
	std::vector<R3MeshEdge*>& array0 = r0->m_EdgeArray;
	std::vector<R3MeshEdge*>& array1 = r1->m_EdgeArray;
	assert(m_RegionMerger);
	m_EdgeArray.clear();
	int p[2] = {0, 0};
	while (p[0] < int(array0.size()) && p[1] < int(array1.size()))
	{
		R3MeshEdge* candidate = NULL;
		if (m_RegionMerger->GetEdgeValue(array0[p[0]]) < m_RegionMerger->GetEdgeValue(array1[p[1]]))
		{
			candidate = array1[p[1]];
			p[1] ++;
		}
		else
		{
			candidate = array0[p[0]];
			p[0] ++;
		}
		if (!m_RegionMerger->IsEdgeOnBorder(candidate, r0, r1))
		{
			m_EdgeArray.push_back(candidate);
		}
	}
	while (p[0] < int(array0.size()))
	{
		R3MeshEdge* candidate = array0[p[0]];
		if (!m_RegionMerger->IsEdgeOnBorder(candidate, r0, r1))
			m_EdgeArray.push_back(candidate);
		p[0] ++;
	}
	while (p[1] < int(array1.size()))
	{
		R3MeshEdge* candidate = array1[p[1]];
		if (!m_RegionMerger->IsEdgeOnBorder(candidate, r0, r1))
			m_EdgeArray.push_back(candidate);
		p[1] ++;
	}
	return 1;
}

int Region::
MergeVertexArrays(Region* r0, Region* r1)
{
	std::vector<R3MeshVertex*>& array0 = r0->m_VertexArray;
	std::vector<R3MeshVertex*>& array1 = r1->m_VertexArray;
	assert(m_RegionMerger);
	m_VertexArray.clear();
	int p[2] = {0, 0};
	while (p[0] < int(array0.size()) && p[1] < int(array1.size()))
	{
		assert (array0[p[0]] != array1[p[1]]);
		if (m_RegionMerger->GetVertexValue(array0[p[0]]) < m_RegionMerger->GetVertexValue(array1[p[1]]))
		{
			m_VertexArray.push_back(array1[p[1]]);
			p[1] ++;
		}
		else
		{
			m_VertexArray.push_back(array0[p[0]]);
			p[0] ++;
		}
	}
	while (p[0] < int(array0.size()))
	{
		m_VertexArray.push_back(array0[p[0]]);
		p[0] ++;
	}
	while (p[1] < int(array1.size()))
	{
		m_VertexArray.push_back(array1[p[1]]);
		p[1] ++;
	}
	return 1;
}
