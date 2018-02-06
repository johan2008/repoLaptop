#include "ANN/ANNMesh.h"
// default constructor
ANNMesh::ANNMesh()
{
	m_mesh = NULL;
	m_kdtree = NULL;
}	

// default xigou
ANNMesh::~ANNMesh()
{
	delete m_kdtree;
	annClose();
	// deallocate an array
	//annDeallocPts(m_data_pts);
}

// constructor
ANNMesh::ANNMesh(R3Mesh* mesh)
{
	int suc = Initialize(mesh);
	if (!suc)
	{
		fprintf(stderr, "Unable to import mesh to ANNMesh\n");
	}
}

int ANNMesh::Initialize(R3Mesh* mesh)
{
	m_mesh = mesh;
	m_num_pts = mesh->NVertices();
	m_data_pts = annAllocPts(m_num_pts, dim);
	// set values for m-data-pts
	for (int i=0; i<m_num_pts; i++)
	{
		R3MeshVertex* v = mesh->Vertex(i);
		R3Point pos = mesh->VertexPosition(v);
		for (int j=0; j<dim; j++)
		{
			m_data_pts[i][j] = pos[j];
		}
	}
	// build kd-tree
	m_kdtree = new ANNkd_tree(m_data_pts, m_num_pts, dim);
	return 1;
}

R3Mesh* ANNMesh::Mesh()
{
	return m_mesh;
}

RNArray<R3MeshVertex*>*
ANNMesh::SearchNearestVertices(R3Point& pts, int k)
{
	ANNpoint queryPt = annAllocPt(dim);
	queryPt[0] = pts[0];
	queryPt[1] = pts[1];
	queryPt[2] = pts[2];
	ANNidxArray nnIdx = new ANNidx[k];
	ANNdistArray dists = new ANNdist[k];
	m_kdtree->annkSearch(queryPt, k, nnIdx, dists, 0.0);
	RNArray<R3MeshVertex*>* vertices = new RNArray<R3MeshVertex*>;
	for (int i=0; i<k; i++)
	{
		R3MeshVertex* v = m_mesh->Vertex(nnIdx[i]);
		vertices->Insert(v);
	}
	delete[] nnIdx;
	delete[] dists;
	annDeallocPt(queryPt);
	return vertices;
}

RNArray<R3MeshVertex*>*
ANNMesh::SearchNearestVertices(R3MeshVertex* vts, int k)
{
	R3Point pts = m_mesh->VertexPosition(vts);
	return SearchNearestVertices(pts, k);
}

vector<R3Point>*
ANNMesh::SearchNearestPoints(R3Point& pts, int k)
{
	RNArray<R3MeshVertex*>* vertices = SearchNearestVertices(pts, k);
	vector<R3Point>* points = new vector<R3Point>;
	for (int i=0; i<vertices->NEntries(); i++)
	{
		R3MeshVertex* v = vertices->Kth(i);
		points->push_back(m_mesh->VertexPosition(v));
	}
	delete vertices;
	return points;
}

vector<R3Point>*
ANNMesh::SearchNearestPoints(R3MeshVertex* vts, int k)
{
	R3Point pts = m_mesh->VertexPosition(vts);
	return SearchNearestPoints(pts, k);
}
