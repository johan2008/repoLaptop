#ifndef __MY_ANNMesh_H
#define __MY_ANNMesh_H
#include "R3Shapes/R3Shapes.h"
#include "ANN.h"
#include <vector>
using namespace std;
class ANNMesh
{
private:
	R3Mesh* m_mesh;
	ANNpointArray m_data_pts;
	static const int dim = 3;
	int m_num_pts;
	ANNkd_tree* m_kdtree;
public:
	ANNMesh();
	~ANNMesh();
	ANNMesh(R3Mesh* mesh);
	R3Mesh* Mesh();
	int Initialize(R3Mesh* mesh);
	RNArray<R3MeshVertex*>* SearchNearestVertices(R3Point& pts, int k);
	RNArray<R3MeshVertex*>* SearchNearestVertices(R3MeshVertex* vts, int k);
	vector<R3Point>* SearchNearestPoints(R3Point& pts, int k);
	vector<R3Point>* SearchNearestPoints(R3MeshVertex* vts, int k);
};
#endif
