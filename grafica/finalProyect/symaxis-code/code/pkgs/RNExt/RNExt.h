#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
// Sparse correspondences
class SparseVertexCorrespondence {
public:
  	R3Mesh *mesh[2];
  	RNArray<R3MeshVertex *> vertices[2];
  	int nvertices;
  	SparseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1) { 
    	mesh[0] = mesh0; 
    	mesh[1] = mesh1; 
    	nvertices = 0;
  	};
};

// Dense correspondences (map)
class DenseVertexCorrespondence {
public:
  	R3Mesh *mesh[2];
  	R3MeshVertex **vertices; 
  	DenseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1) { 
    	mesh[0] = mesh0; 
    	mesh[1] = mesh1; 
    	vertices = new R3MeshVertex * [ mesh0->NVertices() ];
    	for (int i = 0; i < mesh0->NVertices(); i++) vertices[i] = NULL;
  	};
};




// IO
R3Mesh * ReadMesh(char *filename);
R3MeshPropertySet *ReadProperties(R3Mesh *mesh, const char *filename);
R3MeshProperty *ReadProperty(R3Mesh *mesh, const char *filename);
RNArray<R3MeshVertex*>* ReadPoints(R3Mesh* mesh, const char* filename);
SparseVertexCorrespondence *ReadSparseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1, char *filename);
DenseVertexCorrespondence *ReadDenseVertexCorrespondence(R3Mesh* mesh0, R3Mesh* mesh1, char* filename);
int WriteSparseVertexCorrespondence(SparseVertexCorrespondence *correspondence, char *filename);
int WriteDenseVertexCorrespondence(DenseVertexCorrespondence *map, char *filename);
int WritePoints(R3Mesh* mesh, RNArray<R3MeshVertex*> *points, const char* filename);

// Algorithms
// Iterative furthest sampling
RNArray<R3MeshVertex*>*
SelectFurthestVertices(R3Mesh* mesh, const RNArray<R3MeshVertex*>& seeds, int npoints, double min_spacing = 0.0);
RNArray<R3MeshVertex*>*
SelectFurthestVertices(R3Mesh* mesh, int npoints, double min_spacing = 0.0);

// Nearest vertex search
R3MeshVertex* 
NearestNeighbor(R3Mesh* mesh, R3Point pos, R3MeshSearchTree* search_tree);

