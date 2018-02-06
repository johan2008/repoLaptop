#include "R3Shapes/R3Shapes.h"
#include "RNExt/RNExt.h"

struct Point {
  Point(void) 
    : vertex(NULL), position(0,0,0), normal(0,0,0) {};
  Point(const R3Point& position, const R3Vector& normal) 
    : vertex(NULL), position(position), normal(normal) {};
  Point(R3Mesh *mesh, R3MeshVertex *vertex)
    : vertex(vertex), position(mesh->VertexPosition(vertex)), normal(mesh->VertexNormal(vertex)) {};
  R3MeshVertex *vertex;
  R3Point position;
  R3Vector normal;
};

struct VertexData {
  R3MeshVertex *vertex;
  double distance;
  VertexData **heappointer;
};

// Find furthest vertex from a set of existing points
static R3MeshVertex *
FindFurthestVertex(R3Mesh *mesh, const RNArray<R3MeshVertex *>& seeds, VertexData *vertex_data)
{
  // Initialize vertex data
  for (int i = 0; i < mesh->NVertices(); i++) {
    VertexData *data = &vertex_data[i];
    data->distance = FLT_MAX;
    data->heappointer = NULL;
  }

  // Initialize heap
  VertexData tmp;
  RNHeap<VertexData *> heap(&tmp, &(tmp.distance), &(tmp.heappointer));
  for (int i = 0; i < seeds.NEntries(); i++) {
    R3MeshVertex *vertex = seeds.Kth(i);
    VertexData *data = &vertex_data[mesh->VertexID(vertex)];
    data->distance = 0;
    heap.Push(data);
  }

  // Iteratively pop off heap until find furthest vertex
  R3MeshVertex *vertex = NULL;
  while (!heap.IsEmpty()) {
    VertexData *data = heap.Pop();
    vertex = data->vertex;
    for (int j = 0; j < mesh->VertexValence(vertex); j++) {
      R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, j);
      R3MeshVertex *neighbor_vertex = mesh->VertexAcrossEdge(edge, vertex);
      VertexData *neighbor_data = (VertexData *) mesh->VertexData(neighbor_vertex);
      RNScalar old_distance = neighbor_data->distance;
      RNScalar new_distance = mesh->EdgeLength(edge) + data->distance;
      if (new_distance < old_distance) {
        neighbor_data->distance = new_distance;
        if (old_distance < FLT_MAX) heap.Update(neighbor_data);
        else heap.Push(neighbor_data);
      }
    }
  }

  // Return furthest vertex
  return vertex;
}



static RNArray<Point *> *
SelectFurthestPoints(R3Mesh *mesh, const RNArray<R3MeshVertex *>& seeds, int npoints, double min_spacing)
{
  // Check/adjust number of points
  if (npoints > mesh->NVertices()) npoints = mesh->NVertices();

  // Allocate array of points
  RNArray<Point *> *points = new RNArray<Point *>();
  if (!points) {
    fprintf(stderr, "Unable to allocate array of points\n");
    return NULL;
  }

  // Allocate vertex data
  VertexData *vertex_data = new VertexData [ mesh->NVertices() ];
  if (!vertex_data) {
    fprintf(stderr, "Unable to allocate vertex data\n");
    return NULL;
  }

  // Initialize vertex data
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    VertexData *data = &vertex_data[i];
    mesh->SetVertexData(vertex, data);
    data->vertex = vertex;
    data->distance = FLT_MAX;
    data->heappointer = NULL;
  }

  // Copy seeds 
  RNArray<R3MeshVertex *> vertices;
  for (int i = 0; i < seeds.NEntries(); i++) {
    R3MeshVertex *vertex = seeds.Kth(i);
    Point *point = new Point(mesh, vertex);
    points->Insert(point);
    vertices.Insert(vertex);
  }

  // Create random starting vertex, if none provided
  RNBoolean remove_first_vertex = FALSE;
  if (seeds.IsEmpty()) {
    int i = (int) (RNRandomScalar() * mesh->NVertices());
    R3MeshVertex *vertex = mesh->Vertex(i);
    vertices.Insert(vertex);
    remove_first_vertex = TRUE;
  }

  // Iteratively find furthest vertex
  while (vertices.NEntries() < npoints) {
    R3MeshVertex *vertex = FindFurthestVertex(mesh, vertices, vertex_data);
    if (!vertex) break;
    VertexData *data = &vertex_data[mesh->VertexID(vertex)];
    if (remove_first_vertex) { vertices.Truncate(0); remove_first_vertex = FALSE; }
    else if ((min_spacing > 0) && (data->distance < min_spacing)) break;
    Point *point = new Point(mesh, vertex);
    points->Insert(point);
    vertices.Insert(vertex);
  }

  // Delete vertex data
  delete [] vertex_data;

  // Return points
  return points;
}



/*static RNArray<Point *> **/
/*SelectFurthestPoints(R3Mesh *mesh, const RNArray<Point *>& seeds, int npoints, double min_spacing)*/
/*{*/
  /*// Seeed with vertices closest to points*/
  /*RNArray<R3MeshVertex *> seed_vertices;*/
  /*for (int i = 0; i < seeds.NEntries(); i++) {*/
    /*R3MeshVertex *vertex = seeds[i]->vertex;*/
    /*if (vertex) seed_vertices.Insert(vertex);*/
    /*else { fprintf(stderr, "Not implemented\n"); return NULL; }*/
  /*}*/

  /*// Select furthest points*/
  /*return SelectFurthestPoints(mesh, seed_vertices, npoints, min_spacing);*/
/*} */



static RNArray<Point *> *
SelectFurthestPoints(R3Mesh *mesh, int npoints, double min_spacing)
{
  // Select furthest points from scratch
  RNArray<R3MeshVertex *> seeds;
  return SelectFurthestPoints(mesh, seeds, npoints, min_spacing);
}


RNArray<R3MeshVertex*>*
SelectFurthestVertices(R3Mesh* mesh, const RNArray<R3MeshVertex*>& seeds, int npoints, double min_spacing)
{
	RNArray<Point *>* points = SelectFurthestPoints(mesh, seeds, npoints, min_spacing);
	RNArray<R3MeshVertex*>* vertices = new RNArray<R3MeshVertex*>;
	for (int i=0; i<points->NEntries(); i++)
	{
		R3MeshVertex* vertex = points->Kth(i)->vertex;
		vertices->Insert(vertex);
	}
    for (int i=0; i<points->NEntries(); i++)
	{
		delete points->Kth(i);
	}
	delete points;
	return vertices;
}

RNArray<R3MeshVertex*>*
SelectFurthestVertices(R3Mesh* mesh, int npoints, double min_spacing)
{
	RNArray<Point *>* points = SelectFurthestPoints(mesh, npoints, min_spacing);
	RNArray<R3MeshVertex*>* vertices = new RNArray<R3MeshVertex*>;
	for (int i=0; i<points->NEntries(); i++)
	{
		R3MeshVertex* vertex = points->Kth(i)->vertex;
		vertices->Insert(vertex);
	}
    for (int i=0; i<points->NEntries(); i++)
	{
		delete points->Kth(i);
	}
	delete points;
	return vertices;
}
