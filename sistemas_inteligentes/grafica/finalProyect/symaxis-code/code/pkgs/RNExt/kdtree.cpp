#include "R3Shapes/R3Shapes.h"
#include "RNExt/RNExt.h"

R3MeshVertex* 
NearestNeighbor(R3Mesh* mesh, R3Point pos, R3MeshSearchTree* search_tree)
{
	R3MeshIntersection closest;
	search_tree->FindClosest(pos, closest, 0, 100.0);
	if (closest.type == R3_MESH_NULL_TYPE)
	{
		return NULL;
	}

//	RNLength closest_distance_to_surface = closest.t;
	R3Point closest_point_on_surface = closest.point;

	// Determine closest vertex
	R3MeshVertex *closest_vertex = NULL;
   	RNLength closest_distance_to_vertex = FLT_MAX;
   	if (closest.type == R3_MESH_VERTEX_TYPE) {
    	 // Closest point was on vertex
     	closest_vertex = closest.vertex;
     	closest_distance_to_vertex = closest.t;
   	}
   	else if (closest.type == R3_MESH_EDGE_TYPE) {
     	// Closest point was on edge
     	R3MeshVertex *vertex0 = mesh->VertexOnEdge(closest.edge, 0);
     	RNLength d0 = R3Distance(pos, mesh->VertexPosition(vertex0));
     	if (d0 < closest_distance_to_vertex) {
       		closest_vertex = vertex0;
       		closest_distance_to_vertex = d0;
     	}
     	R3MeshVertex *vertex1 = mesh->VertexOnEdge(closest.edge,1);
     	RNLength d1 = R3Distance(pos, mesh->VertexPosition(vertex1));
     	if (d1 < closest_distance_to_vertex) {
       		closest_vertex = vertex1;
       		closest_distance_to_vertex = d1;
     	}
   	}
   	else if (closest.type == R3_MESH_FACE_TYPE) {
     	// Closest point was in middle of face
     	R3MeshVertex *vertex0 = mesh->VertexOnFace(closest.face, 0);
     	RNLength d0 = R3Distance(pos, mesh->VertexPosition(vertex0));
     	if (d0 < closest_distance_to_vertex) {
       		closest_vertex = vertex0;
       		closest_distance_to_vertex = d0;
     	}
     	R3MeshVertex *vertex1 = mesh->VertexOnFace(closest.face,1);
     	RNLength d1 = R3Distance(pos, mesh->VertexPosition(vertex1));
     	if (d1 < closest_distance_to_vertex) {
       		closest_vertex = vertex1;
       		closest_distance_to_vertex = d1;
     	}
     	R3MeshVertex *vertex2 = mesh->VertexOnFace(closest.face,2);
     	RNLength d2 = R3Distance(pos, mesh->VertexPosition(vertex2));
     	if (d2 < closest_distance_to_vertex) {
       		closest_vertex = vertex2;
       		closest_distance_to_vertex = d2;
     	}
   	}

	return closest_vertex;
}

