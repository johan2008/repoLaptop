#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Program: alignedmesh2halveskdtree
// Input 1: mesh
// Output:  dist_x dist_y dist_z
// Use kd-tree to accelerate the nearest neighbor search
// ////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name = NULL;
char * g_output_property_name[3] = {NULL, NULL, NULL};
char * g_output_error_name[3] = {NULL, NULL, NULL};
char * g_output_map_name[3] = {NULL, NULL, NULL};
char * g_output_avg_error_name[3] = {NULL, NULL, NULL};
char * g_output_mesh_name[3] = {NULL, NULL, NULL};

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * g_mesh = NULL;
R3MeshProperty* g_property[3] = {NULL, NULL, NULL};
R3MeshProperty* g_error[3] = {NULL, NULL, NULL};
struct DenseVertexCorrespondence {
  	R3Mesh *mesh[2];
  	R3MeshVertex **vertices; 
  	DenseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1) { 
    	mesh[0] = mesh0; 
    	mesh[1] = mesh1; 
    	vertices = new R3MeshVertex * [ mesh0->NVertices() ];
    	for (int i = 0; i < mesh0->NVertices(); i++) vertices[i] = NULL;
  	};
};
DenseVertexCorrespondence* g_map[3] = {NULL, NULL};
R3MeshSearchTree *g_search_tree = NULL;
////////////////////////////////////////////////////////////////////////
// Flags
////////////////////////////////////////////////////////////////////////

int print_verbose = 0;
int print_debug = 0;

////////////////////////////////////////////////////////////////////////
// Input functions
////////////////////////////////////////////////////////////////////////

// Mesh
static R3Mesh *
ReadMesh(char *filename)
{
  	// Start statistics
  	RNTime start_time;
  	start_time.Read();

  	// Allocate mesh
  	R3Mesh *mesh = new R3Mesh();
  	assert(mesh);

  	// Read mesh from file
  	if (!mesh->ReadFile(filename)) {
    	fprintf(stderr, "Unable to read mesh from %s\n", filename);
    	return NULL;
  	}

  	// Print statistics
  	if (print_verbose) {
    	printf("Read mesh from %s ...\n", filename);
    	printf("  Time = %.2f seconds\n", start_time.Elapsed());
    	printf("  # Faces = %d\n", mesh->NFaces());
    	printf("  # Edges = %d\n", mesh->NEdges());
    	printf("  # Vertices = %d\n", mesh->NVertices());
    	fflush(stdout);
  	}

  	// Return success
  	return mesh;
}

// Dense Correspondence
static int
WriteDenseVertexCorrespondence(DenseVertexCorrespondence *map, char *filename)
{
  	// Start statistics
  	RNTime start_time;
  	start_time.Read();

  	// Open file
  	FILE *fp = fopen(filename, "w");
  	if (!fp) {
    	fprintf(stderr, "Unable to open map file %s\n", filename);
    	return 0;
	}

	// Write correspondences
	for (int i = 0; i < map->mesh[0]->NVertices(); i++) {
		if (map->vertices[i] == NULL) fprintf(fp, "-1\n");
    	else fprintf(fp, "%d\n", map->mesh[1]->VertexID(map->vertices[i]));
  	}

  	// Close file
  	fclose(fp);

  	// Print statistics
  	if (print_verbose) {
    	printf("Wrote dense correspondences to %s ...\n", filename);
    	printf("  Time = %.2f seconds\n", start_time.Elapsed());
    	printf("  # Correspondences = %d\n", map->mesh[0]->NVertices());
    	fflush(stdout);
  	}

  	// Return success
  	return 1;
}

////////////////////////////////////////////////////////////////////////
// Argument Parser
////////////////////////////////////////////////////////////////////////

static int
ParseArgs(int argc, char ** argv)
{
	// Parse arguments
    argc--; argv++;
    while (argc > 0) {
        if ((*argv)[0] == '-') {
            if (!strcmp(*argv, "-v")) print_verbose = 1;
        	else if (!strcmp(*argv, "-debug")) print_debug = 1;
			else if (!strcmp(*argv, "-maps")) 
			{
				for (int i=0; i<3; i++)
				{
					argc--; argv++;
					g_output_map_name[i] = *argv;
				}
			}
			else if (!strcmp(*argv, "-errors"))
			{
				for (int i=0; i<3; i++)
				{
					argc--; argv++;
					g_output_error_name[i] = *argv;
				}
			}
			else if (!strcmp(*argv, "-avg_errors"))
			{
				for (int i=0; i<3; i++)
				{
					argc--; argv++;
					g_output_avg_error_name[i] = *argv;
				}
			}
			else if (!strcmp(*argv, "-output_scores"))
			{
				for (int i=0; i<3; i++)
				{
					argc--; argv++;
					g_output_avg_error_name[i] = *argv;
				}
			}
			else if (!strcmp(*argv, "-output_meshes"))
			{
				for (int i=0; i<3; i++)
				{
					argc--; argv++;
					g_output_mesh_name[i] = *argv;
				}
			}
			else if (!strcmp(*argv, "-error_properties"))
			{
				for (int i=0; i<3; i++)
				{
					argc--; argv++;
					g_output_property_name[i] = *argv;
				}
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!g_input_mesh_name) g_input_mesh_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name) {
    	fprintf(stderr, "Usage: ca_extrinsic_sym_fast [-output_meshes m0 m1 m2 -output_scores s0 s1 s2]\n");
    	return 0;
  	}

	return 1;
}

static void
sort_(int order[], RNScalar val[])
{
//	cout<<val[0]<<" "<<val[1]<<" "<<val[2]<<endl;
	if (val[order[0]] > val[order[1]])
	{
		int temp = order[1];
		order[1] = order[0];
		order[0] = temp;
	}
	if (val[order[1]] > val[order[2]])
	{
		int temp = order[2];
		order[2] = order[1];
		order[1] = temp;
	}
	if (val[order[0]] > val[order[1]])
	{
		int temp = order[1];
		order[1] = order[0];
		order[0] = temp;
	}
//	cout<<order[0]<<" "<<order[1]<<" "<<order[2]<<endl;
}

/*static void*/
/*sort_(int order[], RNScalar val[])*/
/*{*/
	/*if (val[0]<val[1]&&val[0]<val[2])*/
	/*{*/
		/*order[0] = 0;*/
		/*order[1] = 1;*/
		/*order[2] = 2;*/
	/*}*/
	/*else if (val[1]<val[2])*/
	/*{*/
		/*order[0] = 1;*/
		/*order[1] = 2;*/
		/*order[2] = 0;*/
	/*}*/
	/*else*/
	/*{*/
		/*order[0] = 2;*/
		/*order[1] = 0;*/
		/*order[2] = 1;*/
	/*}*/
/*}*/

static void
sort_(int order[], int x)
{
	if (x == 0)
	{
		order[0] = 0;
		order[1] = 1;
		order[2] = 2;
	}
	else if (x == 1)
	{
		order[0] = 1;
		order[1] = 2;
		order[2] = 0;
	}
	else
	{
		order[0] = 2;
		order[1] = 0;
		order[2] = 1;
	}
}

/*static void*/
/*sort_(int order[], RNScalar val[])*/
/*{*/
	/*if (val[0]<val[1]&&val[0]<val[2])*/
	/*{*/
		/*order[0] = 0;*/
		/*if (val[1]<val[2])*/
		/*{*/
			/*order[1] = 1;*/
			/*order[2] = 2;*/
		/*}*/
		/*else*/
		/*{*/
			/*order[2] = 1;*/
			/*order[1] = 2;*/
		/*}*/
	/*}*/
	/*else if (val[1]<val[2])*/
	/*{*/
		/*order[0] = 1;*/
		/*if (val[0]<val[2])*/
		/*{*/
			/*order[1] = 0;*/
			/*order[2] = 2;*/
		/*}*/
		/*else*/
		/*{*/
			/*order[1] = 2;*/
			/*order[2] = 0;*/
		/*}*/
	/*}*/
	/*else*/
	/*{*/
		/*order[0] = 2;*/
		/*if (val[0]<val[1])*/
		/*{*/
			/*order[1] = 0;*/
			/*order[2] = 1;*/
		/*}*/
		/*else*/
		/*{*/
			/*order[1] = 1;*/
			/*order[2] = 0;*/
		/*}*/
	/*}*/
/*}*/

/*static R3MeshVertex* */
//NearestNeighbor(R3Mesh* mesh, R3Point pos, RNScalar* ret_error, RNScalar dx, RNScalar dy, RNScalar dz)
//{
	//RNScalar error = FLT_MAX;
	//R3MeshVertex* best_v = NULL;
	//for (int i=0; i<mesh->NVertices(); i++)
	//{
		//R3MeshVertex* vertex = mesh->Vertex(i);
		//R3Point _p = mesh->VertexPosition(vertex);
		//if (dx < 0 && _p.X() > 0 || dx > 0 && _p.X() < 0) continue;
		//if (dy < 0 && _p.Y() > 0 || dy > 0 && _p.Y() < 0) continue;
		//if (dz < 0 && _p.Z() > 0 || dz > 0 && _p.Z() < 0) continue;
		//R3Vector dist = _p - pos;
		//if (dist.Length() < error)
		//{
			//error = dist.Length();
			//best_v = vertex;
		//}
	//}
	//if (ret_error)
	//{
		//*ret_error = error;
	//}
	//return best_v;
/*}*/
static R3MeshVertex* 
NearestNeighbor(R3Mesh* mesh, R3Point pos, RNScalar* ret_error, RNScalar dx, RNScalar dy, RNScalar dz)
{
	R3MeshIntersection closest;
	g_search_tree->FindClosest(pos, closest, 0, 1.0);
	if (closest.type == R3_MESH_NULL_TYPE)
	{
		*ret_error = 0.2;
		return NULL;
	}

	RNLength closest_distance_to_surface = closest.t;
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

	*ret_error = closest_distance_to_vertex;
	return closest_vertex;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	g_mesh = ReadMesh(g_input_mesh_name);

	// Search tree
	g_search_tree = new R3MeshSearchTree(g_mesh);
	if (!g_search_tree)
	{
		fprintf(stderr, "Unable to create search tree!\n");
		return 0;
	}
	for (int i=0; i<3; i++)
	{
		g_property[i] = new R3MeshProperty(g_mesh);
		g_map[i] = new DenseVertexCorrespondence(g_mesh, g_mesh);
		g_error[i] = new R3MeshProperty(g_mesh);
	}
	RNTime start_time;
	start_time.Read();
	for (int i=0; i<g_mesh->NVertices(); i++)
	{
		R3MeshVertex* v = g_mesh->Vertex(i);
		R3Point pos = g_mesh->VertexPosition(v);
		g_property[0]->SetVertexValue(v, fabs(pos.X()));
		g_property[1]->SetVertexValue(v, fabs(pos.Y()));
		g_property[2]->SetVertexValue(v, fabs(pos.Z()));
		
		R3MeshVertex* _v[3];
		RNScalar error[3];
	   	_v[0] = NearestNeighbor(g_mesh, R3Point(-pos.X(), pos.Y(), pos.Z()), &error[0], -pos.X(), 0.0, 0.0);
		_v[1] = NearestNeighbor(g_mesh, R3Point(pos.X(), -pos.Y(), pos.Z()), &error[1], 0.0, -pos.Y(), 0.0);
		_v[2] = NearestNeighbor(g_mesh, R3Point(pos.X(), pos.Y(), -pos.Z()), &error[2], 0.0, 0.0, -pos.Z());

		for (int j=0; j<3; j++)
		{
			g_map[j]->vertices[i] = _v[j];
			g_error[j]->SetVertexValue(v, error[j]);
		}
	}
	cout<<"***** Time for finding Nearest neighbors : "<<start_time.Elapsed()<<" seconds"<<endl;
	
	// Sort 3 directions by error
	int order[3] = {0, 1, 2};
	RNScalar avg_errors[3];
	for (int i=0; i<3; i++)
		avg_errors[i] = g_error[i]->Mean();

	sort_(order, avg_errors);
//	if (g_error[0]->Mean() < g_error[1]->Mean() && g_error[0]->Mean() < g_error[2]->Mean())
//	{
//		if (g_error[1]->Mean() < g_error[2]->Mean())
//		{
//
//		}
//	}
	for (int p=0; p<3; p++)
	{
		int i = order[p];
		if (g_output_property_name[p])
			g_property[i]->Write(g_output_property_name[p]);
		if (g_output_map_name[p])
			WriteDenseVertexCorrespondence(g_map[i], g_output_map_name[p]);
		if (g_output_error_name[p])
			g_error[i]->Write(g_output_error_name[p]);
		if (g_output_avg_error_name[p])
		{
			ofstream fout(g_output_avg_error_name[p]);
			fout<<g_error[i]->Mean()<<endl;
			fout.close();
		}
		cout<<"Error "<<p<<" : "<<g_error[i]->Mean()<<endl;
	}
	int _order[3];
	_order[0] = order[0]; _order[1] = order[1]; _order[2] = order[2];
	// re-arrange the coordinates of vertices
	vector<R3Point> old_positions;
	for (int i=0; i<g_mesh->NVertices(); i++)
		old_positions.push_back(g_mesh->VertexPosition(g_mesh->Vertex(i)));

	if (g_output_mesh_name[0] && g_output_mesh_name[1] && g_output_mesh_name[2])
	{
		for (int kk=0; kk<3; kk++)
		{
			int considered = _order[kk];
			sort_(order, considered);
			vector<R3Point> new_positions;
			for (int i=0; i<g_mesh->NVertices(); i++)
			{
				new_positions.push_back(R3Point(0,0,0));
			}
			for (int i=0; i<g_mesh->NVertices(); i++)
			{
				for (int p=0; p<3; p++)
				{
					new_positions[i][p] = old_positions[i][order[p]];
				}
			}
			for (int i=0; i<g_mesh->NVertices(); i++)
			{
				g_mesh->SetVertexPosition(g_mesh->Vertex(i), new_positions[i]);
			}
			g_mesh->WriteFile(g_output_mesh_name[kk]);
		}
	}
	return 0;
}
