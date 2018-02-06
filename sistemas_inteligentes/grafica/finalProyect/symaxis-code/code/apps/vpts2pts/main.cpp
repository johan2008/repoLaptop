#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Program: v_pts2pts
// Input 1: mesh
// Input 2: v_pts_name
// Output: *.pid
// Convert vova-style points to *pid
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name = NULL;
char * g_input_v_pts_name = NULL;
char * g_output_pts_name = NULL;

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * g_mesh = NULL;
RNArray<R3MeshVertex*>* g_points = NULL;

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

// Points
static int
WritePoints(R3Mesh* mesh, RNArray<R3MeshVertex*> *points, const char* filename)
{
	ofstream fout(filename);
	if (!fout)
	{
		printf("Cannot open file for writing Points : %s\n", filename);
		return 0;
	}
	for (int i=0; i<points->NEntries(); i++)
		fout<<mesh->VertexID(points->Kth(i))<<endl;

	// Close file
	fout.close();
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
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!g_input_mesh_name) g_input_mesh_name = *argv;
			else if (!g_input_v_pts_name) g_input_v_pts_name = *argv;
			else if (!g_output_pts_name) g_output_pts_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name || !g_input_v_pts_name || !g_output_pts_name) {
    	fprintf(stderr, "Usage: v_pts2pts mesh v_pts pts [-v, -debug]\n");
    	return 0;
  	}

	return 1;
}

static RNArray<R3MeshVertex*> *
ReadVPoints(R3Mesh* mesh, const char* v_pts_name)
{
	ifstream fin(v_pts_name);
	if (!fin)
	{
		cout<<"Unable to open file "<<v_pts_name<<endl;
		return NULL;
	}
	RNArray<R3MeshVertex*>* points = new RNArray<R3MeshVertex*>;
	int num = 0;
	fin >> num;
	for (int i=0; i<num; i++)
	{
		char t_chars[200];
		int faceidx = 0;
		double a, b, c;
		fin >> t_chars >> faceidx >> a >> b >> c;
		R3MeshFace* face = mesh->Face(faceidx);
		R3MeshVertex* v[3];
		v[0] = mesh->VertexOnFace(face, 0);
		v[1] = mesh->VertexOnFace(face, 1);
		v[2] = mesh->VertexOnFace(face, 2);
		if (a > b && a > c)
			points->Insert(v[0]);
		else if ( b > c)
			points->Insert(v[1]);
		else
			points->Insert(v[2]);
	}
	fin.close();
	return points;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	g_mesh = ReadMesh(g_input_mesh_name);
	g_points = ReadVPoints(g_mesh, g_input_v_pts_name);
	WritePoints(g_mesh, g_points, g_output_pts_name);
	return 0;
}
