#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Program: axislength
// Input 1: mesh
// Input 2: axis
// Output: length.txt
// Just a template provided for desiging other programs
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name = NULL;
char * g_input_axis_name = NULL;
char * g_output_length_name = NULL;

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * g_mesh = NULL;
RNArray<R3MeshVertex*>* g_axis = NULL;
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
/*  	if (print_verbose) {*/
        //printf("Read mesh from %s ...\n", filename);
        //printf("  Time = %.2f seconds\n", start_time.Elapsed());
        //printf("  # Faces = %d\n", mesh->NFaces());
        //printf("  # Edges = %d\n", mesh->NEdges());
        //printf("  # Vertices = %d\n", mesh->NVertices());
        //fflush(stdout);
      //}

  	// Return success
  	return mesh;
}

// Points
static RNArray<R3MeshVertex*>*
ReadPoints(R3Mesh* mesh, const char* filename)
{
	RNArray<R3MeshVertex*>* axis = new RNArray<R3MeshVertex*>;
	ifstream fin(filename);
	if (!fin)
	{
		cout<<"Cannot open point file "<<filename<<endl;
		exit(-1);
	}
	int temp;
	while (fin>>temp)
	{
		axis->Insert(mesh->Vertex(temp));
	}
	fin.close();
	return axis;
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
			else if (!g_input_axis_name) g_input_axis_name = *argv;
			else if (!g_output_length_name) g_output_length_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name || !g_input_axis_name || !g_output_length_name) {
    	fprintf(stderr, "Usage: axislength mesh axis length\n");
    	return 0;
  	}

	return 1;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	g_mesh = ReadMesh(g_input_mesh_name);
	g_axis = ReadPoints(g_mesh, g_input_axis_name);

	RNLength l = 0;
	for (int i=0; i<g_axis->NEntries(); i++)
	{
		R3MeshVertex* vertex[2];
	    vertex[0] = g_axis->Kth(i);
		vertex[1] = g_axis->Kth((i+1)%g_axis->NEntries());
		R3MeshEdge* edge = g_mesh->EdgeBetweenVertices(vertex[0], vertex[1]);
		if (!edge)
		{
			cout<<"Axis is disconnected!"<<endl;
			exit(-1);
		}
		l += g_mesh->EdgeLength(edge);
	}

	ofstream fout(g_output_length_name);
	if (!fout)
	{
		cout<<"Unable to open file "<<g_output_length_name<<endl;
		exit(-1);
	}
	fout<<l<<endl;
	fout.close();
	return 0;
}
