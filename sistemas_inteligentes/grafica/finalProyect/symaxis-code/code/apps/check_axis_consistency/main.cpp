#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Program: check_axis_consistency
// Input 1: mesh
// Input 2: tips
// Input 3: axis
// Input 4: symmetry.val
// Output: consistency.txt (contains 'yes' or 'no')
// Check whether the axis is consistent with the symmetry.val (only for tips!)
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name = NULL;
char * g_input_tips_name = NULL;
char * g_input_axis_name = NULL;
char * g_input_sym_name = NULL;
char * g_output_ans_name = NULL;

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * g_mesh = NULL;
RNArray<R3MeshVertex*>* g_tips = NULL;
RNArray<R3MeshVertex*>* g_axis = NULL;
R3MeshProperty* g_sym = NULL;

////////////////////////////////////////////////////////////////////////
// Flags
////////////////////////////////////////////////////////////////////////

int print_verbose = 0;
int print_debug = 0;

RNScalar g_sym_thres = 0.05;
RNScalar g_axis_thres = 0.25;
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

static R3MeshProperty *
ReadProperty(R3Mesh *mesh, const char *filename)
{
	  // Start statistics
  	RNTime start_time;
  	start_time.Read();

  	// Allocate properties
 	R3MeshProperty *property = new R3MeshProperty(mesh);
  	if (!property) {
    	fprintf(stderr, "Unable to allocate property for %s\n", filename);
    	return NULL;
  	}

  	// Read properties from file
  	if (!property->Read(filename)) {
    	delete property;
    	return NULL;
  	}

  	// Print statistics
  	if (print_verbose) {
    	printf("Read property from %s ...\n", filename);
    	printf("  Time = %.2f seconds\n", start_time.Elapsed());
    	printf("  # Vertices = %d\n", property->Mesh()->NVertices());
    	fflush(stdout);
  	}

  	// Return property set
  	return property;
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
			else if (!strcmp(*argv, "-thres_axis"))
			{
				argc--; argv++;
				g_axis_thres = atof(*argv);
			}
			else if (!strcmp(*argv, "-thres_sym"))
			{
				argc--; argv++;
				g_sym_thres = atof(*argv);
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!g_input_mesh_name) g_input_mesh_name = *argv;
			else if (!g_input_tips_name) g_input_tips_name = *argv;
			else if (!g_input_axis_name) g_input_axis_name = *argv;
			else if (!g_input_sym_name) g_input_sym_name = *argv;
			else if (!g_output_ans_name) g_output_ans_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name || !g_input_tips_name || !g_input_axis_name 
			|| !g_input_sym_name || !g_output_ans_name) {
    	fprintf(stderr, "Usage: check_axis_consistency mesh tips axis sym_val ans [-v, -debug, -vicinity]\n");
    	return 0;
  	}

	return 1;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	g_mesh = ReadMesh(g_input_mesh_name);
	g_tips = ReadPoints(g_mesh, g_input_tips_name);
	g_axis = ReadPoints(g_mesh, g_input_axis_name);
	g_sym = ReadProperty(g_mesh, g_input_sym_name);

	int off_axis_tips_counter = 0;
	RNLength * dists = g_mesh->DijkstraDistances(*g_axis);
	// Pick the tips near to the 0-region in sym_val, to see if it is close to axis
	for (int i=0; i<g_tips->NEntries(); i++)
	{
		R3MeshVertex* tip = g_tips->Kth(i);
		if (g_sym->VertexValue(tip) < 2*g_sym_thres)
		{
			if (print_verbose)
			{
				cout<<"Vertex "<<g_mesh->VertexID(tip)<<" [ Distance to axis: "<<dists[g_mesh->VertexID(tip)]
					<<" , Symmetry value: "<<g_sym->VertexValue(tip)<<" ]"<<endl;
			}
			// Considered ought to be close to axis
			if (dists[g_mesh->VertexID(tip)] > g_axis_thres) 
			{
				cout<<"		Inconsistent!"<<endl;
				off_axis_tips_counter ++;
			}
			else
			{
				cout<<"		Consistent!"<<endl;
			}
		}
	}
	ofstream fout(g_output_ans_name);
	if (off_axis_tips_counter == 0)
	{
		fout<<"yes"<<endl;
	}
	else
	{
		fout<<"no"<<endl;
	}
	fout.close();
	return 0;
}
