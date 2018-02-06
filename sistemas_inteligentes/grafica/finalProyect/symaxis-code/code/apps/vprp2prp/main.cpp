#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Program: v_feat2prp
// Input 1: mesh
// Input 2: v_features
// Output: *.val
// A program that converts a vova's style feature to *.val
// Format: featurename, featurename2, #values, values
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name = NULL;
char * g_input_v_feature_name = NULL;
char * g_output_prp_name = NULL;

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * g_mesh = NULL;
R3MeshProperty* g_property = NULL;

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
			else if (!g_input_v_feature_name) g_input_v_feature_name = *argv;
			else if (!g_output_prp_name) g_output_prp_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name || !g_input_v_feature_name || !g_output_prp_name) {
    	fprintf(stderr, "Usage: v_feat2prp mesh v_feat prp\n");
    	return 0;
  	}

	return 1;
}

static R3MeshProperty*
ReadVFeature(R3Mesh* mesh, const char* v_feat)
{
	ifstream fin(v_feat);
	if (!fin)
	{
		cout<<"Unable to open file "<<v_feat<<endl;
		return NULL;
	}
	char feat_name[200];
	fin >> feat_name;
	R3MeshProperty* property = new R3MeshProperty(mesh, feat_name);
	fin >> feat_name; // Skip one word
	int n_vals;
	fin >> n_vals;
	if (n_vals != mesh->NVertices())
	{
		cout<<"Num of values does not match with the number of vertices!"<<endl;
		fin.close();
		delete property;
		return NULL;
	}
	double temp;
	for (int i=0; i<n_vals; i++)
	{
		fin >> temp;
		property->SetVertexValue(i, temp);
	}
	return property;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);

	g_mesh = ReadMesh(g_input_mesh_name);

	g_property = ReadVFeature(g_mesh, g_input_v_feature_name);

	g_property->Write(g_output_prp_name);
	return 0;
}
