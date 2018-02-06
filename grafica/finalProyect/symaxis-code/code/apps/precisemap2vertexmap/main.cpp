#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////
char * input_mesh_name = NULL;
char * input_precise_map_name = NULL;
char * output_vertex_map_name = NULL;
////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * Mesh = NULL;

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
			if (!input_mesh_name) input_mesh_name = *argv;
			else if (!input_precise_map_name) input_precise_map_name = *argv;
			else if (!output_vertex_map_name) output_vertex_map_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!input_mesh_name || !input_precise_map_name || !output_vertex_map_name) {
    	fprintf(stderr, "Usage: PreciseMap2VertexMap mesh precise_map vertex_map [-v, -debug]\n");
    	return 0;
  	}

	return 1;
}

int main(int argc, char ** argv)
{

	if (!ParseArgs(argc, argv)) exit(1);

	Mesh = ReadMesh(input_mesh_name);
	
	if (!Mesh)
	{
		cout<<"Cannot load mesh!"<<endl;
		exit(-1);
	}

	RNArray<R3MeshVertex*> vertex_map;
	ifstream fin(input_precise_map_name);
	if (!fin)
	{
		cout<<"Cannot open precise map "<<input_precise_map_name<<endl;
		exit(-1);
	}
	for (int i=0; i<Mesh->NVertices(); i++)
	{
		int idx;
		RNScalar weights[3];
		fin >> idx >> weights[0] >> weights[1] >> weights[2];
		R3MeshFace* face = Mesh->Face(idx);
		R3MeshVertex* v[3];
		v[0] = Mesh->VertexOnFace(face, 0);
		v[1] = Mesh->VertexOnFace(face, 1);
		v[2] = Mesh->VertexOnFace(face, 2);
		if (weights[0] > weights[1] && weights[0] > weights[2])
			vertex_map.Insert(v[0]);
		else if (weights[1] > weights[2])
			vertex_map.Insert(v[1]);
		else
			vertex_map.Insert(v[2]);
	}
	fin.close();

	ofstream fout(output_vertex_map_name);
	if (!fout)
	{
		cout<<"Cannot open vertex map "<<output_vertex_map_name<<endl;
		exit(-1);
	}
	for (int i=0; i<vertex_map.NEntries(); i++)
	{
		R3MeshVertex* v = vertex_map.Kth(i);
		fout<<Mesh->VertexID(v)<<endl;
	}
	fout.close();
	return 0;
}
