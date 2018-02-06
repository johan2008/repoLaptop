#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include "thinning.h"
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * input_mesh_name = NULL;
char * input_sym_name = NULL;
char * input_binary_name = NULL;
char * output_axis_name = NULL;
////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * Mesh = NULL;
R3MeshProperty* Property = NULL;
R3MeshProperty* extern_binary_property = NULL;
int MinLength = 10;
////////////////////////////////////////////////////////////////////////
// Flags
////////////////////////////////////////////////////////////////////////

int print_verbose = 0;
int print_debug = 0;
int print_time = 0;

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
			else if (!strcmp(*argv, "-time"))
			{
				print_time = 1;
			}
			else if (!strcmp(*argv, "-extern_binary"))
			{
				argv++; argc--;
				input_binary_name = *argv;
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!input_mesh_name) input_mesh_name = *argv;
			else if (!input_sym_name) input_sym_name = *argv;
			else if (!output_axis_name) output_axis_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!input_mesh_name || !input_sym_name || !output_axis_name) {
    	fprintf(stderr, "Usage: sym2axis mesh sym axis [ -v, -debug]\n");
    	return 0;
  	}

	return 1;
}

int main(int argc, char ** argv)
{

	if (!ParseArgs(argc, argv)) exit(1);

	Mesh = ReadMesh(input_mesh_name);

	Property = new R3MeshProperty(Mesh);
	Property->Read(input_sym_name);

	if (input_binary_name)
	{
		extern_binary_property = new R3MeshProperty(Mesh);
		extern_binary_property->Read(input_binary_name);
		// make center to be 1
		for (int i=0; i<Mesh->NVertices(); i++)
		{
			RNScalar val = extern_binary_property->VertexValue(i);
			val = 1.0 - val;
			extern_binary_property->SetVertexValue(i, val);
		}
	}
	RNTime start_time;
	start_time.Read();
	RNArray<R3MeshVertex*>* axis = Thinning(Property);
	if (print_time)
	{
		cout<<"Time: "<<start_time.Elapsed()<<endl;
	}

	WritePoints(Mesh, axis, output_axis_name);
	return 0;
}
