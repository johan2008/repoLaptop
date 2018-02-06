#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
#include <string>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * input_mesh_name = NULL;
char * input_map_name = NULL;
char * output_sym_name = NULL;
////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * Mesh = NULL;
R3MeshProperty* Property = NULL;

RNLength max_distance = 0.3;

////////////////////////////////////////////////////////////////////////
// Flags
////////////////////////////////////////////////////////////////////////

int print_verbose = 0;
int print_debug = 0;
int print_time = 0;
////////////////////////////////////////////////////////////////////////
// Helper
////////////////////////////////////////////////////////////////////////
static void
PrintProgressBar(int percent)
{
    string bar;
    for (int i=0; i<50; i++)
    {
        if (i < (percent/2))
        {
            bar.replace(i, 1, "=");
        }
        else if ( i == (percent/2))
        {
            bar.replace(i, 1, ">");
        }
        else
        {
            bar.replace(i, 1, " ");
        }
    }
    printf("\r        [%s] ",bar.c_str());
    printf("%3d%%", percent);
    fflush(stdout);
}

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
			else if (!strcmp(*argv, "-time"))
			{
				print_time = 1;
			}
			else if (!strcmp(*argv, "-max_dist"))
			{
				argc--; argv++;
				max_distance = atof(*argv);
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!input_mesh_name) input_mesh_name = *argv;
			else if (!input_map_name) input_map_name = *argv;
			else if (!output_sym_name) output_sym_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!input_mesh_name || !input_map_name || !output_sym_name) {
    	fprintf(stderr, "Usage: map2sym mesh input_map output_sym\n");
    	return 0;
  	}

	return 1;
}

static RNArray<R3MeshVertex*>*
ReadMap(R3Mesh* mesh, char* filename)
{
	ifstream fin(filename);
	if (!fin)
	{
		cout<<"Cannot open map file "<<filename<<endl;
		exit(-1);
	}
	RNArray<R3MeshVertex*>* array = new RNArray<R3MeshVertex*>;
	int idx;
	while (fin >> idx)
	{
		array->Insert(mesh->Vertex(idx));
	}
	fin.close();
	return array;
}

int main(int argc, char ** argv)
{

	if (!ParseArgs(argc, argv)) exit(1);
	Mesh = ReadMesh(input_mesh_name);

	RNArray<R3MeshVertex*>* dense_map = ReadMap(Mesh, input_map_name);

    cout<<"Computing geodesic distances between correspondences from mapping ..."<<endl;
	RNTime start_time;
	start_time.Read();
	Property = new R3MeshProperty(Mesh);
	for (int i=0; i<Mesh->NVertices(); i++)
	{
		R3MeshVertex* src = Mesh->Vertex(i);
		R3MeshVertex* dst = dense_map->Kth(i);

	//	RNScalar dist = Mesh->DijkstraDistance(src, dst);
		RNScalar* dist = Mesh->DijkstraDistances(src, max_distance);
		RNLength _d = dist[Mesh->VertexID(dst)];
		if (_d>max_distance)
			_d = max_distance;
		Property->SetVertexValue(src, _d);
        delete[] dist;
        PrintProgressBar((i+1)*100/Mesh->NVertices());
	}
	cout<<"******Time: "<<start_time.Elapsed()<<" seconds\n";

	Property->Write(output_sym_name);
	return 0;
}
