#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Program: AxisRefine
// Input 1: meshname
// Input 2: axisname
// Output: refined axisname
// Generate a connected axis, free of repetition
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name = NULL;
char * g_input_axis_name = NULL;
char * g_output_axis_name = NULL;

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
int flag_check_only = 0;

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
			else if (!strcmp(*argv, "-check_only"))
			{
				flag_check_only = 1;
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!g_input_mesh_name) g_input_mesh_name = *argv;
			else if (!g_input_axis_name) g_input_axis_name = *argv;
			else if (!g_output_axis_name) g_output_axis_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name || !g_input_axis_name || !g_output_axis_name) {
    	fprintf(stderr, "Usage: AxisRefine meshname input_axis output_axis [-v, -debug, -check_only]\n");
    	return 0;
  	}

	return 1;
}

// Check the connectivity and repetition
static int
CheckAxis(R3Mesh* mesh, RNArray<R3MeshVertex*>* axis)
{
	// Connectivity
	for (int i=0; i<axis->NEntries(); i++)
	{
		R3MeshVertex* v[2];
		v[0] = axis->Kth(i);
		v[1] = axis->Kth((i+1)%axis->NEntries());
		R3MeshEdge* edge = mesh->EdgeBetweenVertices(v[0], v[1]);
		if (!edge)
		{
			cout<<"Axis is disconnected between vertex "<<mesh->VertexID(v[0])<<" and vertex "
				<<mesh->VertexID(v[1])<<endl;
			return 0;
		}
	}
	cout<<"Pass the test of connectivity."<<endl;
	// Free of repetition
	R3mesh_mark ++;
	for (int i=0; i<axis->NEntries(); i++)
	{
		R3MeshVertex* v = axis->Kth(i);
		int mark = mesh->VertexMark(v);
		if (mark == R3mesh_mark)
		{
			cout<<"Vertex "<<mesh->VertexID(v)<<" is repeated."<<endl;
			return 0;
		}
		mesh->SetVertexMark(v, R3mesh_mark);
	}
	cout<<"Pass the test of repetition."<<endl;
	return 1;
}

// Refinement of axis
static int
RefineAxis(R3Mesh* mesh, RNArray<R3MeshVertex*>* axis)
{
	// Connecting
	for (int i=0; i<axis->NEntries(); i++)
	{
		R3MeshVertex* v[2];
		v[0] = axis->Kth(i);
		v[1] = axis->Kth((i+1)%axis->NEntries());
		R3MeshEdge* edge = mesh->EdgeBetweenVertices(v[0], v[1]);
		if (!edge)
		{
			RNArray<R3MeshEdge*> edges;
			mesh->DijkstraDistance(v[1], v[0], &edges);
			R3MeshVertex* prev = v[0];
			for (int j=0; j<edges.NEntries()-1; j++)
			{
				prev = mesh->VertexAcrossEdge(edges.Kth(j), prev);
				axis->InsertKth(prev, (i+1)%axis->NEntries());
				i++;
			}
		}
	}
	// Remove repetition
	int * pos_axis = new int[mesh->NVertices()];
	while (1)
	{
		for (int i=0; i<mesh->NVertices(); i++)
		{
			pos_axis[i] = -1;
		}
		int start = -1, end = -1;
		for (int i=0; i<axis->NEntries(); i++)
		{
			int idx = mesh->VertexID(axis->Kth(i));
			if (pos_axis[idx] == -1)
			{
				pos_axis[idx] = i;
			}
			else
			{
				start = pos_axis[idx];
				end = i;
			}
		}
		if (start >= 0 && end >= 0)
		{
			if (end - start < axis->NEntries() / 2)
			{
				// Remove start to end-1
				for (int i=start; i<end; i++)
				{
					axis->RemoveKth(start);
				}
			}
			else
			{
				// Remove end to actual end
				while (axis->NEntries() > end)
				{
					axis->RemoveKth(axis->NEntries()-1);
				}
				// Remove 0 to start-1
				for (int i=0; i<start; i++)
				{
					axis->RemoveKth(0);
				}
			}
		}
		else
		{
			break;
		}
	}
	delete[] pos_axis;
	return 1;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	g_mesh = ReadMesh(g_input_mesh_name);
	g_axis = ReadPoints(g_mesh, g_input_axis_name);

	if (!flag_check_only)
	{
		cout<<"#Vertices before refinement : "<<g_axis->NEntries()<<endl;
		int suc_refine = RefineAxis(g_mesh, g_axis);
		if (!suc_refine)
		{
			cout<<"Error during RefineAxis!"<<endl;
			exit(-1);
		}
		cout<<"#Vertices after refinement : "<<g_axis->NEntries()<<endl;
	}
	else
	{
		cout<<"In the\'check-only\' mode"<<endl;
	}
	// Check the connetivity and repetition
	int suc_check = CheckAxis(g_mesh, g_axis);
	if (!suc_check)
	{
		cout<<"Fail to pass the check!"<<endl;
		exit(-1);
	}
	if (!flag_check_only)
	{
		WritePoints(g_mesh, g_axis, g_output_axis_name);
	}
	return 0;
}
