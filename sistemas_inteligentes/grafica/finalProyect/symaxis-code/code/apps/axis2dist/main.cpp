#include "R3Shapes/R3Shapes.h"
#include "RNExt/RNExt.h"

////////////////////////////////////////////////////////////////////////
// Program     : axis2dist
// Usage	   : axis2dist mesh axis prp(dist) [-euclidean -geodesic]
// Description : compute geodesic/euclidean distance from the axis
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name = NULL;
char * g_input_axis_name = NULL;
char * g_output_dist_name = NULL;

int g_is_geodesic = 1;

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
			else if (!strcmp(*argv, "-geodesic"))
			{
				g_is_geodesic = 1;
			}
			else if (!strcmp(*argv, "-euclidean"))
			{
				g_is_geodesic = 0;
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!g_input_mesh_name) g_input_mesh_name = *argv;
			else if (!g_input_axis_name) g_input_axis_name = *argv;
			else if (!g_output_dist_name) g_output_dist_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name || !g_input_axis_name || !g_output_dist_name) {
    	fprintf(stderr, "Usage: axis2dist mesh axis prp(dist)\n");
    	return 0;
  	}

	return 1;
}

static RNLength*
MinimumEuclideanDistances(R3Mesh* mesh, vector<R3Point>* points)
{
	RNLength* dists = new RNLength[mesh->NVertices()];
	for (int i=0; i<mesh->NVertices(); i++)
	{
		R3MeshVertex* vertex = mesh->Vertex(i);
		R3Point pos = mesh->VertexPosition(vertex);
		RNLength min_dist = FLT_MAX;
		for (int j=0; j<int(points->size()); j++)
		{
			R3Point _p = (*points)[j];
			RNLength d = R3Distance(pos, _p);
			if (d < min_dist)
				min_dist = d;
		}
		dists[i] = min_dist;
	}
	return dists;
}

static vector<R3Point>*
ReadPositions(R3Mesh* mesh, const char* filename)
{
	ifstream fin(filename);
	if (!fin)
	{
		fprintf(stderr, "Unable to open %s!\n", filename);
		return NULL;
	}
	vector<R3Point>* array = new vector<R3Point>;
	RNScalar t[3];
	while (fin >> t[0] >> t[1] >> t[2])
	{
		R3Point p(t[0], t[1], t[2]);
		array->push_back(p);
	}
	fin.close();
	return array;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	g_mesh = ReadMesh(g_input_mesh_name);

	if (g_is_geodesic)
	{
		RNArray<R3MeshVertex*>* axis = ReadPoints(g_mesh, g_input_axis_name);

		RNLength* distances = g_mesh->DijkstraDistances(*axis);

		g_property = new R3MeshProperty(g_mesh, "djikstra", distances);
		g_property->Write(g_output_dist_name);
	}
	else
	{
		vector<R3Point>* axis = ReadPositions(g_mesh, g_input_axis_name);
		RNLength* distances = MinimumEuclideanDistances(g_mesh, axis);

		g_property = new R3MeshProperty(g_mesh, "euclidean", distances);
		g_property->Write(g_output_dist_name);
	}
	return 0;
}
