#include "R3Shapes/R3Shapes.h"
#include "RNExt/RNExt.h"
#include <string>
#include <sstream>

////////////////////////////////////////////////////////////////////////
// Program     : map2error
// Usage       : map2error mesh0 mesh1 map pts0 pts1 error
// Description : compute the euclidean error
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name[2] = {NULL, NULL};
char * g_input_map_name = NULL;
char * g_input_pts_name[2] = {NULL, NULL};
char * g_output_error_name = NULL;

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * g_mesh[2] = {NULL, NULL};
DenseVertexCorrespondence* g_map = NULL;
RNArray<R3MeshVertex*>* g_pts[2] = {NULL, NULL};
int g_geodesic = 0;
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
				g_geodesic = 1;
			}
			else if (!strcmp(*argv, "-euclidean"))
			{
				g_geodesic = 0;
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!g_input_mesh_name[0]) g_input_mesh_name[0] = *argv;
			else if (!g_input_mesh_name[1]) g_input_mesh_name[1] = *argv;
			else if (!g_input_map_name) g_input_map_name = *argv;
			else if (!g_input_pts_name[0]) g_input_pts_name[0] = *argv;
			else if (!g_input_pts_name[1]) g_input_pts_name[1] = *argv;
			else if (!g_output_error_name) g_output_error_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name[0] || !g_input_mesh_name[1] || !g_input_map_name ||
			!g_input_pts_name[0] || !g_input_pts_name[1] ||
			!g_output_error_name) {
    	fprintf(stderr, "Usage: map2error mesh0 mesh1 map pts0 pts1 error\n");
    	return 0;
  	}

	return 1;
}

static RNScalar
Average(vector<RNScalar>& array)
{
	RNScalar avg = 0;
	for (unsigned i=0; i<array.size(); i++)
	{
		avg += array[i];
	}
	avg /= int(array.size());
	return avg;
}

static RNArray<R3MeshVertex*>*
ReadPoints2(R3Mesh* mesh, const char* filename)
{
	RNArray<R3MeshVertex*>* vertices = new RNArray<R3MeshVertex*>;
	ifstream fin(filename);
	string temp;
	while (getline(fin, temp))
	{
		int idx;
		istringstream is(temp);
		is >> idx;
		vertices->Insert(mesh->Vertex(idx));	
//		cout<<idx<<endl;
	}
	fin.close();
	return vertices;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	g_mesh[0] = ReadMesh(g_input_mesh_name[0]);
	g_mesh[1] = ReadMesh(g_input_mesh_name[1]);
	g_map = ReadDenseVertexCorrespondence(g_mesh[0], g_mesh[1], g_input_map_name);
	g_pts[0] = ReadPoints2(g_mesh[0], g_input_pts_name[0]);
	g_pts[1] = ReadPoints2(g_mesh[1], g_input_pts_name[1]);

	RNLength mesh_radius = g_mesh[1]->AverageRadius();
	if (g_pts[0]->NEntries() != g_pts[1]->NEntries())
	{
		fprintf(stderr, "Number of points don't match!\n");
		exit(-1);
	}

	vector<RNScalar> errors;
	for (int i=0; i<g_pts[0]->NEntries(); i++)
	{
		int idx = g_mesh[0]->VertexID(g_pts[0]->Kth(i));
		R3MeshVertex* v = g_map->vertices[idx];
		if (g_geodesic)
		{
			errors.push_back(g_mesh[1]->DijkstraDistance(v, g_pts[1]->Kth(i)));
		}
		else
		{
			R3Point pos[2];
			pos[0] = g_mesh[1]->VertexPosition(v);
			pos[1] = g_mesh[1]->VertexPosition(g_pts[1]->Kth(i));
			errors.push_back(R3Distance(pos[0], pos[1])/mesh_radius);
		}
	}
	
	RNScalar avg = Average(errors);
	ofstream fout(g_output_error_name);
	// first line
	fout<<"Stats:"<<endl;
	fout<<"    Average ["<<avg<<"]"<<endl;
	fout<<"    Folder  ["<<0<<"]"<<endl;
	fout<<"    Guessed ["<<0.0<<"]"<<endl;
	fout<<"    ObjFn   ["<<0<<"]"<<endl;
	fout<<"PerSampleErrors [";
	for (unsigned i=0; i<errors.size()-1; i++)
	{
		fout<<errors[i]<<", ";
	}
	fout<<errors[errors.size()-1];
	fout<<"]"<<endl;
	fout.close();
	return 0;
}
