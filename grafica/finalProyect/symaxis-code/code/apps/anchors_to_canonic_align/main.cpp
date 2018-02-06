#include "R3Shapes/R3Shapes.h"
#include "RNExt/RNExt.h"

////////////////////////////////////////////////////////////////////////
// Program     : samples_anchors_to_canonic_align
// Usage       : samples_anchors_to_canonic_align mesh0 mesh1 samples1 samples2 anchors1 anchors2 canonic_alignment
// Description : given a sequence of samples on each mesh (positions, in the same orientations), and anchors (in correspondence), output a interpolated alignment
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name[2] = {NULL, NULL};
char * g_input_axis_name[2] = {NULL, NULL};
char * g_input_anchors_name[2] = {NULL, NULL};
char * g_output_alignment_name = NULL;

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * g_mesh[2] = {NULL, NULL};
//vector<R3Point>* g_samples[2] = {NULL, NULL};
RNArray<R3MeshVertex*>* g_axis[2] = {NULL, NULL};
RNArray<R3MeshVertex*>* g_anchors[2] = {NULL, NULL};
//vector<int> g_align[2];
SparseVertexCorrespondence* g_cor = NULL;

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
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!g_input_mesh_name[0]) g_input_mesh_name[0] = *argv;
			else if (!g_input_mesh_name[1]) g_input_mesh_name[1] = *argv;
			else if (!g_input_axis_name[0]) g_input_axis_name[0] = *argv;
			else if (!g_input_axis_name[1]) g_input_axis_name[1] = *argv;
			else if (!g_input_anchors_name[0]) g_input_anchors_name[0] = *argv;
			else if (!g_input_anchors_name[1]) g_input_anchors_name[1] = *argv;
			else if (!g_output_alignment_name) g_output_alignment_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name[0] || !g_input_mesh_name[1] ||
			!g_input_axis_name[0] || !g_input_axis_name[1] ||
			!g_input_anchors_name[0] || !g_input_anchors_name[1] ||
			!g_output_alignment_name) {
    	fprintf(stderr, "Usage: samples_anchors_to_canonic_align mesh0 mesh1 samples1 samples2 anchors1 anchors2 canonic_alignment\n");
    	return 0;
  	}

	return 1;
}

static vector<R3Point>*
ReadSamples(const char* filename)
{
	ifstream fin(filename);
	if (!fin)
	{
		printf("Unable to open file %s!\n", filename);
		return NULL;
	}
	vector<R3Point>* samples = new vector<R3Point>;
	RNScalar x, y, z;
	while (fin >> x>> y >> z)
	{
		R3Point pos(x, y, z);
		samples->push_back(pos);
	}
	fin.close();
	printf("Get %d samples from %s.\n", samples->size(), filename);
	return samples;	
}

static RNArray<R3MeshVertex*>*
ReadAnchors(R3Mesh* mesh, const char* filename)
{
	ifstream fin(filename);
	if (!fin)
	{
		printf("Unable to open file %s!\n", filename);
		return NULL;
	}
	RNArray<R3MeshVertex*>* anchors = new RNArray<R3MeshVertex*>;
	int idx;
	while (fin >> idx)
	{
		anchors->Insert(mesh->Vertex(idx));
	}
	fin.close();
	printf("Get %d anchors from %s.\n", anchors->NEntries(), filename);
	return anchors;
}

static int
ClosestSampleFromAnchor(vector<R3Point>* samples, const R3Point& anchor)
{
	int idx = -1;
	RNScalar dist = FLT_MAX;
	for (int i=0; i<samples->size(); i++)
	{
		R3Point p0 = (*samples)[i];
		R3Vector v = anchor - p0;
		if (v.Length() < dist)
		{
			dist = v.Length();
			idx = i;
		}
	}
	return idx;
}

static int
ClosestSampleFromAnchor(R3Mesh* mesh, RNArray<R3MeshVertex*>* samples, R3MeshVertex* anchor)
{
	int idx = -1;
	RNScalar dist = FLT_MAX;
	R3Point pos = mesh->VertexPosition(anchor);
	for (int i=0; i<samples->NEntries(); i++)
	{
		R3Point p0 = mesh->VertexPosition(samples->Kth(i));
		R3Vector v = pos - p0;
		if (v.Length() < dist)
		{
			dist = v.Length();
			idx = i;
		}
	}
	return idx;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	g_mesh[0] = ReadMesh(g_input_mesh_name[0]);
	g_mesh[1] = ReadMesh(g_input_mesh_name[1]);
	g_axis[0] = ReadAnchors(g_mesh[0], g_input_axis_name[0]);
	g_axis[1] = ReadAnchors(g_mesh[1], g_input_axis_name[1]);
	g_anchors[0] = ReadAnchors(g_mesh[0], g_input_anchors_name[0]);
	g_anchors[1] = ReadAnchors(g_mesh[1], g_input_anchors_name[1]);

	// find the closest sample to each anchor
	vector<int> anchor_samples[2];
	for (int k=0; k<2; k++)
	{
		for (int i=0; i<g_anchors[k]->NEntries(); i++)
		{
			R3MeshVertex* v = g_anchors[k]->Kth(i);
//			anchor_samples[k].push_back(ClosestSampleFromAnchor(g_axis[k], g_mesh[k]->VertexPosition(v)));
			anchor_samples[k].push_back(ClosestSampleFromAnchor(g_mesh[k], g_axis[k],v));

		}
	}

	vector<int> g_align[2];
	// first mesh: every vertex appears once
	vector<int> nums;

	for (int i=0; i<anchor_samples[0].size(); i++)
	{
		int x = anchor_samples[0][i];
		int y = anchor_samples[0][(i+1)%anchor_samples[0].size()];
		if (x > y)
		{
			y = y + g_axis[0]->NEntries();
		}
		nums.push_back(y - x);
		for (int j=x; j<y; j++)
		{
			g_align[0].push_back(j%g_axis[0]->NEntries());
		}
	}	

	// second mesh: interpolation
	for (int i=0; i<anchor_samples[1].size(); i++)
	{
		int x = anchor_samples[1][i];
		int y = anchor_samples[1][(i+1)%anchor_samples[1].size()];
		if (x > y)
		{
			y = y + g_axis[1]->NEntries();
		}
		for (int j=0; j<nums[i]; j++)
		{
			int idx = x + j*(y-x)/nums[i];
			g_align[1].push_back(idx%g_axis[1]->NEntries());
		}
	}

	// write to file
	ofstream fout(g_output_alignment_name);
	g_cor = new SparseVertexCorrespondence(g_mesh[0], g_mesh[1]);
	for (int i=0; i<g_align[0].size(); i++)
	{
		g_cor->nvertices++;
		g_cor->vertices[0].Insert(g_axis[0]->Kth(g_align[0][i]));
		g_cor->vertices[1].Insert(g_axis[1]->Kth(g_align[1][i]));
	}
	fout.close();
	
	WriteSparseVertexCorrespondence(g_cor, g_output_alignment_name);	
	return 0;
}
