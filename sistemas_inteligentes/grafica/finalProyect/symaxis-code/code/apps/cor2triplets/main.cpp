#include <vector>
#include <iostream>
#include <fstream>
#include "R3Shapes/R3Shapes.h" 
using namespace std;
////////////////////////////////////////////////////////////////////////
// Program: cor2triplets
// Input 1: mesh0 mesh1
// Input 2: tipscorr
// Output:  axiscorr
// Just a template provided for desiging other programs
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name[2] = {NULL, NULL};
char * g_input_tipscorr_name = NULL;
char * g_input_axiscorr_name = NULL;
char * g_output_corr_name = NULL;
char * g_output_triplets_name = NULL;

struct SparseVertexCorrespondence {
 	R3Mesh *mesh[2];
  	RNArray<R3MeshVertex *> vertices[2];
  	int nvertices;
  	SparseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1) { 
    	mesh[0] = mesh0; 
    	mesh[1] = mesh1; 
    	nvertices = 0;
  	};
};

SparseVertexCorrespondence * g_axiscorr = NULL;
SparseVertexCorrespondence * g_tipscorr = NULL;
SparseVertexCorrespondence * g_totalcorr = NULL;

R3Mesh* g_mesh[2] = {NULL, NULL};
struct Triplet
{
	int x, y, z;
	double w;
};

RNArray<Triplet*> * g_triplets = NULL;

int print_verbose = 0;
int print_debug = 0;

int flag_ttt = 1;
int flag_tta = 1;
int flag_taa = 1;
int flag_aaa = 1;
// Please remove tips very near to the axis
int g_thres_vicinity = 0.1;
// how many samples we get from the axis alignemnt
// handle the case where axis is very short
int g_num_axis_samples = 30; // always adjacent
int g_num_axis_smaples_skip_inverse = 6;
// For tip-axis-axis
// how many triplets we generate for each tip
// -1: all 
// other: nearest 10 samples
// handle the case where axis is very short
int g_num_triplets_each_tip = 5;

// For axis-axis-axis
// always pick triplets composed of adjacent correspondences


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

  // Return success
  return mesh;
}

static int
ParseArgs(int argc, char ** argv)
{
	// Parse arguments
    argc--; argv++;
    while (argc > 0) {
        if ((*argv)[0] == '-') {
            if (!strcmp(*argv, "-v")) print_verbose = 1;
        	else if (!strcmp(*argv, "-debug")) print_debug = 1;
			else if (!strcmp(*argv, "-no_ttt"))
			{
				// no tip+tip+tip
				flag_ttt = 0;
			}
			else if (!strcmp(*argv, "-no_tta"))
			{
				// no tip+tip+axis
				flag_tta = 0;
			}
			else if (!strcmp(*argv, "-no_taa"))
			{
				// no tip+axis+axis
				flag_taa = 0;
			}
			else if (!strcmp(*argv, "-no_aaa"))
			{
				// include axis+axis+axis
				flag_aaa = 0;
			}
			else if (!strcmp(*argv, "-naxissamples"))
			{
				argc--; argv++;
				g_num_axis_samples = atoi(*argv);	
			}
			else if (!strcmp(*argv, "-naxissamples_skip"))
			{
				argc--; argv++;
				g_num_axis_smaples_skip_inverse = atoi(*argv);
			}
			else if (!strcmp(*argv, "-ntriplets_each_tip"))
			{
				argc--; argv++;
				g_num_triplets_each_tip = atoi(*argv);
			}
			else if (!strcmp(*argv, "-tipscorr"))
			{
				argc--; argv++;
				g_input_tipscorr_name = *argv;
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
			if (!g_input_mesh_name[0]) g_input_mesh_name[0] = *argv;
			else if (!g_input_mesh_name[1]) g_input_mesh_name[1] = *argv;
			else if (!g_input_axiscorr_name) g_input_axiscorr_name = *argv;
			else if (!g_output_corr_name) g_output_corr_name = *argv;
			else if (!g_output_triplets_name) g_output_triplets_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_axiscorr_name || !g_output_corr_name || 
			!g_output_triplets_name || !g_input_mesh_name[0] || !g_input_mesh_name[1] ) {
    	fprintf(stderr, "Usage: cor2triplets mesh0 mesh1 axis_corr output_corr output_triplets [-no_ttt -no_tta -no_taa -no_aaa" 
				" -naxissamples -tipscorr -ntriplets_each_tip]\n");
    	return 0;
  	}

	return 1;
}

// Read Correspondences
static SparseVertexCorrespondence *
ReadSparseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1, char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Parse filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .cor)\n", filename);
    return 0;
  }

  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open correspondence file %s\n", filename);
    return NULL;
  }

  // Create correspondence
  SparseVertexCorrespondence *correspondence = new SparseVertexCorrespondence(mesh0, mesh1);
  if (!correspondence) {
    fprintf(stderr, "Unable to allocate correspondence for %s\n", filename); 
    return NULL; 
  }

  // Check filename extension
  if (!strcmp(extension, ".cor")) {
    // Read correspondences
    int id0, id1;
    while (fscanf(fp, "%d%d", &id0, &id1) == (unsigned int) 2) {
      R3MeshVertex *vertex0 = mesh0->Vertex(id0);
      R3MeshVertex *vertex1 = mesh1->Vertex(id1);
      correspondence->vertices[0].Insert(vertex0);
      correspondence->vertices[1].Insert(vertex1);
      correspondence->nvertices++;
    }
  }
  else if (!strcmp(extension, ".vk")) {
    // Read header
    double fdummy;
    char format[256], dummy[256];
    int nmeshes, ncorrespondences, idummy;
    if (fscanf(fp, "%s%s%s%s%s%d%s%s%s%d%s%d%d%d%d%d", dummy, format,  dummy, dummy,  dummy, &nmeshes,
               dummy, dummy,  dummy, &idummy,   dummy, &idummy, &idummy, &idummy, &idummy, &ncorrespondences) != (unsigned int) 16) {
      fprintf(stderr, "Unable to read %s\n", filename);
      return NULL;
    }

    // Read correspondences
    int id0, id1;
    while (fscanf(fp, "%d%lg%d%lg", &id0, &fdummy, &id1, &fdummy) == (unsigned int) 4) {
      R3MeshVertex *vertex0 = mesh0->Vertex(id0);
      R3MeshVertex *vertex1 = mesh1->Vertex(id1);
      correspondence->vertices[0].Insert(vertex0);
      correspondence->vertices[1].Insert(vertex1);
      correspondence->nvertices++;
    }

    // Check number of correspondences
    if (correspondence->nvertices != ncorrespondences) {
      fprintf(stderr, "Mismatching number of correspondences in %s\n", filename);
      return NULL;
    }
  }
  else {
    fprintf(stderr, "Unrecognized correspondence file extension: %s\n", extension);
    return NULL;
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Read sparse correspondences from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Correspondences = %d\n", correspondence->nvertices);
    fflush(stdout);
  }

  // Return correspondence
  return correspondence;
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

struct AxisSamplePair
{
	int idx[2]; // index in AxisCorr
	double value;
};

int CompareAxisSamplePair(const void* a, const void* b)
{
	AxisSamplePair** pa = (AxisSamplePair**)a;
	AxisSamplePair** pb = (AxisSamplePair**)b;
	if ((*pa)->value < (*pb)->value)
		return -1;
	else if ((*pa)->value == (*pb)->value)
		return 0;
	else
		return 1;
}

static int *
Sampling(int total, int num)
{
	int* samples = new int[num];
	for (int i=0; i<num; i++)
	{
		samples[i] = i*total/num;
	}
	return samples;
}

static int *
my_sort(RNScalar* vals, int num)
{
	int * idx = new int[num];
	for (int i=0; i<num; i++)
	{
		idx[i] = i;
	}

	for (int i=0; i<num-1; i++)
	{
		for (int j=0; j<num-1-i; j++)
		{
			if (vals[idx[j]]>vals[idx[j+1]])
			{
				// swap idx[j] and idx[j+1]
				int temp = idx[j];
				idx[j] = idx[j+1];
				idx[j+1] = temp;
			}
		}
	}
	return idx;
}
int main(int argc, char ** argv)
{

	if (!ParseArgs(argc, argv)) exit(1);

	g_mesh[0] = ReadMesh(g_input_mesh_name[0]);
	g_mesh[1] = ReadMesh(g_input_mesh_name[1]);

	SparseVertexCorrespondence* axiscorr 
		= ReadSparseVertexCorrespondence(g_mesh[0], g_mesh[1], g_input_axiscorr_name);
	if (g_num_triplets_each_tip<0)
	{
		g_num_triplets_each_tip = g_num_axis_samples;
	}

	int skipped_axis_samples = 1;
	// Handle short axes: add all correspondences
	g_axiscorr = new SparseVertexCorrespondence(g_mesh[0], g_mesh[1]);
	if (axiscorr->nvertices < g_num_axis_samples)
	{
		if (print_verbose)
		{
			cout<<"Processing short axis alignment ..."<<endl;
			cout<<"Length = "<<axiscorr->nvertices<<endl;
		}
		g_num_axis_samples = axiscorr->nvertices;
		for (int i=0; i<axiscorr->nvertices; i++)
		{
			g_axiscorr->nvertices ++;
			g_axiscorr->vertices[0].Insert(axiscorr->vertices[0].Kth(i));
			g_axiscorr->vertices[1].Insert(axiscorr->vertices[1].Kth(i));
		}
		if (axiscorr->nvertices< g_num_triplets_each_tip)
		{
			g_num_triplets_each_tip = axiscorr->nvertices;
		}
		skipped_axis_samples = g_num_axis_samples/g_num_axis_smaples_skip_inverse;	
		if (print_verbose)
		{
			cout<<"g_num_axis_samples: "<<g_num_axis_samples<<endl;
			cout<<"g_num_triplets_each_tip: "<<g_num_triplets_each_tip<<endl;
			cout<<"skipped_axis_samples: "<<skipped_axis_samples<<endl;
		}
	}
	else
	{
		skipped_axis_samples = g_num_axis_samples/g_num_axis_smaples_skip_inverse;
		if (print_verbose)
		{
			cout<<"Processing normal axis alignment ..."<<endl;
			cout<<"Compress "<<axiscorr->nvertices<<" to "<<g_num_axis_samples<<endl;
		}
		int * samples_idx = Sampling(axiscorr->nvertices, g_num_axis_samples);

		for (int i=0; i<g_num_axis_samples; i++)
		{
			g_axiscorr->nvertices ++;
			g_axiscorr->vertices[0].Insert(axiscorr->vertices[0].Kth(samples_idx[i]));
			g_axiscorr->vertices[1].Insert(axiscorr->vertices[1].Kth(samples_idx[i]));
		}
		delete[] samples_idx;
		if (print_verbose)
		{
			cout<<"g_num_axis_samples: "<<g_num_axis_samples<<endl;
			cout<<"g_num_triplets_each_tip: "<<g_num_triplets_each_tip<<endl;
			cout<<"skipped_axis_samples: "<<skipped_axis_samples<<endl;
		}

	}
	if (skipped_axis_samples == 0)
		skipped_axis_samples = 1;

	if (g_input_tipscorr_name)
	{
		SparseVertexCorrespondence*	tipscorr
	   		= ReadSparseVertexCorrespondence(g_mesh[0], g_mesh[1], g_input_tipscorr_name);
		g_tipscorr = new SparseVertexCorrespondence(g_mesh[0], g_mesh[1]);
		if (tipscorr->nvertices > 1)
		{
			// Remove the tips at the vicinity of the axis
			RNLength* distance_to_axes[2];
			for (int i=0; i<2; i++)
			{
				distance_to_axes[i]
				   	= g_mesh[i]->DijkstraDistances(g_axiscorr->vertices[i], g_thres_vicinity);
			}
			for (int i=0; i<tipscorr->nvertices; i++)
			{
				R3MeshVertex* v[2];
				v[0] = tipscorr->vertices[0].Kth(i);
				v[1] = tipscorr->vertices[1].Kth(i);
				if (distance_to_axes[0][g_mesh[0]->VertexID(v[0])]>g_thres_vicinity
						&& distance_to_axes[1][g_mesh[1]->VertexID(v[1])]>g_thres_vicinity)
				{
					g_tipscorr->nvertices ++;
					g_tipscorr->vertices[0].Insert(v[0]);
					g_tipscorr->vertices[1].Insert(v[1]);
				}
			}
			delete[] distance_to_axes[0];
			if (print_verbose)
			{
				cout<<"Num of tip correspondences off vicinity of axis : "<<g_tipscorr->nvertices<<endl;
			}
		}
		else
		{
			if (print_verbose)
			{
				cout<<"No tip correspondences detected!"<<endl;
			}
		}

	}
	else
	{
		g_tipscorr = new SparseVertexCorrespondence(g_mesh[0], g_mesh[1]);
	}
	int ntips = g_tipscorr->nvertices;
	int naxis = g_axiscorr->nvertices;
	// Exceptional cases
	// no tips
	if (ntips == 0)
	{
		flag_ttt = 0;
		flag_tta = 0;
		flag_taa = 0;
		cout<<"flag_aaa is forced to be on!"<<endl;
		flag_aaa = 1;
	}
	else if (ntips == 1)
	{
		flag_ttt = 0;
		flag_tta = 0;
	}
	else if (ntips == 2)
	{
		flag_ttt = 0;
	}
	
	if (naxis == 0)
	{
		flag_aaa = 0;
		flag_tta = 0;
		flag_taa = 0;
		cout<<"flag_ttt is forced to be on!"<<endl;
		if (ntips >=3)
		{
			flag_ttt = 1;
		}
	}
	else if (naxis == 1)
	{
		flag_aaa = 0;
		flag_taa = 0;
		cout<<"flag_ttt is forced to be on!"<<endl;
		if (ntips >=3)
		{
			flag_ttt = 1;
		}
	}
	else if (naxis == 2)
	{
		flag_aaa = 0;
	}

	g_triplets = new RNArray<Triplet*>;
	// tip + tip + tip
    if (flag_ttt)
	{
		RNTime ttt_time;
		ttt_time.Read();
		int ttt_counter = 0;
		for (int i=0; i<ntips; i++)
		{
			for (int j=i+1; j<ntips; j++)
			{
				for (int k=j+1; k<ntips; k++)
				{
					Triplet* triplet = new Triplet;
					triplet->x = i;
					triplet->y = j;
					triplet->z = k;
					triplet->w = 1.0;
					g_triplets->Insert(triplet);
					ttt_counter ++;
				}
			}
		}

		cout<<"# ttt triplets: "<<ttt_counter<<endl;
		cout<<"    Time: "<<ttt_time.Elapsed()<<endl;
	}

	// tip + tip + axis
	if (flag_tta)
	{
		cout<<"Tip pairs are not supported!"<<endl;
	}

	// tip + axis + axis
	if (flag_taa)
	{
		RNTime taa_time;
		taa_time.Read();
		int taa_counter = 0;
		for (int i=0; i<ntips; i++)
		{
			R3MeshVertex* v[2];
			v[0] = g_tipscorr->vertices[0].Kth(i);
			v[1] = g_tipscorr->vertices[1].Kth(i);
			RNLength* dists[2];
			dists[0] = g_mesh[0]->DijkstraDistances(v[0]);
			dists[1] = g_mesh[1]->DijkstraDistances(v[1]);
			RNScalar * dists_to_axes = new RNScalar[naxis];
			for (int j=0; j<naxis; j++)
			{
				R3MeshVertex* u[2];
				u[0] = g_axiscorr->vertices[0].Kth(j);
				u[1] = g_axiscorr->vertices[1].Kth(j);
				int idx[2];
				idx[0] = g_mesh[0]->VertexID(u[0]);
				idx[1] = g_mesh[1]->VertexID(u[1]);
				dists_to_axes[j] = dists[0][idx[0]]+dists[1][idx[1]];
			}
			// Use the pair of axis correspondences nearest to the tips
			int * idx_sorted = my_sort(dists_to_axes, naxis);

			for (int j=0; j<g_num_triplets_each_tip; j++)
			{
				Triplet* triplet = new Triplet;
				triplet->x = i;
				triplet->y = (idx_sorted[j]+naxis-skipped_axis_samples/2)%naxis+ntips;
				triplet->z = (idx_sorted[j]+skipped_axis_samples/2)%naxis+ntips;
				triplet->w = 1.0;
				g_triplets->Insert(triplet);
				taa_counter ++;
			}
			delete[] idx_sorted;
			delete[] dists_to_axes;
		}
		cout<<"# taa triplets: "<<taa_counter<<endl;
		cout<<"    Time: "<<taa_time.Elapsed()<<endl;

	}

	// axis + axis + axis
	if (flag_aaa)
	{
		RNTime aaa_time;
		aaa_time.Read();
		int aaa_counter = 0;
		for (int i=0; i<naxis; i++)
		{
			Triplet* triplet = new Triplet;
			triplet->x = ntips+i;
			triplet->y = ntips+(i+skipped_axis_samples)%naxis;
			triplet->z = ntips+(i+skipped_axis_samples*2)%naxis;
			triplet->w = 1.0;
			g_triplets->Insert(triplet);
			aaa_counter ++;
		}
		cout<<"# aaa triplets: "<<aaa_counter<<endl;
		cout<<"    Time: "<<aaa_time.Elapsed()<<endl;

	}

	cout<<"Generate "<<g_triplets->NEntries()<<" triplets!"<<endl;

	ofstream fout_corr(g_output_corr_name);
	// tips
	for (int i=0; i<ntips; i++)
	{
		fout_corr<<g_mesh[0]->VertexID(g_tipscorr->vertices[0].Kth(i))<<" "
			<<g_mesh[1]->VertexID(g_tipscorr->vertices[1].Kth(i))<<endl;
	}
	// Axiscorr
	for (int i=0; i<naxis; i++)
	{
		fout_corr<<g_mesh[0]->VertexID(g_axiscorr->vertices[0].Kth(i))<<" "
			<<g_mesh[1]->VertexID(g_axiscorr->vertices[1].Kth(i))<<endl;
	}
	fout_corr.close();

	ofstream fout_triplets(g_output_triplets_name);
	// triplets
	for (int i=0; i<g_triplets->NEntries(); i++)
	{
		fout_triplets<<g_triplets->Kth(i)->x<<" "<<g_triplets->Kth(i)->y<<" "
			<<g_triplets->Kth(i)->z<<" "<<g_triplets->Kth(i)->w<<endl;
	}
	fout_triplets.close();
	return 0;
}
