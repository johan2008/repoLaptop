#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
#include <algorithm>
using namespace std;
// mesh name
char * input_mesh_name[2] = {NULL, NULL};
// axis name
char * input_axis_name[2] = {NULL, NULL};
// feature name
char * input_prp_name[2] = {NULL, NULL};
// missing penalty name
char * input_missing_prp_name[2] = {NULL, NULL};
// correspondence name (output)
char * output_corr_name = NULL;
// Score file
char * output_score = NULL; 
// Length file
char * output_length = NULL;
// Samples file
char * output_samples[2] = {NULL, NULL};
char * output_samples_cor = NULL;
// vertex index correspondence
char * output_index_cor = NULL;
//// Input samples file
//char * input_samples[2] = {NULL, NULL};
// Individual cost
char * output_ind_cost = NULL;
// Feature vectors on samples only
char * output_sample_features[2] = {NULL, NULL};
// Sequence file
//char * output_sequence[2] = {NULL, NULL};

R3Mesh * Mesh[2] = {NULL, NULL};
R3MeshPropertySet * FeatureProperties[2] = {NULL, NULL};
R3MeshProperty* MissingProperties[2] = {NULL, NULL};
// Missing Penalty is very large
static RNScalar MissingPenalty = 1e3;
static RNScalar StretchingPenalty = 0;
static int NProperties = 0;
// To be normalized
// Properties from two meshes normalized at the same time (mutiply and shift by a same value)
static int Normalize_Both_Flag = 0;
// Normalize separately
static int Normalize_Separate_Flag = 0;
// Rule out Properties
static int StartPropertyIdx = 0;
static int EndPropertyIdx = -1;
static int Shifting_Flag = 0;
// Samples
struct Sample
{
	RNScalar * Features;
	RNScalar GapPenalty;
	R3Point    Position;
	R3Mesh *   MyMesh;
	R3MeshVertex*  Vertex[2];
	RNScalar   Alpha;
	int 	   Idx;
};

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

RNArray<Sample*>* SrcAxisSamples;
RNArray<Sample*>* DstAxisSamples;

// Axis
RNArray<R3MeshVertex*>* SrcAxis;
RNArray<R3MeshVertex*>* DstAxis;

int print_verbose = 0;
int print_debug = 0;

// How many samples 
int NSamples = 200;

// Read Descriptor
static RNScalar*
ReadDescriptorsFromArff(R3MeshPropertySet* propertyset, R3MeshVertex* vertex)
{
	int ndescriptors = propertyset->NProperties();
	RNScalar * descriptors = new RNScalar [ndescriptors];
	for (int j=0; j<ndescriptors; j++)
	{
		descriptors[j] = propertyset->Property(j)->VertexValue(vertex);
	}

	return descriptors;
}

////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////
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

// Property
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
/*  	if (print_verbose) {*/
        //printf("Read property from %s ...\n", filename);
        //printf("  Time = %.2f seconds\n", start_time.Elapsed());
        //printf("  # Vertices = %d\n", property->Mesh()->NVertices());
        //fflush(stdout);
  	/*}*/

  	// Return property set
  	return property;
}

static R3MeshPropertySet *
ReadProperties(R3Mesh *mesh, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate properties
  R3MeshPropertySet *properties = new R3MeshPropertySet(mesh);
  if (!properties) {
    fprintf(stderr, "Unable to allocate properties for %s\n", filename);
    return NULL;
  }

  // Read properties from file
  if (!properties->Read(filename)) {
    delete properties;
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read properties from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Properties = %d\n", properties->NProperties());
    printf("  # Vertices = %d\n", properties->Mesh()->NVertices());
    fflush(stdout);
  }

  // Return property set
  return properties;
}

static RNArray<R3MeshVertex*>*
ReadAxis(R3Mesh* mesh, char* filename)
{
	ifstream fin(filename);
	if (!fin)
	{
		cout<<"Cannot open axis file "<<filename<<endl;
		exit(-1);
	}
	RNArray<R3MeshVertex*> * axis = new RNArray<R3MeshVertex*>;
	int temp;
	while (fin>>temp)
	{
		if (temp >= mesh->NVertices())
		{
			printf("mesh nvertices = %d, temp = %d\n", mesh->NVertices(), temp);
			exit(-1);
		}
		axis->Insert(mesh->Vertex(temp));
	}
	fin.close();
	return axis;
}

// Sparse Correspondence
static int
WriteSparseVertexCorrespondence(SparseVertexCorrespondence *correspondence, char *filename)
{
	// Start statistics
  	RNTime start_time;
	start_time.Read();

	// Open file
	FILE *fp = fopen(filename, "w");
	if (!fp) {
		fprintf(stderr, "Unable to open sparse file %s\n", filename);
		return 0;
	}

	// Write correspondences
	for (int i=0; i < correspondence->nvertices; i++) {
		fprintf(fp, "%6d\t%6d\n", correspondence->mesh[0]->VertexID(correspondence->vertices[0].Kth(i)), 
				correspondence->mesh[1]->VertexID(correspondence->vertices[1].Kth(i)));
	}

	// Close file
	fclose(fp);

	// Print statistics
	if (print_verbose) {
    	printf("Wrote sparse correspondences to %s ...\n", filename);
    	printf("  Time = %.2f seconds\n", start_time.Elapsed());
    	printf("  # Correspondences = %d\n", correspondence->nvertices);
    	fflush(stdout);
  	}

	// Return success
	return 1;	
}

////////////////////////////////////////////////////////////////////////
// Parse Arguments
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
			else if (!strcmp(*argv, "-property"))
			{
				argv++; argc--;
				input_prp_name[0] = *argv;
				argv++; argc--;
				input_prp_name[1] = *argv;
			}
			else if (!strcmp(*argv, "-nsamples"))
			{
				argv++; argc--;
				NSamples = atoi(*argv);
			}
			else if (!strcmp(*argv, "-score_file"))
			{
				argv++; argc--;
				output_score = *argv;
			}
			else if (!strcmp(*argv, "-missing_property"))
			{
				argv++; argc--;
				input_missing_prp_name[0] = *argv;
				argv++; argc--;
				input_missing_prp_name[1] = *argv;	
			}
			else if (!strcmp(*argv, "-normalize_both"))
			{
				Normalize_Both_Flag = 1;
				Normalize_Separate_Flag = 0;
			}
			else if (!strcmp(*argv, "-normalize_separate"))
			{
				Normalize_Both_Flag = 0;
				Normalize_Separate_Flag = 1;
			}
			else if (!strcmp(*argv, "-missing"))
			{
				argv++; argc--;
				MissingPenalty = atof(*argv);
			}
			else if (!strcmp(*argv, "-stretching"))
			{
				argv++; argc--;
				StretchingPenalty = atof(*argv);
			}
			else if (!strcmp(*argv, "-start"))
			{
				argv++; argc--;
				StartPropertyIdx = atoi(*argv);
			}
			else if (!strcmp(*argv, "-end"))
			{
				argv++; argc--;
				EndPropertyIdx = atoi(*argv);
			}
			else if (!strcmp(*argv, "-shift"))
			{
				Shifting_Flag = 1;
			}
			else if (!strcmp(*argv, "-length_file"))
			{
				argv++; argc--;
				output_length = *argv;
			}
			else if (!strcmp(*argv, "-samples"))
			{
				argv++; argc--;
				output_samples[0] = *argv;
				argv++; argc--;
				output_samples[1] = *argv;
			}
   /*         else if (!strcmp(*argv, "-input_samples"))*/
			/*{*/
				/*argv++; argc--;*/
				/*input_samples[0] = *argv;*/
				/*argv++; argc--;*/
				/*input_samples[1] = *argv;*/
			/*}*/
			else if (!strcmp(*argv, "-samples_cor"))
			{
				argv++; argc--;
				output_samples_cor = *argv;
			}
			else if (!strcmp(*argv, "-individual_cost"))
			{
				argv++; argc--;
				output_ind_cost = *argv;
			}
			else if (!strcmp(*argv, "-sample_features"))
			{
				argv++; argc--;
				output_sample_features[0] = *argv;
				argv++; argc--;
				output_sample_features[1] = *argv;
			}
			else if (!strcmp(*argv, "-index_cor"))
			{
				argv++; argc--;
				output_index_cor = *argv;
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!input_mesh_name[0]) input_mesh_name[0] = *argv;
			else if (!input_mesh_name[1]) input_mesh_name[1] = *argv;
			else if (!input_axis_name[0]) input_axis_name[0] = *argv;
			else if (!input_axis_name[1]) input_axis_name[1] = *argv;
			else if (!output_corr_name) output_corr_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!input_mesh_name[0] || !input_mesh_name[1] || !input_axis_name[0] ||
		   	!input_axis_name[1] || !output_corr_name 
			|| !input_prp_name[0] || !input_prp_name[1]) {
    	fprintf(stderr, "Usage: axis2align meshname0 meshname1 axis0 axis1 outputcor -property [-start int, -end int, -normalize_both, -normalize_separate, -score_file, -nsamples, -samples, -samples_cor, -individual_cost]\n");
    	return 0;
  	}
	return 1;
}



////////////////////////////////////////////////////////////////////////
// Memory
//////////////////////////////////////////////////////////////////////////

static double**
Create2DArray(int h, int w)
{
	double ** array = new double*[h];
	for (int i=0; i<h; i++)
		array[i] = new double[w];
	return array;
}

static void
Delete2DArray(double** array, int h)
{
	for (int i=0; i<h; i++)
		delete[] array[i];
	delete[] array;
}

static int**
Create2DIntArray(int h, int w)
{
	int ** array = new int*[h];
	for (int i=0; i<h; i++)
		array[i] = new int[w];
	return array;
}

static void
Delete2DIntArray(int** array, int h)
{
	for (int i=0; i<h; i++)
		delete[] array[i];
	delete[] array;
}


////////////////////////////////////////////////////////////////////////
// Difference of Descriptor
////////////////////////////////////////////////////////////////////////
static double
ComputeDescriptorDifference(RNScalar* des1, RNScalar* des2, int n)
{
	double diff = 0;
//	cout<<"[";
	for (int i=StartPropertyIdx; i<EndPropertyIdx; i++)
	{
		double tmp = des1[i] - des2[i];
//		cout<<tmp<<",";
//		cout<<i<<" ["<<des1[i]<<" "<<des2[i]<<"]"<<endl;
		diff += tmp*tmp;
	}
//	getchar();
	diff /= EndPropertyIdx - StartPropertyIdx;
	diff = sqrt(diff);
	return diff;
}

////////////////////////////////////////////////////////////////////////
// Sampling on axis
////////////////////////////////////////////////////////////////////////
static RNArray<Sample*>*
SamplingOnAxis(RNArray<R3MeshVertex*>* axis, R3Mesh* mesh, R3MeshPropertySet* propertyset, int nsamples, R3MeshProperty* missing_prp)
{
	RNArray<Sample*>* samples = new RNArray<Sample*>;
	RNLength* distances = new RNLength[axis->NEntries() + 1];
	R3Point* positions = new R3Point[axis->NEntries() + 1];
	// Positions of all vertices on axis
	for (int i=0; i<axis->NEntries(); i++)
	{
		positions[i] = mesh->VertexPosition(axis->Kth(i));
	}
	positions[axis->NEntries()] = mesh->VertexPosition(axis->Kth(0));

	// Calculate distance of adjacent points on the axis
	distances[0] = 0;
	for (int i=1; i<axis->NEntries()+1; i++)
	{
		R3Vector v = positions[i] - positions[i-1];
		distances[i] = v.Length();
	}
	// Calculate distance from origin
	for (int i=1; i<axis->NEntries()+1; i++)
	{
		distances[i] = distances[i] + distances[i-1];
	}
	RNLength length = distances[axis->NEntries()];
	RNLength unit = length/double(nsamples);
	// Create Samples
	// Current distance
	RNLength dist = 0;
	// Current position in distances
	int pos = 0;
	for (int i=0; i<nsamples; i++)
	{
		while (pos <= axis->NEntries() && distances[pos] <= dist) pos ++;
		// pos - 1 and pos is correct
		if ( pos < 1)
		{
			cout<< "Pos < 1"<<endl;
			exit(-1);
		}
		if ( pos > axis->NEntries() )
		{
			cout<<" pos > axis->NEntries() "<<endl;
			exit(-1);
		}
		Sample* sample = new Sample;
		sample->Vertex[0] = axis->Kth(pos-1);
		sample->Vertex[1] = axis->Kth(pos % axis->NEntries());
		sample->MyMesh = mesh;
		sample->Alpha = (dist - distances[pos-1])/(distances[pos] - distances[pos-1]);
		// Because Sample Sequence may be changed
		sample->Idx   = i;
		samples->Insert(sample);
		// Increment
		dist += unit;
	}
	if (print_debug)
	{
		cout<<"Sample Initialize Ok!"<<endl;
	}

	// Compute Feature and Position, based on Alpha
	int nprp = propertyset->NProperties();
	for (int i=0; i<samples->NEntries(); i++)
	{
		Sample* sample = samples->Kth(i);
		R3Point pos[2];
		pos[0] = mesh->VertexPosition(sample->Vertex[0]);
		pos[1] = mesh->VertexPosition(sample->Vertex[1]);
		sample->Position = (1.0-sample->Alpha)*pos[0] + sample->Alpha*pos[1];

		sample->Features = new RNScalar[nprp];
		RNScalar* feat[2];
		feat[0] = ReadDescriptorsFromArff(propertyset, sample->Vertex[0]);
		feat[1] = ReadDescriptorsFromArff(propertyset, sample->Vertex[1]);
		for (int j=0; j<nprp; j++)
		{
			sample->Features[j] = (1.0-sample->Alpha)*feat[0][j] + sample->Alpha*feat[1][j];
		}
		delete[] feat[0];
		delete[] feat[1];	
		// compute missing penalty
		RNScalar gaps_penalties[2];
		gaps_penalties[0] = missing_prp->VertexValue(sample->Vertex[0]);
		gaps_penalties[1] = missing_prp->VertexValue(sample->Vertex[1]);
		sample->GapPenalty = (1.0-sample->Alpha)*gaps_penalties[0] + sample->Alpha*gaps_penalties[1];
	}
	if (print_debug)
	{
		cout<<"Features and Position Ok!"<<endl;
	}
	delete[] distances;
	return samples;
}



////////////////////////////////////////////////////////////////////////
// Compute Pairwise Distance
////////////////////////////////////////////////////////////////////////
static double ***
ComputePairwiseEuclideanDistance(RNArray<Sample*>* samples)
{
	double *** table;
	table = new double**;
	int n = samples->NEntries();
	*table = Create2DArray(n, n);
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
		{
			R3Point position[2];
			position[0] = samples->Kth(i)->Position;
			position[1] = samples->Kth(j)->Position;
			R3Vector dist = position[0] - position[1];
			(*table)[i][j] = dist.Length();
		}
	}
	return table;
}

static double ***
ComputePairwiseGeodesicDistance(RNArray<Sample*>* samples)
{
	double *** table;
	table = new double**;
	int n = samples->NEntries();
	*table = Create2DArray(n, n);
	R3Mesh* mesh = samples->Kth(0)->MyMesh;
	for (int i=0; i<n; i++)
	{
		if (print_verbose)
		{
			cout<<".";
			cout.flush();
		}
		for (int j=0; j<n; j++)
		{
			R3MeshVertex * vertex[2][2];
			Sample* sample[2];
			sample[0] = samples->Kth(i);
			sample[1] = samples->Kth(j);
			vertex[0][0] = sample[0]->Vertex[0];
			vertex[0][1] = sample[0]->Vertex[1];
			vertex[1][0] = sample[1]->Vertex[0];
			vertex[1][1] = sample[1]->Vertex[1];
			RNLength distance[4];
			distance[0] = mesh->DijkstraDistance(vertex[0][0], vertex[1][0]);
			distance[1] = mesh->DijkstraDistance(vertex[0][0], vertex[1][1]);
			distance[2] = mesh->DijkstraDistance(vertex[0][1], vertex[1][0]);
			distance[3] = mesh->DijkstraDistance(vertex[0][1], vertex[1][1]);
			double alpha[2];
			alpha[0] = sample[0]->Alpha;
			alpha[1] = sample[1]->Alpha;
			(*table)[i][j] = (1 - alpha[0])*(1 - alpha[1])*distance[0]
						   + (1 - alpha[0])*alpha[1]*distance[1]
						   + alpha[0]*(1 - alpha[1])*distance[2]
						   + alpha[0]*alpha[1]*distance[3];
		}
	}
	if (print_verbose)
		cout<<endl;
	return table;
}

////////////////////////////////////////////////////////////////////////
// Axis Alignment 
//////////////////////////////////////////////////////////////////////////
static double
MatchAxisOffsetFixed(const RNArray<Sample*>* src, const RNArray<Sample*>* _dst, vector<int>* corrs=NULL, int* fold_over=NULL)
{
	// dst is rotating, and we need to duplicate its head node at the end, to make sure the optimal solution
	RNArray<Sample*>* dst = new RNArray<Sample*>;
	for (int i=0; i<_dst->NEntries(); i++)
	{
		dst->Insert(_dst->Kth(i));
	}
	dst->Insert(_dst->Kth(0));

	int nsrc = src->NEntries();
	int ndst = dst->NEntries();
	double ** table;
	// table
   	table = Create2DArray(nsrc+1, ndst+1);
	int **path;
    path = Create2DIntArray(nsrc+1, ndst+1);
	// 0 - 0 : base
	table[0][0] = 0;//2*ComputeDescriptorDifference(src->Kth(0)->Features, dst->Kth(0)->Features, NProperties);
	path[0][0] = -1;

	// Code for path:
	// 0 : A_{1:i-1}B_{1:j-1} + gamma(a_i, b_j)
	// 1 : A_{1:i-1}B_{1:j} + gamma(a_i, b_j)
	// 2 : A_{1:i-1}B_{1:j} + gamma(a_i,\lambda)
	// 3 : A_{1:i}B_{1:j-1} + gamma(a_i,b_j)
	// 4 : A_{1:i}B_{1:j-1} + gamma(\lambda, b_j)
	
	// 0 - X : src[0] and dst[X] can also be matched
	for (int i=1; i<=ndst; i++)
	{
		table[0][i] = table[0][i-1] + MissingPenalty*dst->Kth(i-1)->GapPenalty;
		path[0][i] = 4;
	}
	// X - 0 : src[X] and dst[0] can also be matched
	for (int i=1; i<=nsrc; i++)
	{
		table[i][0] = table[i-1][0] + MissingPenalty*src->Kth(i-1)->GapPenalty;
		path[i][0] = 2;
	}

	// Other Entries in the Matrix
	for (int i=1; i<=nsrc; i++)
	{
		for (int j=1; j<=ndst; j++)
		{
			double dist = ComputeDescriptorDifference(src->Kth(i-1)->Features, dst->Kth(j-1)->Features, NProperties);

			RNScalar _val[5];
			_val[0] = table[i-1][j-1] + dist;
			_val[1] = table[i-1][j] + dist + StretchingPenalty;// stretching
			_val[2] = table[i-1][j] + MissingPenalty*src->Kth(i-1)->GapPenalty;
			_val[3] = table[i][j-1] + dist + StretchingPenalty;// stretching
			_val[4] = table[i][j-1] + MissingPenalty*dst->Kth(j-1)->GapPenalty;

			RNScalar min_val = FLT_MAX;
			int min_idx = -1;
			for (int k=0; k<5; k++)
			{
				if (_val[k] < min_val)
				{
					min_val = _val[k];
					min_idx = k;
				}
			}

			table[i][j] = min_val;
			path[i][j] = min_idx;
		}
	}
	double ret_value = min(table[nsrc][ndst], table[nsrc][ndst-1]);
	
	if (fold_over)
	{
		if (table[nsrc][ndst] < table[nsrc][ndst-1])
		{
			*fold_over = 1;
		}
		else
		{
			*fold_over = 0;
		}
	}
	// back trace
	if (corrs)
	{
		int idx_i = -1, idx_j = -1;
		RNScalar check_cost = 0;
		int gap_counter = 0;
		corrs[0].clear();
		corrs[1].clear();
		if (table[nsrc][ndst] < table[nsrc][ndst-1])
		{
			idx_i = nsrc;
			idx_j = ndst;
		}
		else
		{
			idx_i = nsrc;
			idx_j = ndst - 1;
		}
		// head and tail of src is mapped to the same node on ndst
		while (1)
		{
			int action = path[idx_i][idx_j];
			if (action == -1)
				break;
			int dj = 0, di = 0;
			int skip = 0;
			if (action == 0)
			{
				di = -1;
				dj = -1;
				skip = 0;
			}
			else if (action == 1)
			{
				di = -1;
				dj = 0;
				skip = 0;
			}
			else if (action == 2)
			{
				di = -1;
				dj = 0;
				skip = 1;
			}
			else if (action == 3)
			{
				di = 0;
				dj = -1;
				skip = 0;
			}
			else if (action == 4)
			{
				di = 0;
				dj = -1;
				skip = 1;
			}
			if (!skip)
			{
				corrs[0].push_back(idx_i - 1);
				corrs[1].push_back(idx_j - 1);
				if (print_debug)
				{
					check_cost += ComputeDescriptorDifference(src->Kth(idx_i-1)->Features, dst->Kth(idx_j-1)->Features, NProperties);
				}
			}
			else
			{
				if (print_debug)
				{
//					check_cost += MissingPenalty;
					gap_counter ++;
//					cout<<"Add penalty!"<<endl;
//					getchar();
				}
			}
			idx_i += di;
			idx_j += dj;
		}
		// change all ndst to 0
		for (int i=0; i<corrs[1].size(); i++)
		{
			if (corrs[1][i] == ndst-1)
			{
				corrs[1][i] = 0;
			}
		}
		if (print_debug)
		{
			cout<<"Check cost = "<<check_cost<<endl;
			cout<<"#Gaps      = "<<gap_counter<<endl;
		}
	}
	Delete2DArray(table, nsrc+1);
	Delete2DIntArray(path, nsrc+1);

	return ret_value;
}

static  double
MatchAxis(RNArray<Sample*>* src, RNArray<Sample*>* dst, vector<int>* align)
{
	int ndst = dst->NEntries();
	double cost = FLT_MAX;
	int best_offset = -1;
    for (int i=0; i<dst->NEntries(); i++)
	{
		int fold_over = 0;
		double _cost = MatchAxisOffsetFixed(src, dst, NULL, &fold_over);
//		if (print_debug)
//		{
//			char f = '*';
//			if (!fold_over)
//				f = ' ';
//			printf("Cost (offset = %3d): %4.3f%c\n", i, _cost, f);
//		}
		if (_cost < cost)
		{
			cost = _cost;
			best_offset = i;
		}
		Sample* sample = dst->Kth(0);
		dst->RemoveHead();
		dst->Insert(sample);
	}
	if (print_debug)
	{
		cout<<"Best offset: "<<best_offset<<endl;
	}
	for (int i=0; i<best_offset; i++)
	{
		Sample* sample =dst->Kth(0);
		dst->RemoveHead();
		dst->Insert(sample);
	}
	cost = MatchAxisOffsetFixed(src, dst, align);
	for (int i=0; i<align[1].size(); i++)
	{
		align[1][i] = (align[1][i] + best_offset)%dst->NEntries();
	}
	// recover to the correct order
	for (int i=best_offset; i<dst->NEntries(); i++)
	{
		Sample* sample = dst->Kth(0);
		dst->RemoveHead();
		dst->Insert(sample);
	}
   /* if (print_verbose)*/
	/*{*/
		/*cout<<"#Raw Alignment = "<<align[0].size()<<endl;*/
	/*}*/
   /* if (print_debug)*/
	/*{*/
		/*cout<<"Alignment"<<endl;*/
		/*for (int i=0; i<align[0].size(); i++)*/
		/*{*/
			/*cout<<"["<<align[0][i]<<" , "<<align[1][i]<<"]"<<"\tDist: ["<<ComputeDescriptorDifference(src->Kth(align[0][i])->Features, dst->Kth(align[1][i])->Features, NProperties)<<"]"<<endl;*/
   /*[>        for (int kk=0; kk<NProperties; kk++)<]*/
			/*[>{<]*/
				/*[>cout<<"      ";<]*/
				/*[>cout<<kk<<" : ";<]*/
				/*[>cout<<src->Kth(align[0][i])->Features[kk]<<" "<<src->Kth(align[1][i])->Features[kk]<<endl;<]*/
			/*[>}<]*/
			/*[>getchar();<]*/
		/*}*/
		/*cout<<endl;*/
	/*}*/
	return cost;
}

static int
WriteCorrespondence(RNArray<Sample*>* src, RNArray<Sample*>* dst, int* match, char* filename)
{
	int n = src->NEntries();
	R3Mesh* srcMesh = src->Kth(0)->MyMesh;
	R3Mesh* dstMesh = dst->Kth(0)->MyMesh;
	RNArray<R3MeshVertex*> srcVertices;
	RNArray<R3MeshVertex*> dstVertices;
	for (int i=0; i<n; i++)
	{
		if (match[i] == -1) continue;
		Sample* srcSample = src->Kth(i);
		Sample* dstSample = dst->Kth(match[i]);
		R3MeshVertex* srcVertex = (srcSample->Alpha<0.5)?srcSample->Vertex[0]:srcSample->Vertex[1];
		R3MeshVertex* dstVertex = (dstSample->Alpha<0.5)?dstSample->Vertex[0]:dstSample->Vertex[1];
		if ((!srcVertices.IsEmpty() && srcVertices.Tail() == srcVertex) 
				|| (!dstVertices.IsEmpty() && dstVertices.Tail() == dstVertex))
			continue;
		srcVertices.Insert(srcVertex);
		dstVertices.Insert(dstVertex);
	}
	// Check the last one
	if (srcVertices.Tail() == srcVertices.Head() || dstVertices.Tail() == dstVertices.Head())
	{
		srcVertices.RemoveTail();
		dstVertices.RemoveTail();
	}
	// Write to file
	ofstream fout(filename);
	for (int i=0; i<srcVertices.NEntries(); i++)
	{
		fout<<srcMesh->VertexID(srcVertices.Kth(i))<<" ";
		fout<<dstMesh->VertexID(dstVertices.Kth(i))<<endl;
	}
	fout.close();
	return 1;
}

// Write Correspondence Sequence (in Axis ID)
/*static void*/
//WriteSequence(int* match, char* filename, int n)
//{
	//ofstream fout(filename);
	//for (int i=0; i<n; i++)
	//{
		//int src = i;
		//int dst = match[i];
		//if (dst>=0)
			//fout<<src<<" "<<dst<<endl;
	//}
	//fout.close();
/*}*/

static void
NormalizeScalarsByWeights(RNScalar* val0, RNScalar* val1, RNScalar* weight0, RNScalar* weight1, int dim0, int dim1)
{
	// put in one array
	RNScalar* values = new RNScalar[dim0 + dim1];
	RNScalar* weights= new RNScalar[dim0 + dim1];
	for (int i=0; i<dim0; i++)
	{
		values[i] = val0[i];
		weights[i]= weight0[i];
	}
	for (int i=0; i<dim1; i++)
	{
		values[i+dim0] = val1[i];
		weights[i+dim0]= weight1[i];
	}
	// mean
	RNScalar mean = 0;
	RNScalar w = 0;
	for (int i=0; i<dim0+dim1; i++)
	{
		mean += values[i]*weights[i];
		w += weights[i];
	}
	if (w > 1e-9)
		mean /= w;

	// Standard Deviation
	RNScalar dev = 0;
	for (int i=0; i<dim0+dim1; i++)
	{
		dev += (values[i] - mean)*(values[i] - mean)*weights[i];
	}
	if (w > 1e-9)
		dev /= w;
	dev = sqrt(dev);

	// finalize
	if (dev > 1e-9)
	{
		for (int i=0; i<dim0+dim1; i++)
		{
			values[i] = (values[i] - mean)/dev;
		}
	}
	// Copy back to two arrays
	for (int i=0; i<dim0; i++)
		val0[i] = values[i];
	for (int i=0; i<dim1; i++)
		val1[i] = values[i + dim0];

   /* if (print_debug)*/
	/*{*/
		/*cout<<"mean: "<<mean<<endl;*/
		/*cout<<"Dev:  "<<dev<<endl;*/
	/*}*/
	delete[] weights;
	delete[] values;
}

static void
NormalizeProperties(R3MeshPropertySet* prp0, R3MeshPropertySet* prp1)
{
	R3Mesh* mesh[2];
	mesh[0] = prp0->mesh;
	mesh[1] = prp1->mesh;
	int dim[2];
	dim[0] = mesh[0]->NVertices();
	dim[1] = mesh[1]->NVertices();
	RNScalar* vals[2];
	vals[0] = new RNScalar[dim[0]];
	vals[1] = new RNScalar[dim[1]];
	RNScalar* weights[2];
	weights[0] = new RNScalar[dim[0]];
	weights[1] = new RNScalar[dim[1]];
	// Weights decided by vertex area
	for (int i=0; i<dim[0]; i++)
	{
		R3MeshVertex* vertex = mesh[0]->Vertex(i);
		weights[0][i] = 1.0;//= mesh[0]->VertexArea(vertex);
	}
	for (int i=0; i<dim[1]; i++)
	{
		R3MeshVertex* vertex = mesh[1]->Vertex(i);
		weights[1][i] = 1.0;//= mesh[1]->VertexArea(vertex);
	}

	// Normalize each property
	int nproperties = prp0->NProperties();
	for (int i=0; i<nproperties; i++)
	{
		R3MeshProperty* property[2];
		property[0] = prp0->Property(i);
		property[1] = prp1->Property(i);
		// Copy property to vals
		for (int j=0; j<dim[0]; j++)
			vals[0][j] = property[0]->VertexValue(j);
		for (int j=0; j<dim[1]; j++)
			vals[1][j] = property[1]->VertexValue(j);
		NormalizeScalarsByWeights(vals[0], vals[1], weights[0], weights[1], dim[0], dim[1]);
		// Copy vals back to property
		for (int j=0; j<dim[0]; j++)
			property[0]->SetVertexValue(j, vals[0][j]);
		for (int j=0; j<dim[1]; j++)
			property[1]->SetVertexValue(j, vals[1][j]);
	}

	delete[] weights[0];
	delete[] weights[1];
	delete[] vals[0];
	delete[] vals[1];
}

static SparseVertexCorrespondence*
RoughCorrs2SparseVertexCorrespondence(RNArray<Sample*>* src, RNArray<Sample*>* dst, vector<int>* align)
{
//	SparseVertexCorrespondence* corrs = new SparseVertexCorrespondence(src->Kth(0)->MyMesh, dst->Kth(0)->MyMesh);
	SparseVertexCorrespondence* corrs = new SparseVertexCorrespondence(Mesh[0], Mesh[1]);
	RNArray<Sample*> sample_corrs[2];
	int idx = 0;
	int nalign = align[0].size();
//	for (int i=0; i<nalign; i++)
//	{
//		cout<<i<<". ("<<align[0][i]<<","<<align[1][i]<<")"<<endl;
//	}
	while (idx<nalign)
	{
		if (idx+1 < nalign && align[0][idx] == align[0][idx+1])
		{
			// x cluster
			if (align[1][idx] == align[1][idx+1])
			{
				cout<<"Both idx_i and idx_j are same"<<endl;
				exit(-1);
			}
			int counter = 1;
			while (idx+1 < nalign && align[0][idx] == align[0][idx+1])
			{
				idx ++;
				counter ++;
			}
			sample_corrs[0].Insert(src->Kth(align[0][idx]));
			sample_corrs[1].Insert(dst->Kth(align[1][idx - counter/2]));
			idx ++;
		}
		else if (idx+1 < nalign && align[1][idx] == align[1][idx+1])
		{
			// y cluster
			if (align[0][idx] == align[0][idx+1])
			{
				cout<<"Both idx_j and idx_i are same"<<endl;
				exit(-1);
			}
			int counter = 1;
			while (idx+1 < nalign && align[1][idx] == align[1][idx+1])
			{
				idx ++;
				counter ++;
			}
			sample_corrs[0].Insert(src->Kth(align[0][idx-counter/2]));
			sample_corrs[1].Insert(dst->Kth(align[1][idx]));
			idx ++;
		}
		else
		{
			// single pair
			sample_corrs[0].Insert(src->Kth(align[0][idx]));
			sample_corrs[1].Insert(dst->Kth(align[1][idx]));
			idx ++;
		}
	}
	if (print_debug)
	{
		cout<<"Center samples : "<<sample_corrs[0].NEntries()<<endl;
	}
	RNArray<R3MeshVertex*>& srcVertices = corrs->vertices[0];
	RNArray<R3MeshVertex*>& dstVertices = corrs->vertices[1];
	for (int i=0; i<sample_corrs[0].NEntries(); i++)
	{
		Sample* srcSample = sample_corrs[0].Kth(i);
		Sample* dstSample = sample_corrs[1].Kth(i);
		R3MeshVertex* srcVertex = (srcSample->Alpha<0.5)?srcSample->Vertex[0]:srcSample->Vertex[1];
		R3MeshVertex* dstVertex = (dstSample->Alpha<0.5)?dstSample->Vertex[0]:dstSample->Vertex[1];
		if ((!srcVertices.IsEmpty() && srcVertices.Tail() == srcVertex)
				|| (!dstVertices.IsEmpty() && dstVertices.Tail() == dstVertex))
			continue;
		srcVertices.Insert(srcVertex);
		dstVertices.Insert(dstVertex);
		corrs->nvertices++;
	}
	// Check the last one
	if (srcVertices.Tail() == srcVertices.Head() || dstVertices.Tail() == dstVertices.Head())
	{
		srcVertices.RemoveTail();
		dstVertices.RemoveTail();
		corrs->nvertices--;
	}

	return corrs;
}


static void
WriteScore2File(RNScalar score, char* filename)
{
	ofstream fout(filename);
	fout<<score<<endl;
	fout.close();
}

static void
WriteLength2File(int length, char* filename)
{
	ofstream fout(filename);
	fout<<length<<endl;
	fout.close();
}

static void
NormalizeProperty(R3MeshProperty* property)
{
	R3Mesh* mesh = property->Mesh();
	RNScalar mean = 0;
	for (int i=0; i<mesh->NVertices(); i++)
	{
		R3MeshVertex* v = mesh->Vertex(i);
		mean += property->VertexValue(v) * mesh->VertexArea(v);
	}

	RNScalar std = 0;
	for (int i=0; i<mesh->NVertices(); i++)
	{
		R3MeshVertex* v = mesh->Vertex(i);
		RNScalar diff = property->VertexValue(v) - mean;
		std += diff * diff * mesh->VertexArea(v);
	}
	std = sqrt(std);
	if (std < 1e-4)
	{
		for (int i=0; i<mesh->NVertices(); i++)
		{
			R3MeshVertex* v = mesh->Vertex(i);
			property->SetVertexValue(v, 0);
		}
	}
	else
	{
		for (int i=0; i<mesh->NVertices(); i++)
		{
			R3MeshVertex* v = mesh->Vertex(i);
			RNScalar val = property->VertexValue(v);
			val = (val - mean)/std;
			property->SetVertexValue(v, val);
		}
	}
}

// Write samples to files
static int
WriteSamples(char * filename, RNArray<Sample*>* samples)
{
	ofstream fout(filename);
	for (int i=0; i<samples->NEntries(); i++)
	{
		Sample* sample = samples->Kth(i);
		fout<<sample->Position.X()<<" "<<sample->Position.Y()<<" "<<sample->Position.Z()<<endl;
	}
	fout.close();
	return 1;
}

// Write correspondences of samples to files
static int
WriteSampleCorrespondences(char* filename,vector<int>* corrs)
{
	ofstream fout(filename);
	for (int i=0; i<corrs[0].size(); i++)
	{
		fout<<corrs[0][i]<<" "<<corrs[1][i]<<endl;
	}
	fout.close();
	return 1;
}

static int
WriteSampleFeatures(char* filename, RNArray<Sample*>* samples)
{
	ofstream fout(filename);
	fout<<samples->NEntries()<<" ";
	fout<<NProperties<<endl;
	for (int i=0; i<samples->NEntries(); i++)
	{
		Sample* sample = samples->Kth(i);
		for (int j=StartPropertyIdx; j<EndPropertyIdx; j++)
		{
			fout<<sample->Features[j]<<" ";
		}
		fout<<endl;
	}
	fout.close();
	return 1;
}


int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);

	if (print_verbose)
	{
		cout<<"NSamples = "<<NSamples<<endl;
	}
	// Read Mesh
	Mesh[0] = ReadMesh(input_mesh_name[0]);
	Mesh[1] = ReadMesh(input_mesh_name[1]);

	// Read Axis
	SrcAxis = ReadAxis(Mesh[0], input_axis_name[0]);
	DstAxis = ReadAxis(Mesh[1], input_axis_name[1]);

	// Read Propertyset
	FeatureProperties[0] = ReadProperties(Mesh[0], input_prp_name[0]);
	FeatureProperties[1] = ReadProperties(Mesh[1], input_prp_name[1]);
	NProperties = FeatureProperties[0]->NProperties();

	// Read missing penalties (thin part has lower penalty, easier to be cut off)
	if (input_missing_prp_name[0] && input_missing_prp_name[1])
	{
		MissingProperties[0] = ReadProperty(Mesh[0], input_missing_prp_name[0]);
		MissingProperties[1] = ReadProperty(Mesh[1], input_missing_prp_name[1]);
		for (int i=0; i<Mesh[0]->NVertices(); i++)
		{
			if (MissingProperties[0]->VertexValue(i) < 0.4)
			{
				MissingProperties[0]->SetVertexValue(i, 0.0);
			}
			else
			{
				MissingProperties[0]->SetVertexValue(i, 100.0);
			}
		}
		for (int i=0; i<Mesh[1]->NVertices(); i++)
		{
			if (MissingProperties[1]->VertexValue(i) < 0.4)
			{
				MissingProperties[1]->SetVertexValue(i, 0.0);
			}
			else
			{
				MissingProperties[1]->SetVertexValue(i, 100.0);
			}
		}


	}
	else
	{
		MissingProperties[0] = new R3MeshProperty(Mesh[0]);
		for (int i=0; i<Mesh[0]->NVertices(); i++)
		{
			MissingProperties[0]->SetVertexValue(i, 1.0);
		}
		MissingProperties[1] = new R3MeshProperty(Mesh[1]);
		for (int i=0; i<Mesh[1]->NVertices(); i++)
		{
			MissingProperties[1]->SetVertexValue(i, 1.0);
		}
	}
	if (EndPropertyIdx < 0)
	{
		EndPropertyIdx = NProperties;
	}
	if (print_debug)
		cout<<"Read Property Ok!"<<endl;

	// Normalize
	RNTime prp_norm_time;
	prp_norm_time.Read();
	if (Normalize_Both_Flag)
	{
		if (print_verbose)
		{
			cout<<"Need to normalize descriptors together!"<<endl;
		}
		if (print_verbose)
			cout<<"Normalizing ..."<<endl;
		NormalizeProperties(FeatureProperties[0], FeatureProperties[1]);
	}
	else if (Normalize_Separate_Flag)
	{
		if (print_verbose)
		{
			cout<<"Need to normalize descriptors separately!"<<endl;
		}
		if (print_verbose)
			cout<<"Normalizing ..."<<endl;
		for (int i=0; i<NProperties; i++)
		{
			NormalizeProperty(FeatureProperties[0]->Property(i));
			NormalizeProperty(FeatureProperties[1]->Property(i));
		}
		if (print_debug)
		{
			FeatureProperties[0]->Write("1.arff");
			FeatureProperties[1]->Write("2.arff");
		}
	}
	printf("[Normalizing property time : %.2f seconds]\n", prp_norm_time.Elapsed());
//	if (Shifting_Flag)
//	{
//		if (print_verbose)
//		{
//			cout<<"Shift one propertyset for testing!"<<endl;
//		}
//		// Shift FeatureProperties[0]
//		for (int i=0; i<FeatureProperties[0]->NProperties(); i++)
//		{
//			R3MeshProperty* prp = FeatureProperties[0]->Property(i);
//			for (int j=0; j<Mesh[0]->NVertices(); j++)
//			{
//				RNScalar value = prp->VertexValue(j);
//				value *= 2.0;
//				prp->SetVertexValue(j, value);
//			}
//		}	
//	}

	// Start statistics
  	RNTime sampling_start_time;
  	sampling_start_time.Read();

	// Samples
	if (print_verbose)
		cout<<"Sampling on Axes ..."<<endl;
	SrcAxisSamples = SamplingOnAxis(SrcAxis, Mesh[0], FeatureProperties[0], NSamples, MissingProperties[0]);
	DstAxisSamples = SamplingOnAxis(DstAxis, Mesh[1], FeatureProperties[1], NSamples, MissingProperties[1]);

	if (output_sample_features[0] && output_sample_features[1])
	{
		WriteSampleFeatures(output_sample_features[0], SrcAxisSamples);
		WriteSampleFeatures(output_sample_features[1], DstAxisSamples);
	} 

	// output samples
	if (output_samples[0] && output_samples[1])
	{
		WriteSamples(output_samples[0], SrcAxisSamples);
		WriteSamples(output_samples[1], DstAxisSamples);
	}
	printf("[Sampling Time (2 axes) = %.2f seconds]\n", sampling_start_time.Elapsed());

	double Cost;
	vector<int> RoughCorrs[2];	
	if (print_verbose)
	{
		cout<<"Match Axis ..."<<endl;
	}
	RNTime match_start_time;
	match_start_time.Read();
    Cost = MatchAxis(SrcAxisSamples, DstAxisSamples, RoughCorrs);
	cout<<">>>>>>>>>>>>>>>>>>>>>>> Total Cost : "<<Cost<<endl;
	printf("[Matching Time = %.2f seconds]\n", match_start_time.Elapsed());

	if (output_samples_cor)
	{
		WriteSampleCorrespondences(output_samples_cor, RoughCorrs);
	}

	// if needs to output inidividual costs
	if (output_ind_cost)
	{
		vector<double> ind_cost;
		double sum = 0;
		for (int i=0; i<RoughCorrs[0].size(); i++)
		{
			int x = RoughCorrs[0][i];
			int y = RoughCorrs[1][i];
			Sample* sample[2];
			sample[0] = SrcAxisSamples->Kth(x);
			sample[1] = DstAxisSamples->Kth(y);
			double diff = ComputeDescriptorDifference(sample[0]->Features, sample[1]->Features, NProperties);
			ind_cost.push_back(diff);
			sum += diff;
		}
		cout<<"Sum of individual costs = "<<sum<<endl;
		ofstream fout(output_ind_cost);
		for (int i=0; i<ind_cost.size(); i++)
		{
			fout<<ind_cost[i]<<endl;
		}
		fout.close();
	}


	RNTime generate_corr_start_time;
	generate_corr_start_time.Read();
	SparseVertexCorrespondence* FinalCorrs = NULL;
	FinalCorrs = RoughCorrs2SparseVertexCorrespondence(SrcAxisSamples, DstAxisSamples, RoughCorrs);
	printf("Generating Corr Time = %.2f seconds\n", generate_corr_start_time.Elapsed());

	cout<<"NCorrespondence : "<<FinalCorrs->nvertices<<endl;
	if (output_length)
		WriteLength2File(FinalCorrs->nvertices, output_length);

	if (output_score)
		WriteScore2File(Cost, output_score);
	// Write to files
	WriteSparseVertexCorrespondence(FinalCorrs, output_corr_name);
	if (output_index_cor)
	{
		vector<int> index_cor[2];
		for (int i=0; i<FinalCorrs->nvertices; i++)
		{
			R3MeshVertex* v[2];
			v[0] = FinalCorrs->vertices[0][i];
			v[1] = FinalCorrs->vertices[1][i];
			int idx[2] = {-1, -1};
			for (int j=0; j<SrcAxis->NEntries(); j++)
			{
				if (v[0] == SrcAxis->Kth(j))
				{
					idx[0] = j;
					break;
				}
			}
			for (int j=0; j<DstAxis->NEntries(); j++)
			{
				if (v[1] == DstAxis->Kth(j))
				{
					idx[1] = j;
					break;
				}
			}
			if (idx[0] < 0 || idx[1] < 0)
			{
				cout<<"Unable to find vertex in the input axis!"<<endl;
				break;
			}
			index_cor[0].push_back(idx[0]);
			index_cor[1].push_back(idx[1]);
		}
		ofstream fout(output_index_cor);
		for (int i=0; i<int(index_cor[0].size()); i++)
		{
			fout<<index_cor[0][i]<<" "<<index_cor[1][i]<<endl;
		}
		fout.close();
	}
	return 0;
}

