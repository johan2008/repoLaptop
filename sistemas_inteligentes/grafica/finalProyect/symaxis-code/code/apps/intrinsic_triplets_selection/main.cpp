#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Program: intrinsic_triplets_selection2
// Input 1: mesh tips property (agd)
// Input 2: output_tips output_triplets
// Pruning criteria:
// 	d(z,w)/d(m(z),m(w)) < 0.9
// 	d(z,m(z)) < sqrt(area(M)/32PI)
// 	agd(z)/agd(w) < 0.9
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name = NULL;
char * g_input_tips_name = NULL;
char * g_output_tips_name = NULL;
char * g_output_triplets_name = NULL;
char * g_input_prp_name = NULL;

char * g_output_my_triplets_name = NULL;
char * g_output_my_coarse_name = NULL;
int max_triplets = 1000;
// Triplet
struct Triplet
{
	int idx[3];
};

// Pair of triplets
struct ConformalMap
{
	int idx[6];
};

RNScalar g_thres_dist_diff = 0.9;
RNScalar g_thres_dist = -1;
RNScalar g_thres_prp_diff = 0.9;

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * g_mesh = NULL;
R3MeshProperty* g_property = NULL;
RNArray<R3MeshVertex*>* g_tips = NULL;
RNScalar * g_tip_features = NULL;
RNScalar ** g_distances = NULL;
RNArray<ConformalMap*>* g_conformal_maps = NULL;
////////////////////////////////////////////////////////////////////////
// Flags
////////////////////////////////////////////////////////////////////////

int print_verbose = 0;
int print_debug = 0;

int flag_abc = 1;
int flag_abcd = 0;
int flag_abcde = 1;

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
  	if (print_verbose) {
    	printf("Read property from %s ...\n", filename);
    	printf("  Time = %.2f seconds\n", start_time.Elapsed());
    	printf("  # Vertices = %d\n", property->Mesh()->NVertices());
    	fflush(stdout);
  	}

  	// Return property set
  	return property;
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
			else if (!strcmp(*argv, "-my_triplets"))
			{
				argc--; argv++;
				g_output_my_triplets_name = *argv;
			}
			else if (!strcmp(*argv, "-my_coarse"))
			{
				argc--; argv++;
				g_output_my_coarse_name = *argv;
			}
			else if (!strcmp(*argv, "-property"))
			{
				argc--; argv++;
				g_input_prp_name = *argv;
			}
			else if (!strcmp(*argv, "-thres_dist_diff"))
			{
				argc--; argv++;
				g_thres_dist_diff = atof(*argv);
			}
			else if (!strcmp(*argv, "-thres_dist"))
			{
				argc--; argv++;
				g_thres_dist = atof(*argv);
			}
			else if (!strcmp(*argv, "-thres_prp_diff"))
			{
				argc--; argv++;
				g_thres_prp_diff = atof(*argv);
			}
			else if (!strcmp(*argv, "-abc"))
			{
				flag_abc = 1;
			}
			else if (!strcmp(*argv, "-abcd"))
			{
				flag_abcd = 1;
			}
			else if (!strcmp(*argv, "-abcde"))
			{
				flag_abcde = 1;
			}
			else if (!strcmp(*argv, "-max_triplets"))
			{
				argc--; argv++;
				max_triplets = atoi(*argv);
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!g_input_mesh_name) g_input_mesh_name = *argv;
			else if (!g_input_tips_name) g_input_tips_name = *argv;
			else if (!g_output_tips_name) g_output_tips_name = *argv;
			else if (!g_output_triplets_name) g_output_triplets_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name || !g_input_tips_name || !g_output_tips_name || !g_output_triplets_name || !g_input_prp_name) {
    	fprintf(stderr, "Usage:intrinsic_triplets_selection mesh tips output_tips output_triplets -property prp [-v, -debug, -my_triplets, -my_coarse] \n");
    	return 0;
  	}

	return 1;
}

static RNScalar*
GetTipsFeatures(RNArray<R3MeshVertex*>* tips, R3MeshProperty* property)
{
	RNScalar* features = new RNScalar[tips->NEntries()];
	for (int i=0; i<tips->NEntries(); i++)
	{
		features[i] = property->VertexValue(tips->Kth(i));
	}
	return features;
}

static double**
Create2DArray(int h, int w)
{
	double ** array = new double*[h];
	for (int i=0; i<h; i++)
		array[i] = new double[w];
	return array;
}


static RNScalar**
GetDistanceMatrix(RNArray<R3MeshVertex*>* points, R3Mesh* mesh)
{
	int ntips = points->NEntries();
	RNScalar** matrix = Create2DArray(ntips, ntips);
	for (int i=0; i<ntips; i++)
	{
		matrix[i][i] = 0;
		for (int j=i+1; j<ntips; j++)
		{
			matrix[i][j] = mesh->DijkstraDistance(points->Kth(i), points->Kth(j));
			matrix[j][i] = matrix[i][j];
		}
	}
	if (print_debug)
	{
		for (int i=0; i<ntips; i++)
		{
			for (int j=0; j<ntips; j++)
			{
				cout<<matrix[i][j]<<" ";
			}
			cout<<endl;
		}
	}
	return matrix;
}

// Write to triplets
static void
WriteTriplets(RNArray<ConformalMap*>* maps, const char* filename)
{
	ofstream fout(filename);
	for (int i=0; i<maps->NEntries(); i++)
	{
		ConformalMap* cm = maps->Kth(i);
		fout<<cm->idx[0]<<" "<<cm->idx[1]<<" "
			<<cm->idx[2]<<" "<<cm->idx[3]<<" "
			<<cm->idx[4]<<" "<<cm->idx[5]<<" ";
		fout<<1;
		if (i != maps->NEntries()-1)
			fout<<endl;
	}
	fout.close();
}

// Write to tips_var
static void
WriteTips(RNArray<R3MeshVertex*>* points, R3Mesh* mesh, char* filename)
{
	ofstream fout(filename);
	for (int i=0; i<points->NEntries(); i++)
	{
		fout<<mesh->VertexID(points->Kth(i))<<endl;
	}
	fout.close();
}

static int
GenerateABCPairs()
{
	int num_conformal_maps = g_conformal_maps->NEntries();
	RNTime time_abc;
	time_abc.Read();
	for (int i=0; i<g_tips->NEntries(); i++)
	{
		for (int j=0; j<g_tips->NEntries(); j++)
		{
			for (int k=j+1; k<g_tips->NEntries(); k++)
			{
				if (i == j || j == k || k == i)
				{
					continue;
				}
				// 	d(z,m(z)) < sqrt(area(M)/32PI)
				if (g_distances[j][k] < g_thres_dist)
				{
					continue;
				}
				// 	d(z,w)/d(m(z),m(w)) < 0.9
				if (g_distances[i][j]/g_distances[i][k]<g_thres_dist_diff ||
						g_distances[i][k]/g_distances[i][j]<g_thres_dist_diff)
				{
					continue;
				}
				//	agd(z)/agd(w) < 0.9
				if (g_tip_features[j] < 1e-4 || g_tip_features[k]<1e-4)
				{
					continue;
				}
				if (g_tip_features[j]/g_tip_features[k]<g_thres_prp_diff ||
						g_tip_features[k]/g_tip_features[j]<g_thres_prp_diff)
				{
					continue;
				}
				ConformalMap* cm[2];
			    cm[0] = new ConformalMap;
				cm[0]->idx[0] = i;
				cm[0]->idx[1] = i;
				cm[0]->idx[2] = j;
				cm[0]->idx[3] = k;
				cm[0]->idx[4] = k;
				cm[0]->idx[5] = j;
				g_conformal_maps->Insert(cm[0]);
   /*             cm[1] = new ConformalMap;*/
				//cm[1]->idx[0] = i;
				//cm[1]->idx[1] = i;
				//cm[1]->idx[2] = k;
				//cm[1]->idx[3] = j;
				//cm[1]->idx[4] = j;
				//cm[1]->idx[5] = k;
				/*g_conformal_maps->Insert(cm[1]);*/

			}
		}
	}
	cout<<"Number of (a,b,c)(a,c,b) type :"<<g_conformal_maps->NEntries() - num_conformal_maps<<endl;
	cout<<"    Time: "<<time_abc.Elapsed()<<" seconds. "<<endl;
	return 1;
}

static int
GenerateABCDPairs()
{
	int num_conformal_maps = g_conformal_maps->NEntries();
	RNTime time_abcd;
	time_abcd.Read();
	for (int i=0; i<g_tips->NEntries(); i++)
	{
		// a
		for (int j=i+1; j<g_tips->NEntries(); j++)
		{
			// b
			// check the similarity of a and b
			if (g_distances[i][j] < g_thres_dist)
			{
				continue;
			}
			if (g_tip_features[i] < 1e-4 || g_tip_features[j] < 1e-4)
			{
				continue;
			}
			if (g_tip_features[i]/g_tip_features[j]<g_thres_prp_diff ||
					g_tip_features[j]/g_tip_features[i]<g_thres_prp_diff)
			{
				continue;
			}
			for (int p=0; p<g_tips->NEntries(); p++)
			{
				if (p == i || p == j)
				{
					continue;
				}
				// c
				for (int q=p+1; q<g_tips->NEntries(); q++)
				{
					if (q == i || q == j)
					{
						continue;
					}
					// d
					// check the similarity of c and d
					if (g_distances[p][q] < g_thres_dist)
					{
						continue;
					}
					if (g_tip_features[p]<1e-4 || g_tip_features[q]<1e-4)
					{
						continue;
					}
					if (g_tip_features[p]/g_tip_features[q]<g_thres_prp_diff ||
							g_tip_features[q]/g_tip_features[p]<g_thres_prp_diff)
					{
						continue;
					}
					// ac vs bd
					if (g_distances[i][p]/g_distances[j][q] < g_thres_dist_diff ||
							g_distances[j][q]/g_distances[i][p] < g_thres_dist_diff)
					{
						continue;
					}
					// ad vs bc
					if (g_distances[i][q]/g_distances[j][p] < g_thres_dist_diff ||
							g_distances[j][p]/g_distances[i][q] < g_thres_dist_diff)
					{
						continue;
					}
					ConformalMap* cm[2];
				    cm[0] = new ConformalMap;
					cm[0]->idx[0] = i;
					cm[0]->idx[1] = j;
					cm[0]->idx[2] = j;
					cm[0]->idx[3] = i;
					cm[0]->idx[4] = p;
					cm[0]->idx[5] = q;
					g_conformal_maps->Insert(cm[0]);
   /*                 cm[1] = new ConformalMap;*/
					//cm[1]->idx[0] = i;
					//cm[1]->idx[1] = j;
					//cm[1]->idx[2] = p;
					//cm[1]->idx[3] = q;
					//cm[1]->idx[4] = q;
					//cm[1]->idx[5] = p;
					/*g_conformal_maps->Insert(cm[1]);*/
				}
			}
		}
	}
	cout<<"Number of a<->b, c<->d type :"<<g_conformal_maps->NEntries() - num_conformal_maps<<endl;
	cout<<"    Time: "<<time_abcd.Elapsed()<<" seconds. "<<endl;
	return 1;
}

static int
GenerateABCDEPairs()
{
	int num_conformal_maps = g_conformal_maps->NEntries();
	RNTime time_abcde;
	time_abcde.Read();
	for (int a=0; a<g_tips->NEntries(); a++)
	{
		for (int b=0; b<g_tips->NEntries(); b++)
		{
			if (b == a) continue;
			for (int d=0; d<g_tips->NEntries(); d++)
			{
				if (d == a || d == b) continue;
				// (1) distance of bd
				if (g_distances[b][d] < g_thres_dist)
					continue;
				// (2) distance ratio of ab/ad
				if (g_distances[a][b]/g_distances[a][d]<g_thres_dist_diff 
						|| g_distances[a][d]/g_distances[a][b]<g_thres_dist_diff)
					continue;
				// (3) feature value ratio of f(b)/f(d)
				if (g_tip_features[b]/g_tip_features[d]<g_thres_prp_diff 
						|| g_tip_features[d]/g_tip_features[b]<g_thres_prp_diff)
					continue;

				for (int c=b+1; c<g_tips->NEntries(); c++)
				{
					if (c == a || c == d) continue;
					for (int e=0; e<g_tips->NEntries(); e++)
					{
						if (e == d || e == c || e == b || e == a) continue;
						// a-b <==> a-d, a-c <==> a-e
						// (1) distance of bd, and ce
						if (g_distances[c][e] < g_thres_dist)
							continue;
						// (2) distance ratio of ab/ad and ac/ae
						if (g_distances[a][c]/g_distances[a][e]<g_thres_dist_diff
								|| g_distances[a][e]/g_distances[a][c]<g_thres_dist_diff)
							continue;
						// (3) feature value ratio of f(b)/f(d) and f(c)/f(e)
						if (g_tip_features[c]/g_tip_features[e]<g_thres_prp_diff
								|| g_tip_features[e]/g_tip_features[c]<g_thres_prp_diff)
							continue;
						ConformalMap* cm[2];
					    cm[0] = new ConformalMap;
						cm[0]->idx[0] = a;
						cm[0]->idx[1] = a;
						cm[0]->idx[2] = b;
						cm[0]->idx[3] = d;
						cm[0]->idx[4] = c;
						cm[0]->idx[5] = e;
						g_conformal_maps->Insert(cm[0]);
   /*                     cm[1] = new ConformalMap;*/
						//cm[1]->idx[0] = b;
						//cm[1]->idx[1] = d;
						//cm[1]->idx[2] = e;
						//cm[1]->idx[3] = c;
						//cm[1]->idx[4] = c;
						//cm[1]->idx[5] = e;
						/*g_conformal_maps->Insert(cm[1]);*/

					}
				}
			}
		}
	}
	cout<<"Number of (a,b,c)(a,d,e) type :"<<g_conformal_maps->NEntries() - num_conformal_maps<<endl;
	cout<<"    Time: "<<time_abcde.Elapsed()<<" seconds. "<<endl;
	return 1;
}

static void
TruncateConformalMaps(RNArray<ConformalMap*>* maps, int max_num)
{
	if (maps->NEntries()<max_num)
		return;
	RNArray<ConformalMap*> temp_array;
	for (int i=0; i<max_num; i++)
	{
		int idx = i*maps->NEntries()/max_num;
		temp_array.Insert(maps->Kth(idx));
	}

	maps->Empty();
	for (int i=0; i<temp_array.NEntries(); i++)
	{
		maps->Insert(temp_array.Kth(i));
	}
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);

	g_mesh = ReadMesh(g_input_mesh_name);
	g_tips = ReadPoints(g_mesh, g_input_tips_name);
	g_property = ReadProperty(g_mesh, g_input_prp_name);

	if (g_thres_dist < 0)
	{
		g_thres_dist = sqrt(g_mesh->Area()/(32.0*RN_PI));
	}
	if (print_verbose)
	{
		cout<<"g_thres_dist      = "<<g_thres_dist<<endl;
		cout<<"g_thres_dist_diff = "<<g_thres_dist_diff<<endl;
		cout<<"g_thres_prp_diff  = "<<g_thres_prp_diff<<endl;
	}
	
	g_tip_features = GetTipsFeatures(g_tips, g_property);
	g_distances = GetDistanceMatrix(g_tips, g_mesh);

	cout<<"#Tips = "<<g_tips->NEntries()<<endl;
	g_conformal_maps = new RNArray<ConformalMap*>;
	// Pair of triplets like: (a,b,c) (a,c,b)
	if (flag_abc)
	{
		GenerateABCPairs();
	}
	// Pair of triplets like: a<->b, c<->d
	if (flag_abcd)
	{
		GenerateABCDPairs();
	}
	// Pair of triplets like: (a,b,c) (a,d,e)
	if (flag_abcde)
	{
		GenerateABCDEPairs();
	}
	cout<<"# Passing triplet : "<<g_conformal_maps->NEntries()<<endl;
	// if there is only one triplet, duplicate it (Vova's program cannot handle one triplet)
	if (g_conformal_maps->NEntries() == 1)
	{
		ConformalMap* cm = new ConformalMap;
		for (int i=0; i<6; i++)
			cm->idx[i] = g_conformal_maps->Kth(0)->idx[i];
		g_conformal_maps->Insert(cm);
		cout<<"We duplicate the first triplet (because there is only one triplet)!"<<endl;
	}
	else if (g_conformal_maps->NEntries() == 0)
	{
		cout<<"No plausible triplet!"<<endl;
	}

	if (g_conformal_maps->NEntries()>max_triplets)
	{
		cout<<"Truncate number of triplets ..."<<endl;
		TruncateConformalMaps(g_conformal_maps, max_triplets);
		cout<<"Remaining # conformal maps : "<<g_conformal_maps->NEntries()<<endl;
	}

	WriteTriplets(g_conformal_maps, g_output_triplets_name);
	WriteTips(g_tips, g_mesh, g_output_tips_name);
	if (g_output_my_triplets_name || g_output_my_coarse_name)
	{
		ofstream fout0(g_output_my_triplets_name);
		ofstream fout1(g_output_my_coarse_name);
		int counter = 0;
		for (int i=0; i<g_conformal_maps->NEntries(); i++)
		{
			ConformalMap* cm = g_conformal_maps->Kth(i);
			int a = cm->idx[0];
			int b = cm->idx[1];
			int c = cm->idx[2];
			int d = cm->idx[3];
			int e = cm->idx[4];
			int f = cm->idx[5];
			fout0<<counter<<" "<<counter+1<<" "<<counter+2<<" "<<1<<endl;
			fout1<<g_mesh->VertexID(g_tips->Kth(a))<<" "<<g_mesh->VertexID(g_tips->Kth(b))<<endl;
			fout1<<g_mesh->VertexID(g_tips->Kth(c))<<" "<<g_mesh->VertexID(g_tips->Kth(d))<<endl;
			fout1<<g_mesh->VertexID(g_tips->Kth(e))<<" "<<g_mesh->VertexID(g_tips->Kth(f))<<endl;
			counter += 3;
		}
		fout0.close();
		fout1.close();
	}
	return 0;
}
