#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
#include "matrix.h"
#include "munkres.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Program: AxisAlignPlusTipsPair2TipsCorRatio
// Input 1: mesh0 mesh1
// Input 2: tipscorr (in the form of coarse correspondences)
// Input 3: axisalignment
// Input 4: orientation (0 or 1), to tell the correspondences within pairs
// Output: Correspondences of tips (.cor)
// Use the coupled tips (with some singletons) and alignment of axis to 
// find correspondences between tips across meshes
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * input_mesh_name[2] = {NULL, NULL};
char * input_symtips_name[2] = {NULL, NULL};
char * input_axisalign_name = NULL;
char * output_tipscor_name = NULL;

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * Mesh[2] = {NULL, NULL};

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
SparseVertexCorrespondence* AxisAlign = NULL;
SparseVertexCorrespondence* SymTips[2] = {NULL, NULL};
SparseVertexCorrespondence* TipsCor = NULL;

int Orient = 0;
////////////////////////////////////////////////////////////////////////
// Flags
////////////////////////////////////////////////////////////////////////

int print_verbose = 0;
int print_debug = 0;
int flag_remove_tips_on_axis = 0;
int flag_distance_normalized = 0;
double dummy_node_cost = 0.2;
int assign_type = 1;
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

// Sparse Correspondence
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
			correspondence->nvertices ++;
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
    	int count = 0;
    	while (fscanf(fp, "%d%lg%d%lg", &id0, &fdummy, &id1, &fdummy) == (unsigned int) 4) {
      		R3MeshVertex *vertex0 = mesh0->Vertex(id0);
      		R3MeshVertex *vertex1 = mesh1->Vertex(id1);
      		correspondence->vertices[0].Insert(vertex0);
      		correspondence->vertices[1].Insert(vertex1);
      		count++;
			correspondence->nvertices ++;
    	}

    	// Check number of correspondences
    	if (count != ncorrespondences) {
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
    	printf("  # Correspondences = %d\n", correspondence->vertices[0].NEntries());
    	fflush(stdout);
  	}

  	// Return correspondence
  	return correspondence;
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
			else if (!strcmp(*argv, "-orient"))
			{
				argc--; argv++;
				Orient = atoi(*argv);
			}
			else if (!strcmp(*argv, "-remove_tips_on_axis"))
			{
				flag_remove_tips_on_axis = 1;
			}
			else if (!strcmp(*argv, "-distance_normalized"))
			{
				flag_distance_normalized = 1;
			}
			else if (!strcmp(*argv, "-dummy_node_cost"))
			{
				argc--; argv++;
				dummy_node_cost = atof(*argv);	
			}
			else if (!strcmp(*argv, "-assign_type"))
			{
				argc--; argv++;
				if (!strcmp(*argv, "MCN"))
				{
					assign_type = 0;
				}
				else if (!strcmp(*argv, "Linear"))
				{
					assign_type = 1;
				}
				else 
				{
					fprintf(stderr, "Unrecognized assign type: %s", *argv); exit(1);
				}
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!input_mesh_name[0]) input_mesh_name[0] = *argv;
			else if (!input_mesh_name[1]) input_mesh_name[1] = *argv;
			else if (!input_symtips_name[0]) input_symtips_name[0] = *argv;
			else if (!input_symtips_name[1]) input_symtips_name[1] = *argv;
			else if (!input_axisalign_name) input_axisalign_name = *argv;
			else if (!output_tipscor_name) output_tipscor_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!input_mesh_name[0] || !input_mesh_name[1] || !input_symtips_name[0] ||
			!input_symtips_name[1] || !input_axisalign_name ||
			!output_tipscor_name) {
    	fprintf(stderr, "Usage: AxisAlignPlusTipsPair2TipsCor mesh0 mesh1 symtips0 symtips1"
				"axisalign output_tipscor [-v, -debug, -orient int, -remove_tips_on_axis, -distance_normalized, -dummy_node_cost, -assign_type {Linear, MCN}]\n");
    	return 0;
  	}

	return 1;
}

static void
NormalizeVector(RNScalar* vec, int n)
{

	RNScalar max = 0;
	RNScalar min = FLT_MAX;
	for (int i=0; i<n; i++)
	{
		if (vec[i] > max)
			max = vec[i];
		if (vec[i] < min)
			min = vec[i];
	}
	for (int i=0; i<n; i++)
	{
		vec[i] = (vec[i] - min)/(max - min);
	}
}

static RNScalar
DifferenceOfVector(RNScalar* v1, RNScalar* v2, int n)
{
	RNScalar sum = 0;
	for (int i=0; i<n; i++)
	{
		RNScalar diff = v1[i] - v2[i];
		sum += diff * diff;
	}
	sum /= n;
	return sum;
}

static RNScalar
RatioDifference(RNScalar* v1, RNScalar* v2, int n)
{
	RNScalar sum = 0;
	for (int i=0; i<n; i++)
	{
		RNScalar diff = fabs(v1[i] - v2[i]);
		if (v1[i] > 0)
		{
			sum += diff/v1[i];
		}
		else if (diff > 0)
		{
			sum += RN_INFINITY;
		}
	}
	sum /= n;
	return sum;
}


static RNScalar**
CreateSignatures(R3Mesh* mesh, SparseVertexCorrespondence* sym, RNArray<R3MeshVertex*>* axis)
{
	int nsyms = sym->nvertices;
	int naxis = axis->NEntries();
	RNScalar** signatures = new RNScalar*[nsyms];
	for (int i=0; i<nsyms; i++)
		signatures[i] = new RNScalar[naxis];

	// Enumerate each point on axis
	for (int i=0; i<nsyms; i++)
	{
		R3MeshVertex* v[2];
		v[0] = sym->vertices[0].Kth(i);
		v[1] = sym->vertices[1].Kth(i);
		for (int j=0; j<naxis; j++)
		{
			RNLength d[2];
		    d[0] = mesh->DijkstraDistance(v[0], axis->Kth(j));
			d[1] = mesh->DijkstraDistance(v[1], axis->Kth(j));
		//	d[0] = max(d[0], 0.01);
		//	d[1] = max(d[1], 0.01);
			// Reduce the influence of the far-away points
			signatures[i][j] = d[0] + d[1];//1.0/d[0] + 1.0/d[1];
		}
		// Normalize
		if (flag_distance_normalized)
			NormalizeVector(signatures[i], naxis);
	}
	return signatures;
}

static void
PrintSignature(RNScalar* sig, int dim)
{
	cout<<"[";
	for (int i=0; i<dim; i++)
		cout<<sig[i]<<",";
	cout<<"]";
}

static SparseVertexCorrespondence*
FindMutuallyClosest(R3Mesh* mesh0, R3Mesh* mesh1, SparseVertexCorrespondence* sym0, SparseVertexCorrespondence* sym1, 
		RNScalar** signatures0, RNScalar** signatures1, int dim)
{
	SparseVertexCorrespondence* tipscor = new SparseVertexCorrespondence(mesh0, mesh1);
	RNScalar matrix[4][4];
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			matrix[i][j] = FLT_MAX;

	for (int i=0; i<sym0->nvertices; i++)
	{
		RNScalar* sig0 = signatures0[i];
		if (print_debug)
		{
			cout<<"-----------------"<<endl;
			cout<<"Cand 0 : "<<i<<" ";
			PrintSignature(sig0, dim);
			cout<<endl;
			cout<<endl;
		}
		R3MeshVertex* cand0[2];
		cand0[0] = sym0->vertices[0].Kth(i);
		cand0[1] = sym0->vertices[1].Kth(i);
		int isPair0 = (cand0[0] != cand0[1]);

		RNScalar min_diff = FLT_MAX;
		int min_cand1_idx = -1;
		for (int j=0; j<sym1->nvertices; j++)
		{
			RNScalar* sig1 = signatures1[j];
			R3MeshVertex* cand1[2];
			cand1[0] = sym1->vertices[0].Kth(j);
			cand1[1] = sym1->vertices[1].Kth(j);
			int isPair1 = (cand1[0] != cand1[1]);
			if (isPair0 != isPair1)
				continue;
//			RNScalar diff = DifferenceOfVector(sig0, sig1, dim);
			RNScalar diff = RatioDifference(sig0, sig1, dim);
			if (print_debug)
				matrix[i][j] = diff;
			if (diff < min_diff)
			{
				min_diff = diff;
				min_cand1_idx = j;
			}
		}
		if (min_cand1_idx == -1)
		{
			// Cannot find
			continue;
		}
		RNScalar* sig1 = signatures1[min_cand1_idx];
		if (print_debug)
		{
			cout<<"Cand 1: "<<min_cand1_idx<<" ";
			PrintSignature(sig1, dim);
			cout<<endl;
			cout<<min_diff<<endl;
			cout<<endl;
			
   /*         if (i == 1)*/
			//{
				//cout<<"Best Cand1: "<<1<<" ";
				//PrintSignature(signatures1[1], dim);
				//cout<<endl;
				//cout<<endl;
			/*}*/
		}
		R3MeshVertex* cand1[2];
		cand1[0] = sym1->vertices[0].Kth(min_cand1_idx);
		cand1[1] = sym1->vertices[1].Kth(min_cand1_idx);
		int isPair1 = isPair0;

		min_diff = FLT_MAX;
		int min_cand2_idx = -1;
		for (int j=0; j<sym0->nvertices; j++)
		{
			RNScalar* sig2 = signatures0[j];
			R3MeshVertex* cand2[2];
			cand2[0] = sym0->vertices[0].Kth(j);
			cand2[1] = sym0->vertices[1].Kth(j);
			int isPair2 = (cand2[0] != cand2[1]);
			if (isPair1 != isPair2)
				continue;
//			RNScalar diff = DifferenceOfVector(sig1, sig2, dim);
			RNScalar diff = RatioDifference(sig1, sig2, dim);

			if (diff < min_diff)
			{
				min_diff = diff;
				min_cand2_idx = j;
			}
		}
		if (min_cand2_idx == -1)
			continue;
		if (print_debug)
		{
			cout<<"Cand 2: "<<min_cand2_idx<<" ";
			PrintSignature(signatures0[min_cand2_idx], dim);
			cout<<endl;
			cout<<min_diff<<endl;
			cout<<endl;
			cout<<endl;
			cout<<endl;
		}
		if (min_cand2_idx == i)
		{
			// mutually closest
			if (isPair0)
			{
				// Add a pair
				tipscor->nvertices ++;
				tipscor->vertices[0].Insert(cand0[0]);
				tipscor->vertices[1].Insert(cand1[0]);
				tipscor->nvertices ++;
				tipscor->vertices[0].Insert(cand0[1]);
				tipscor->vertices[1].Insert(cand1[1]);
			}
			else
			{
				// Add a singleton
				tipscor->nvertices ++;
				tipscor->vertices[0].Insert(cand0[0]);
				tipscor->vertices[1].Insert(cand1[0]);
			}
		}
	}
	if (print_debug)
	{
		for (int i=0; i<4; i++)
		{
			for (int j=0; j<4; j++)
			{
				printf("%5.3f\t", matrix[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}
	if (print_verbose)
	{
		cout<<"# Tips Correspondences: "<<tipscor->nvertices<<endl;
	}
	return tipscor;
}

static SparseVertexCorrespondence*
FindLinearAssignmentWithDummyNodes(R3Mesh* mesh0, R3Mesh* mesh1, 
		SparseVertexCorrespondence* sym0, SparseVertexCorrespondence* sym1, 
		RNScalar** signatures0, RNScalar** signatures1, int dim)
{
	SparseVertexCorrespondence* tipscor = new SparseVertexCorrespondence(mesh0, mesh1);
	Matrix<double> matrix(sym0->nvertices+sym1->nvertices, sym1->nvertices+sym0->nvertices);
	if (print_debug)
	{
		cout<<"Matrix:"<<endl;
	}
	for (int i=0; i<sym0->nvertices; i++)
	{
		RNScalar * sig0 = signatures0[i];
		for (int j=0; j<sym1->nvertices; j++)
		{
			RNScalar* sig1 = signatures1[j];
//			RNScalar diff = DifferenceOfVector(sig0, sig1, dim);
			RNScalar diff = (RatioDifference(sig0, sig1, dim)+RatioDifference(sig1, sig0, dim))*0.5;
			matrix(i, j) = diff;
			if (print_debug)
			{
				cout<<diff<<" ";
			}
		}
		// Dummy nodes
		for (int j=sym1->nvertices; j<sym1->nvertices+sym0->nvertices; j++)
		{
			matrix(i, j) = dummy_node_cost;
			if (print_debug)
			{
				cout<<matrix(i,j)<<"d ";
			}
		}
		if (print_debug)
		{
			cout<<endl;
		}
	}
	for (int i=sym0->nvertices; i<sym0->nvertices+sym1->nvertices; i++)
	{
		for (int j=0; j<sym1->nvertices+sym0->nvertices; j++)
		{
			matrix(i, j) = dummy_node_cost;
			if (print_debug)
			{
				cout<<matrix(i,j)<<"d ";
			}
		}
		if (print_debug)
		{
			cout<<endl;
		}
	}

	Munkres m;
	m.solve(matrix);
   if (!flag_remove_tips_on_axis)
	{
		cout<<"Do not support tips on axis!"<<endl;
		return NULL;
	}
	for (int i=0; i<sym0->nvertices; i++)
	{
		int rowcount = 0;
		for (int j=0; j<sym1->nvertices; j++)
		{
			if (matrix(i, j) == 0)
			{
				tipscor->nvertices = tipscor->nvertices + 2;
				tipscor->vertices[0].Insert(sym0->vertices[0].Kth(i));
				tipscor->vertices[0].Insert(sym0->vertices[1].Kth(i));
				tipscor->vertices[1].Insert(sym1->vertices[0].Kth(j));
				tipscor->vertices[1].Insert(sym1->vertices[1].Kth(j));
				rowcount ++;
			}
		}
		if (rowcount > 1)
		{
			cout<<"Row "<<i<<" has "<<rowcount<<" columns matched!"<<endl;
		}
	}
//	if (print_verbose)
//	{
		cout<<"# Tips Correspondences: "<<tipscor->nvertices<<endl;
//	}
	return tipscor;
}


static SparseVertexCorrespondence*
FindLinearAssignment(R3Mesh* mesh0, R3Mesh* mesh1, 
		SparseVertexCorrespondence* sym0, SparseVertexCorrespondence* sym1, 
		RNScalar** signatures0, RNScalar** signatures1, int dim)
{
	SparseVertexCorrespondence* tipscor = new SparseVertexCorrespondence(mesh0, mesh1);
	Matrix<double> matrix(sym0->nvertices, sym1->nvertices);
	if (print_debug)
	{
		cout<<"Matrix:"<<endl;
	}
	for (int i=0; i<sym0->nvertices; i++)
	{
		RNScalar * sig0 = signatures0[i];
		for (int j=0; j<sym1->nvertices; j++)
		{
			RNScalar* sig1 = signatures1[j];
//			RNScalar diff = DifferenceOfVector(sig0, sig1, dim);
			RNScalar diff = RatioDifference(sig0, sig1, dim);

			matrix(i, j) = diff;
			if (print_debug)
			{
				cout<<diff<<" ";
			}
		}
		if (print_debug)
		{
			cout<<endl;
		}
	}

	Munkres m;
	m.solve(matrix);
   /* if (!flag_remove_tips_on_axis)*/
	/*{*/
		/*cout<<"Do not support tips on axis!"<<endl;*/
		/*return NULL;*/
	/*}*/
	for (int i=0; i<sym0->nvertices; i++)
	{
		int rowcount = 0;
		for (int j=0; j<sym1->nvertices; j++)
		{
			if (matrix(i, j) == 0)
			{
				tipscor->nvertices = tipscor->nvertices + 2;
				tipscor->vertices[0].Insert(sym0->vertices[0].Kth(i));
				tipscor->vertices[0].Insert(sym0->vertices[1].Kth(i));
				tipscor->vertices[1].Insert(sym1->vertices[0].Kth(j));
				tipscor->vertices[1].Insert(sym1->vertices[1].Kth(j));
				rowcount ++;
			}
		}
		if (rowcount > 1)
		{
			cout<<"Row "<<i<<" has "<<rowcount<<" columns matched!"<<endl;
		}
	}
//	if (print_verbose)
//	{
		cout<<"# Tips Correspondences: "<<tipscor->nvertices<<endl;
//	}
	return tipscor;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	Mesh[0] = ReadMesh(input_mesh_name[0]);
	Mesh[1] = ReadMesh(input_mesh_name[1]);

	SymTips[0] = ReadSparseVertexCorrespondence(Mesh[0], Mesh[0], input_symtips_name[0]);
	SymTips[1] = ReadSparseVertexCorrespondence(Mesh[1], Mesh[1], input_symtips_name[1]);
	if (flag_remove_tips_on_axis)
	{
		for (int i=0; i<2; i++)
		{
			for (int j=0; j<SymTips[i]->nvertices; j++)
			{
				if (SymTips[i]->vertices[0].Kth(j) == SymTips[i]->vertices[1].Kth(j))
				{
					SymTips[i]->nvertices --;
					SymTips[i]->vertices[0].RemoveKth(j);
					SymTips[i]->vertices[1].RemoveKth(j);
					j --;
				}
			}
		}
		cout<<"#TipsPairs in 0: "<<SymTips[0]->nvertices<<endl;
		cout<<"#TipsPairs in 1: "<<SymTips[1]->nvertices<<endl;
	}
	// Swapping
	if (Orient)
	{
		for (int i=0; i<SymTips[1]->nvertices; i++)
		{
			R3MeshVertex* temp[2];
		    temp[0]	= SymTips[1]->vertices[0].Kth(i);
			temp[1] = SymTips[1]->vertices[1].Kth(i);
			SymTips[1]->vertices[0].RemoveKth(i);
			SymTips[1]->vertices[1].RemoveKth(i);
			SymTips[1]->vertices[0].InsertKth(temp[1], i);
			SymTips[1]->vertices[1].InsertKth(temp[0], i);
		}
	}

	AxisAlign = ReadSparseVertexCorrespondence(Mesh[0], Mesh[1], input_axisalign_name);

	RNScalar ** signatures[2];
	RNTime start_sig;
	start_sig.Read();
	if (print_verbose)
	{
		cout<<"Creating signatures ..."<<endl;
	}
	signatures[0] = CreateSignatures(Mesh[0], SymTips[0], &(AxisAlign->vertices[0]));
	signatures[1] = CreateSignatures(Mesh[1], SymTips[1], &(AxisAlign->vertices[1]));
	cout<<"********Timing creating signatures: "<<start_sig.Elapsed()<<" seconds"<<endl;

//	TipsCor = FindMutuallyClosest(Mesh[0], Mesh[1], SymTips[0], SymTips[1], signatures[0], signatures[1], AxisAlign->nvertices);
//	TipsCor = FindLinearAssignment(Mesh[0], Mesh[1], SymTips[0], SymTips[1], signatures[0], signatures[1], AxisAlign->nvertices);
		// mutually closest
	if (assign_type == 0)
	{
		TipsCor = FindMutuallyClosest(Mesh[0], Mesh[1], SymTips[0], SymTips[1], signatures[0], signatures[1], AxisAlign->nvertices);
	}
	// linear assignment
	else if (assign_type == 1)
	{
		cout<<"Linear assginment ( Dummy node cost: "<<dummy_node_cost<<")...."<<endl;
		RNTime start;
		start.Read();
		TipsCor = FindLinearAssignmentWithDummyNodes(Mesh[0], Mesh[1], SymTips[0], SymTips[1], signatures[0], signatures[1], AxisAlign->nvertices);
		cout<<"TipsCor Size : "<<TipsCor->nvertices<<endl;
		cout<<"******** Timing linear assignment : "<<start.Elapsed()<<endl;
	}
	WriteSparseVertexCorrespondence(TipsCor, output_tipscor_name);
	return 0;
}
