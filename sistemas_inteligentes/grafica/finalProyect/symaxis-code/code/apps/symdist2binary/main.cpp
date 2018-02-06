#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * input_mesh_name = NULL;
char * input_prp_name = NULL;
char * output_binary_name = NULL;
////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * Mesh = NULL;
R3MeshProperty* Property = NULL;
// It should be Max or Min on the center?
// 0 : Max
// 1 : Min
int Center_Is_High = 1;


RNScalar Initial_Threshold = 0;
int Round_Counter = 1;

////////////////////////////////////////////////////////////////////////
// Flags
////////////////////////////////////////////////////////////////////////
int print_verbose = 0;
int print_debug = 0;
int print_time = 0;


////////////////////////////////////////////////////////////////////////
// New Struct
////////////////////////////////////////////////////////////////////////
class Component
{
public:
	RNArray<R3MeshVertex*> * vertices;
	// 1 = Above thres (always in the center)
	// 0 = below thres
	int mark;

	Component()
	{
		vertices = new RNArray<R3MeshVertex*>;
	}
	~Component()
	{
		delete vertices;
	}
};

static void
DeleteComponents(RNArray<Component*>** p_components)
{
	for (int i=0; i<(*p_components)->NEntries(); i++) delete (*p_components)->Kth(i);
	delete (*p_components);
	*p_components = NULL;
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
// Input functions
////////////////////////////////////////////////////////////////////////

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
			else if (!strcmp(*argv, "-low_center"))
			{
				Center_Is_High = 0;
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!input_mesh_name) input_mesh_name = *argv;
			else if (!input_prp_name) input_prp_name = *argv;
			else if (!output_binary_name) output_binary_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!input_mesh_name || !input_prp_name || !output_binary_name) {
    	fprintf(stderr, "Usage: IRSA2Binary mesh prp binary_prp\n");
    	return 0;
  	}

	return 1;
}

// Partition Property given a threshold
static R3MeshProperty*
ThresholdProperty(R3MeshProperty* prp, RNScalar thres)
{
	R3Mesh* mesh = prp->Mesh();
	R3MeshProperty* _prp = new R3MeshProperty(mesh);
	for (int i=0; i<mesh->NVertices(); i++)
	{
		if (prp->VertexValue(i) > thres)
			_prp->SetVertexValue(i, 1.0);
		else
			_prp->SetVertexValue(i, 0.0);
	}
	return _prp;
}

static RNArray<Component*> *
PartitionProperty(R3MeshProperty* binary_prp)
{
	RNArray<Component*>* components = new RNArray<Component*>;
	R3Mesh* mesh = binary_prp->Mesh();
	int nvertices = mesh->NVertices();
	R3mesh_mark ++;

	// Find a vertex not been visited
	int last_visited_idx = -1;	
	while (1)
	{
		last_visited_idx ++;
		// Find the next unvisited vertex
		while (last_visited_idx < nvertices &&
				mesh->VertexMark(mesh->Vertex(last_visited_idx)) == R3mesh_mark) last_visited_idx ++;
		// If all vertices are visited
		if (last_visited_idx >= nvertices) break;

		Component* component = new Component;
		
		// Mark
		RNScalar val = binary_prp->VertexValue(last_visited_idx);
		if (val > 0.5) component->mark = 1;
		else component->mark = 0;

		RNArray<R3MeshVertex*>* array = component->vertices;
		array->Insert(mesh->Vertex(last_visited_idx));
		// Set Mark
		mesh->SetVertexMark(mesh->Vertex(last_visited_idx), R3mesh_mark);

		int pointer = 0;

		while (pointer < array->NEntries())
		{
			R3MeshVertex* vertex = array->Kth(pointer);
			pointer ++;

			for (int neighbor=0; neighbor<mesh->VertexValence(vertex); neighbor++)
			{
				R3MeshEdge* edge = mesh->EdgeOnVertex(vertex, neighbor);
				R3MeshVertex* _v = mesh->VertexAcrossEdge(edge, vertex);
				if (val>0.5 && binary_prp->VertexValue(_v) < 0.5) continue;
				if (val<0.5 && binary_prp->VertexValue(_v) > 0.5) continue;
			   	if (mesh->VertexMark(_v) == R3mesh_mark) continue;
				array->Insert(_v);
				mesh->SetVertexMark(_v, R3mesh_mark);	
			}
		}
		components->Insert(component);
	}
	
	return components;
}

static RNArea
Area(Component* component, R3Mesh* mesh)
{
	RNArea area = 0;
	RNArray<R3MeshVertex*>* vertices = component->vertices;
	for (int i=0; i<vertices->NEntries(); i++)
	{
		area += mesh->VertexArea(vertices->Kth(i));
	}
	return area;
}

static RNArea
Area(RNArray<R3MeshVertex*>* component, R3Mesh* mesh)
{
	RNArea area = 0;
	RNArray<R3MeshVertex*>* vertices = component;
	for (int i=0; i<vertices->NEntries(); i++)
	{
		area += mesh->VertexArea(vertices->Kth(i));
	}
	return area;
}


static RNArray<Component*> *
PartitionProperty(R3MeshProperty* binary_prp, RNScalar min_area)
{
	RNArray<Component*> * components = PartitionProperty(binary_prp);
	while (1)
	{
		int changed_flag = 0;
		for (int i=0; i<components->NEntries(); i++)
		{
			RNArea area = Area(components->Kth(i), binary_prp->Mesh());
			if (area < min_area)
			{
				changed_flag = 1;
				RNArray<R3MeshVertex*>* c = (components->Kth(i))->vertices;
				for (int j=0; j<c->NEntries(); j++)
				{
					RNScalar val = binary_prp->VertexValue(c->Kth(j));
					if (val == 1.0) val = 0.0;
					else val = 1.0;
					binary_prp->SetVertexValue(c->Kth(j), val);				
				}
			}
		}
		if (!changed_flag) break;
		DeleteComponents(&components);
		components = PartitionProperty(binary_prp);
	}
	// Put the component with mark = 1 at the beginning
	RNArray<Component*> * ordered_components = new RNArray<Component*>;
	for (int i=0; i<components->NEntries(); i++)
	{
		if (components->Kth(i)->mark == 1)
			ordered_components->InsertHead(components->Kth(i));
		else
			ordered_components->InsertTail(components->Kth(i));
	}
	delete components;
	return ordered_components;
}

// Count Number of positive components (in the middle)
static int
CountPositiveComponents(RNArray<Component*>* components)
{
	int counter = 0;
	for (int i=0; i<components->NEntries(); i++)
	{
		if (components->Kth(i)->mark == 1)
			counter ++;
	}
	return counter;
}

// Count Number of Negative components
static int
CountNegativeComponents(RNArray<Component*>* components)
{
	int counter = 0;
	for (int i=0; i<components->NEntries(); i++)
	{
		if (components->Kth(i)->mark == 0)
			counter ++;
	}
	return counter;
}

struct PartitionData
{
	int Percentile;
	int NPositive;
	int NNegative;
	RNScalar Ratio;
	PartitionData(int p, int np, int nn, RNScalar r=0){ Percentile = p; NPositive = np; NNegative = nn; Ratio = r;}
};

// Find the Initial Binary Partition (three part, A C B, A and B are very close to each other)
static RNArray<Component*> *
InitialPartition(R3MeshProperty* prp)
{
	R3Mesh* mesh = prp->Mesh();
	R3MeshProperty* binary_prp = NULL;
	RNArray<Component*> * components = NULL;
	// Linear Search
	int percentile = 1;
	RNArray<PartitionData*> all_partitions;
	while (percentile < 100)
	{
		RNScalar thres = prp->Percentile(percentile);
		binary_prp = ThresholdProperty(prp, thres);
		components = PartitionProperty(binary_prp, 0.10);
		int nnegative = CountNegativeComponents(components);
		int npositive = CountPositiveComponents(components);
		if (print_debug)
			cout<<percentile<<" "<<nnegative<<" "<<npositive<<endl;
		if (nnegative == 2 && npositive == 1)
		{
			RNArea area[2];
			area[0] = Area(components->Kth(1), mesh);
			area[1] = Area(components->Kth(2), mesh);
			// ratio of two nnegative regions
			RNScalar r = min(area[0], area[1])/max(area[0], area[1]);
			all_partitions.Insert(new PartitionData(percentile, npositive, nnegative, r));	
		}
		else
		{
			all_partitions.Insert(new PartitionData(percentile, npositive, nnegative));
		}
		percentile ++;
		delete binary_prp;
		DeleteComponents(&components);
	}

	// Find the critical percentile (on the border of n:p = 2:1 and n:p = 1:1) with the highest ratio
	RNScalar max_ratio = 0.0;
	int max_idx = -1;
	for (int i=0; i<all_partitions.NEntries()-1; i++)
	{
		PartitionData* data = all_partitions.Kth(i);
		PartitionData* _data = all_partitions.Kth(i+1);
		if (data->NNegative != 2 || data->NPositive != 1) continue;
		if (_data->NNegative != 1) continue;
		if (data->Ratio > max_ratio)
		{
			max_ratio = data->Ratio;
			max_idx = i;
		}
	}
	// If unable to find the critical percentile,
	if (max_idx == -1)
	{
		cout<<"Cannot find the critial percentile!"<<endl;
		cout<<"Set Percentile = 30"<<endl;
		percentile = 30;
	}
	else
	{
		percentile = all_partitions.Kth(max_idx)->Percentile;
		if (print_debug)
		{
			RNScalar ratio = all_partitions.Kth(max_idx)->Ratio;
			int p = all_partitions.Kth(max_idx)->Percentile;
			cout<<"Best Initial Partition: Ratio("<<ratio<<") Percentile("<<p<<")"<<endl;
		}
	}

	RNScalar thres = prp->Percentile(percentile);
	binary_prp = ThresholdProperty(prp, thres);
	components = PartitionProperty(binary_prp, 0.02);
	delete binary_prp;
	Initial_Threshold = percentile;
	return components;
}

static void
PrintComponents(const char* prefix, RNArray<Component*>* components)
{
	R3MeshProperty* property = new R3MeshProperty(Mesh);
	for (int i=0; i<components->NEntries(); i++)
	{
		Component* component = components->Kth(i);
		RNScalar mark = 0.5;
		if (component->mark == 1)	mark = 1.0;
		else if (component->mark == -1) mark = 0.0;
		for (int j=0; j<component->vertices->NEntries(); j++)
		{
			R3MeshVertex* vertex = component->vertices->Kth(j);
			property->SetVertexValue(vertex, mark);
		}
	}
	cout<<prefix<<"\t"<<components->NEntries()<<endl;
	char filename[100];
	sprintf(filename, "./%s.val", prefix);
	property->Write(filename);
}

int main(int argc, char ** argv)
{

	if (!ParseArgs(argc, argv)) exit(1);

	Mesh = ReadMesh(input_mesh_name);
	Property = ReadProperty(Mesh, input_prp_name);
	if (!Center_Is_High)
	{
		Property->Negate();
	}
	if (print_verbose)
		cout<<"Read property Ok!"<<endl;
	RNTime start_time;
	start_time.Read();
	RNArray<Component*> * components = InitialPartition(Property);
	if (print_verbose)
	{
		cout<<"Initial Partition Ok!"<<endl;
		cout<<components->NEntries()<<endl;
	}
	R3MeshProperty* property = new R3MeshProperty(Mesh);
	for (int i=0; i<components->NEntries(); i++)
	{
		Component* component = components->Kth(i);
		RNScalar mark = 0.0;
		if (component->mark == 1)
			mark = 1.0;
		for (int j=0; j<component->vertices->NEntries(); j++)
		{
			R3MeshVertex* vertex = component->vertices->Kth(j);
			property->SetVertexValue(vertex, mark);
		}
	}

	if (!Center_Is_High)
	{
		for (int i=0; i<Mesh->NVertices(); i++)
		{
			RNScalar val = property->VertexValue(i);
			val = 1.0 - val;
			property->SetVertexValue(i, val);
		}
	}
	if (print_time)
	{
		cout<<"Time: "<<start_time.Elapsed()<<endl;
	}
	property->Write(output_binary_name);
	return 0;
}
