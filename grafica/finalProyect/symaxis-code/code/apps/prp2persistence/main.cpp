#include "R3Shapes/R3Shapes.h"
#include "RNExt/RNExt.h"
#include "RegionMerger.h"
#include "ANN/ANNMesh.h"

////////////////////////////////////////////////////////////////////////
// Program     : prp2persistence
// Usage	   : prp2persistence mesh property [-lambda, -points, -persistence, -segments]
// Description : Comopute persistence on each vertex from property
//				 Implementation of paper
//				 "Persistent Heat Signature for Pose-oblivious Matching of Incomplete Models"
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name = NULL;
char * g_input_property_name = NULL;

char * g_output_segments_name = NULL;
char * g_output_persistence_name = NULL;
char * g_output_feature_pts_name = NULL;

RNScalar g_lambda = FLT_MAX;
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
			else if (!strcmp(*argv, "-lambda"))
			{
				argc--; argv++;
				g_lambda = atof(*argv);
			}
			else if (!strcmp(*argv, "-points"))
			{
				argc--; argv++;
				g_output_feature_pts_name = *argv;
			}
			else if (!strcmp(*argv, "-segments"))
			{
				argc--; argv++;
				g_output_segments_name = *argv;
			}
			else if (!strcmp(*argv, "-persistence"))
			{
				argc--; argv++;
				g_output_persistence_name = *argv;
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!g_input_mesh_name) g_input_mesh_name = *argv;
			else if (!g_input_property_name) g_input_property_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name || !g_input_property_name) {
    	fprintf(stderr, "Usage: prp2persistence mesh property [-lambda, -points, -persistence, -segments]\n");
    	return 0;
  	}

	return 1;
}

int WriteVertexVector(R3Mesh* mesh, std::vector<R3MeshVertex*>& points, const char* filename)
{
	ofstream fout(filename);
	for (unsigned i=0; i<points.size(); i++)
	{
		R3MeshVertex* v = points[i];
		fout<<mesh->VertexID(v)<<endl;
	}
	fout.close();
	return 1;
}

static int
WriteSegments(std::vector<int>& segments, const char* filename)
{
	ofstream fout(filename);
	for (unsigned i=0; i<segments.size(); i++)
	{
		fout<<segments[i]<<endl;
	}
	fout.close();
	return 1;
}

static int NumConnectedComponents(R3Mesh *mesh, int * flag_array)
{
	int num = mesh->NVertices();
	for (int i=0; i<num; i++)
		flag_array[i] = -1;

  	// Iterate finding connected components
  	int count = 0;
  	int start = 0;
  	for (;;) {
    	// Find an unmarked face
    	R3MeshVertex *seed = NULL;
    	for (; start < mesh->NVertices(); start++) {
			R3MeshVertex* vertex = mesh->Vertex(start);
      		if (flag_array[start] == -1) {
        		seed = vertex;
        		break;
      		}
    	}

    	// Check if found a new component
    	if (!seed) break;
    	else count++;

    	// Mark connected component 
    	RNArray<R3MeshVertex *> stack;
    	stack.InsertTail(seed);
    	while (!stack.IsEmpty()) {
     		R3MeshVertex *vertex = stack.Tail();
      		stack.RemoveTail();
			flag_array[mesh->VertexID(vertex)] = count - 1;
      		for (int i = 0; i<mesh->VertexValence(vertex); i++) {
				R3MeshEdge* edge = mesh->EdgeOnVertex(vertex, i);
        		R3MeshVertex* neighbor = mesh->VertexAcrossEdge(edge, vertex);
				assert(neighbor);
				if (flag_array[mesh->VertexID(neighbor)] == count - 1) continue;
				assert(flag_array[mesh->VertexID(neighbor)] == -1);
          		stack.InsertTail(neighbor);
     	 	}
    	}
  	}
  	// Return number of connected components
  	return count;
}

static int
Mesh2SingleComponent(R3Mesh* mesh, int start_num_neighbors)
{
	printf("Convert mesh to single-component mesh...\n");
	// initialize component flag
	int * component_flag = new int [mesh->NVertices()];

	int component_num = NumConnectedComponents(mesh, component_flag);
	
	ANNMesh ann_mesh(mesh);

	while (component_num > 1)
	{
		printf("    #Components = %d    ", component_num);
		// connecting
		int count = 0;
		for (int i=0; i<mesh->NVertices(); i++)
		{
			R3MeshVertex* v = mesh->Vertex(i);
			int cv = component_flag[mesh->VertexID(v)];
			RNArray<R3MeshVertex*>* cand = ann_mesh.SearchNearestVertices(v, start_num_neighbors);
			for (int j=0;j <cand->NEntries(); j++)
			{
				R3MeshVertex* neighbor = cand->Kth(j);
				int cn = component_flag[mesh->VertexID(neighbor)];
				if (cv == cn) continue;
				if (mesh->EdgeBetweenVertices(v, neighbor)) continue;
				// create an edge connecting them
				mesh->CreateEdge(v, neighbor);
				count ++;
			}
			delete cand;
		}
		printf(" : creating %d edges! (neighborsize = %d)\n", count, start_num_neighbors);
		start_num_neighbors += 3;
		// re-calculate number of components
		component_num = NumConnectedComponents(mesh, component_flag);
	}
	printf("    Ok!\n");
	return 1;
}


int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	R3Mesh* mesh = ReadMesh(g_input_mesh_name);
	// added for handling non-connected mesh (polygon soup)
	Mesh2SingleComponent(mesh, 3);
	R3MeshProperty* property = ReadProperty(mesh, g_input_property_name);
	RegionMerger region_merger(mesh, property, g_lambda);
	region_merger.Run();
	std::vector<R3MeshVertex*> centrals;
	std::vector<RNScalar> persistence;
	region_merger.GetCentralVertices(centrals, persistence);
	R3MeshProperty* prp_consistence = region_merger.GetPersistence();
	std::vector<int> region_assign;
	region_merger.GetRegionAssignment(region_assign);
	if (g_output_persistence_name)
	{
		prp_consistence->Write(g_output_persistence_name);
	}
	if (g_output_feature_pts_name)
	{
		WriteVertexVector(mesh, centrals, g_output_feature_pts_name);
	}
	if (g_output_segments_name)
	{
		WriteSegments(region_assign, g_output_segments_name);
	}
	return 0;
}
