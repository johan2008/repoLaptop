#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Program: alignedmesh2halveskdtree
// Input 1: mesh
// Output:  dist_x dist_y dist_z
// Use kd-tree to accelerate the nearest neighbor search
// ////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name = NULL;
char * g_input_map_name = NULL;
char * g_input_axis_name = NULL;
char * g_output_map_name = NULL;

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////
struct DenseVertexCorrespondence {
  	R3Mesh *mesh[2];
  	R3MeshVertex **vertices; 
  	DenseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1) { 
    	mesh[0] = mesh0; 
    	mesh[1] = mesh1; 
    	vertices = new R3MeshVertex * [ mesh0->NVertices() ];
    	for (int i = 0; i < mesh0->NVertices(); i++) vertices[i] = NULL;
  	};
};
////////////////////////////////////////////////////////////////////////
// Flags
////////////////////////////////////////////////////////////////////////

int print_verbose = 0;
int print_debug = 0;

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

// Dense Correspondence
static DenseVertexCorrespondence *
ReadDenseVertexCorrespondence(R3Mesh* mesh0, R3Mesh* mesh1, char* filename)
{
	// Create Dense Correspondence
	DenseVertexCorrespondence* map = new DenseVertexCorrespondence(mesh0, mesh1);
    
	// Open file
	ifstream fin(filename);
	if (!fin)
	{
		printf("Cannot open file for reading dense correspondence : %s\n", filename);
		delete map;
		return NULL;
	}
	int temp;
	int counter = 0;
	while (fin >> temp)
	{
		if (counter >= mesh0->NVertices())
		{
			printf("Map file does not match to Mesh : %s\n", filename);
			delete map;
			fin.close();
			return NULL;
		}
		map->vertices[counter] = mesh1->Vertex(temp);
		counter ++;
	}
    
	if (counter != mesh0->NVertices())
	{
		printf("Map file does not match Mesh : %s\n", filename);
		delete map;
		fin.close();
		return NULL;
	}
	
	fin.close();
	return map;
}

static int
WriteDenseVertexCorrespondence(DenseVertexCorrespondence *map, char *filename)
{
  	// Start statistics
  	RNTime start_time;
  	start_time.Read();

  	// Open file
  	FILE *fp = fopen(filename, "w");
  	if (!fp) {
    	fprintf(stderr, "Unable to open map file %s\n", filename);
    	return 0;
	}

	// Write correspondences
	for (int i = 0; i < map->mesh[0]->NVertices(); i++) {
		if (map->vertices[i] == NULL) fprintf(fp, "-1\n");
    	else fprintf(fp, "%d\n", map->mesh[1]->VertexID(map->vertices[i]));
  	}

  	// Close file
  	fclose(fp);

  	// Print statistics
  	if (print_verbose) {
    	printf("Wrote dense correspondences to %s ...\n", filename);
    	printf("  Time = %.2f seconds\n", start_time.Elapsed());
    	printf("  # Correspondences = %d\n", map->mesh[0]->NVertices());
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
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!g_input_mesh_name) g_input_mesh_name = *argv;
            else if (!g_input_axis_name) g_input_axis_name = *argv;
            else if (!g_input_map_name) g_input_map_name = *argv;
            else if (!g_output_map_name) g_output_map_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name || !g_input_axis_name || !g_input_map_name
        || !g_output_map_name) {
    	fprintf(stderr, "Usage: symmaprefine mesh axis org_map new_map [-v, -debug]\n");
    	return 0;
  	}

	return 1;
}

static R3MeshVertex* 
NearestNeighbor(R3MeshSearchTree* search_tree, R3Point query)
{
    R3Mesh* mesh = search_tree->mesh;
	R3MeshIntersection closest;
    search_tree->FindClosest(query, closest, 0, 1.0);
	if (closest.type == R3_MESH_NULL_TYPE)
	{
        assert(FALSE);
	}

	RNLength closest_distance_to_surface = closest.t;
	R3Point closest_point_on_surface = closest.point;

	// Determine closest vertex
	R3MeshVertex *closest_vertex = NULL;
   	RNLength closest_distance_to_vertex = FLT_MAX;
   	if (closest.type == R3_MESH_VERTEX_TYPE) {
    	 // Closest point was on vertex
     	closest_vertex = closest.vertex;
     	closest_distance_to_vertex = closest.t;
   	}
   	else if (closest.type == R3_MESH_EDGE_TYPE) {
     	// Closest point was on edge
     	R3MeshVertex *vertex0 = mesh->VertexOnEdge(closest.edge, 0);
     	RNLength d0 = R3Distance(query, mesh->VertexPosition(vertex0));
     	if (d0 < closest_distance_to_vertex) {
       		closest_vertex = vertex0;
       		closest_distance_to_vertex = d0;
     	}
     	R3MeshVertex *vertex1 = mesh->VertexOnEdge(closest.edge,1);
     	RNLength d1 = R3Distance(query, mesh->VertexPosition(vertex1));
     	if (d1 < closest_distance_to_vertex) {
       		closest_vertex = vertex1;
       		closest_distance_to_vertex = d1;
     	}
   	}
   	else if (closest.type == R3_MESH_FACE_TYPE) {
     	// Closest point was in middle of face
     	R3MeshVertex *vertex0 = mesh->VertexOnFace(closest.face, 0);
     	RNLength d0 = R3Distance(query, mesh->VertexPosition(vertex0));
     	if (d0 < closest_distance_to_vertex) {
       		closest_vertex = vertex0;
       		closest_distance_to_vertex = d0;
     	}
     	R3MeshVertex *vertex1 = mesh->VertexOnFace(closest.face,1);
     	RNLength d1 = R3Distance(query, mesh->VertexPosition(vertex1));
     	if (d1 < closest_distance_to_vertex) {
       		closest_vertex = vertex1;
       		closest_distance_to_vertex = d1;
     	}
     	R3MeshVertex *vertex2 = mesh->VertexOnFace(closest.face,2);
     	RNLength d2 = R3Distance(query, mesh->VertexPosition(vertex2));
     	if (d2 < closest_distance_to_vertex) {
       		closest_vertex = vertex2;
       		closest_distance_to_vertex = d2;
     	}
   	}
	return closest_vertex;
}

static int*
CreateSideFunction(R3Mesh* mesh, RNArray<R3MeshVertex*>* axis)
{
	int * sides = new int [mesh->NVertices()];
	for (int i=0; i<mesh->NVertices(); i++) sides[i] = -1;
    
	RNArray<R3MeshEdge*> axis_edges;
	for (int i=0; i<axis->NEntries(); i++)
	{
		R3MeshVertex* v[2];
		v[0] = axis->Kth(i);
		v[1] = axis->Kth((i+1)%axis->NEntries());
        
		R3MeshEdge* edge = mesh->EdgeBetweenVertices(v[0], v[1]);
		if (edge)
			axis_edges.Insert(edge);
		else
		{
			cout<<"Axis is disconnected"<<endl;
			exit(-1);
		}
	}
    
	// Prepare : set all vertices to be initial
	R3mesh_mark ++;
	int init_mark = R3mesh_mark;
	for (int i=0; i<mesh->NVertices(); i++)
		mesh->SetVertexMark(mesh->Vertex(i), init_mark);
	// Preparation: Mark the vertices on axis
	R3mesh_mark ++;
	int axis_mark = R3mesh_mark;
	for (int i=0; i<axis->NEntries(); i++)
	{
		mesh->SetVertexMark(axis->Kth(i), axis_mark);
	}
	R3mesh_mark ++;
	int cw_mark = R3mesh_mark;
	R3mesh_mark ++;
	int ccw_mark = R3mesh_mark;
	int idx = 0;
	while (idx < axis->NEntries())
	{
		R3MeshVertex* seed = NULL;
		R3MeshVertex* vertex = axis->Kth(idx);
		R3MeshEdge* edge = axis_edges.Kth(idx);
		int visited_counter = 0;
		while (1)
		{
			edge = mesh->EdgeOnVertex(vertex, edge, RN_CCW);
			if (!edge) {/*idx++;continue;*/break;}
			seed = mesh->VertexAcrossEdge(edge, vertex);
			if (!seed) {/*idx++;continue;*/break;}
			if (mesh->VertexMark(seed) == axis_mark) // another axis edge
			{ /*idx++;continue;*/ break; }
			if (mesh->VertexMark(seed) == ccw_mark) // alread marked
			{/* idx++;*/continue;}
			else if (mesh->VertexMark(seed) == cw_mark)
			{
				cout<<"ccw to cw conflict (init)!"<<endl;
				cout<<idx<<endl;
				if (print_debug)
					getchar();
				/*idx++;*/
				/*continue;*/
				break;
			}
			RNArray<R3MeshVertex*> stack;
			stack.Insert(seed);
			mesh->SetVertexMark(seed, ccw_mark);
            
			while (!stack.IsEmpty())
			{
				R3MeshVertex* vertex = stack.Tail();
				stack.RemoveTail();
                //				cout<<mesh->VertexID(vertex)<<endl;
                //				visited_counter ++;
				for (int i=0; i<mesh->VertexValence(vertex); i++)
				{
					R3MeshEdge* edge = mesh->EdgeOnVertex(vertex, i);
					R3MeshVertex* neighbor = mesh->VertexAcrossEdge(edge, vertex);
					if (!neighbor) continue;
					if (mesh->VertexMark(neighbor) == ccw_mark)
						continue;
					if (mesh->VertexMark(neighbor) == cw_mark)
					{
						cout<<"ccw to cw conflict!"<<endl;
						cout<<idx<<endl;
						if (print_debug)
							getchar();
					}
					// set new mark
					if (mesh->VertexMark(neighbor) == init_mark)
					{
						stack.Insert(neighbor);
						mesh->SetVertexMark(neighbor, ccw_mark);
					}
				}
			}
		}
		// cw
		vertex = axis->Kth(idx);
		edge = axis_edges.Kth(idx);
		while (1)
		{
			edge = mesh->EdgeOnVertex(vertex, edge, RN_CW);
			if (!edge) {/*idx++;continue;*/break;}
			seed = mesh->VertexAcrossEdge(edge, vertex);
			if (!seed) {/*idx++;continue;*/break;}
			if (mesh->VertexMark(seed) == axis_mark) // another axis edge
			{/*idx++;continue;*/break;}
			if (mesh->VertexMark(seed) == cw_mark) // alread marked
			{/*idx++;*/continue;}
			else if (mesh->VertexMark(seed) == ccw_mark)
			{
				cout<<"cw to ccw conflict (init)!"<<endl;
				cout<<idx<<endl;
				if (print_debug)
					getchar();
				/*idx++;
                 continue;*/
				break;
			}
			RNArray<R3MeshVertex*> stack;
			stack.Empty();
			stack.Insert(seed);
			mesh->SetVertexMark(seed, cw_mark);
			while (!stack.IsEmpty())
			{
				R3MeshVertex* vertex = stack.Tail();
				stack.RemoveTail();
				for (int i=0; i<mesh->VertexValence(vertex); i++)
				{
					R3MeshEdge* edge = mesh->EdgeOnVertex(vertex, i);
					R3MeshVertex* neighbor = mesh->VertexAcrossEdge(edge, vertex);
					if (!neighbor) continue;
					if (mesh->VertexMark(neighbor) == cw_mark)
						continue;
					if (mesh->VertexMark(neighbor) == ccw_mark)
					{
						cout<<"cw to ccw conflict!"<<endl;
						cout<<idx<<endl;
						if (print_debug)
							getchar();
					}
					// set new mark
					if (mesh->VertexMark(neighbor) == init_mark)
					{
						stack.Insert(neighbor);
						mesh->SetVertexMark(neighbor, cw_mark);
					}
				}
			}
		}
		idx ++;
	}
    
	for (int i=0; i<mesh->NVertices(); i++)
	{
		R3MeshVertex* vertex = mesh->Vertex(i);
		if (mesh->VertexMark(vertex) == ccw_mark)
		{
			sides[i] = 0;
		}
		else if (mesh->VertexMark(vertex) == cw_mark)
		{
			sides[i] = 1;
		}
		else if (mesh->VertexMark(vertex) == axis_mark)
		{
			sides[i] = 2;
		}
		else if (mesh->VertexMark(vertex) == init_mark)
		{
            sides[i] = 2;
//			cout<<"still init_mark!"<<endl;
//			cout<<i<<endl;
//			getchar();
		}
		else
		{
			cout<<"Should not be such a mark!"<<endl;
			cout<<i<<endl;
			getchar();
		}
	}
	return sides;
}

static int*
CheckConnectivity(R3Mesh* mesh, int* sides, int flag)
{
	int * array = new int[mesh->NVertices()];
	for (int i=0; i<mesh->NVertices(); i++)
		array[i] = 0;
	for (int i=0; i<mesh->NFaces(); i++)
	{
		R3MeshFace* f = mesh->Face(i);
		// Check whether the face is flagged
		R3MeshVertex* v[3];
		v[0] = mesh->VertexOnFace(f, 0);
		v[1] = mesh->VertexOnFace(f, 1);
		v[2] = mesh->VertexOnFace(f, 2);
		if (!v[0] || !v[1] || !v[2])
			continue;
		int idx[3];
		idx[0] = mesh->VertexID(v[0]);
		idx[1] = mesh->VertexID(v[1]);
		idx[2] = mesh->VertexID(v[2]);
		if (!flag)
		{
			if (sides[idx[0]] == 1 ||
                sides[idx[1]] == 1 ||
                sides[idx[2]] == 1)
				continue;
		}
		else
		{
			if (sides[idx[0]] == 0 ||
                sides[idx[1]] == 0 ||
                sides[idx[2]] == 0 )
				continue;
		}
        
		array[idx[0]] = 1;
		array[idx[1]] = 1;
		array[idx[2]] = 1;
	}
	return array;
}


static R3Mesh*
CreateMeshOnOneSide(R3Mesh* mesh, int* sides, int flag,
                    DenseVertexCorrespondence** mapping)
{
	R3Mesh* _mesh = new R3Mesh;
	*mapping = new DenseVertexCorrespondence(mesh, _mesh);
	// Return an array with length = mesh->NVertices()
	int* connectivity_flag = CheckConnectivity(mesh, sides, flag);
	// Vertices
	int counter = 0;
	for (int i=0; i<mesh->NVertices(); i++)
	{
		R3MeshVertex* v = mesh->Vertex(i);
		int idx = mesh->VertexID(v);
		if (flag == 0 && sides[idx] == 1)
			continue;
		if (flag == 1 && sides[idx] == 0)
			continue;
		_mesh->CreateVertex(mesh->VertexPosition(v));
		(*mapping)->vertices[i] = _mesh->Vertex(counter);
		counter ++;
	}
	// Faces
	for (int i=0; i<mesh->NFaces(); i++)
	{
		R3MeshFace* f = mesh->Face(i);
		R3MeshVertex* v[3];
		v[0] = mesh->VertexOnFace(f, 0);
		v[1] = mesh->VertexOnFace(f, 1);
		v[2] = mesh->VertexOnFace(f, 2);
		if (!v[0] || !v[1] || !v[2])
			continue;
		int idx[3];
		idx[0] = mesh->VertexID(v[0]);
		idx[1] = mesh->VertexID(v[1]);
		idx[2] = mesh->VertexID(v[2]);
		if (!flag)
		{
			if (sides[idx[0]] == 1 ||
                sides[idx[1]] == 1 ||
                sides[idx[2]] == 1 )
				continue;
		}
		else
		{
			if (sides[idx[0]] == 0 ||
                sides[idx[1]] == 0 ||
                sides[idx[2]] == 0 )
				continue;
		}
		_mesh->CreateFace((*mapping)->vertices[idx[0]], (*mapping)->vertices[idx[1]],
                          (*mapping)->vertices[idx[2]]);
	}
	delete[] connectivity_flag;
	return _mesh;
}

// forward mapping -> backward mapping
static DenseVertexCorrespondence*
CreateInverseMapping(DenseVertexCorrespondence* mapping)
{
	R3Mesh* mesh0 = mapping->mesh[1];
	R3Mesh* mesh1 = mapping->mesh[0];
	DenseVertexCorrespondence* inv_mapping = new DenseVertexCorrespondence(mesh0, mesh1);
	for (int i=0; i<mesh1->NVertices(); i++)
	{
		R3MeshVertex* v = mapping->vertices[i];
		if (v != NULL)
		{
			int idx = mesh0->VertexID(v);
			inv_mapping->vertices[idx] = mesh1->Vertex(i);
		}
	}
	return inv_mapping;
}


// return side function
static int*
CutMeshAlongAxes(R3Mesh* whole, RNArray<R3MeshVertex*>* axis,
                 R3Mesh** half0, R3Mesh** half1,
                 DenseVertexCorrespondence** half2whole_map0,
                 DenseVertexCorrespondence** half2whole_map1)
{
    int* side_function = CreateSideFunction(whole, axis);
    // forward maps: whole -> halves
    DenseVertexCorrespondence* forward_maps[2];
    *half0 = CreateMeshOnOneSide(whole, side_function, 0, &forward_maps[0]);
    *half1 = CreateMeshOnOneSide(whole, side_function, 1, &forward_maps[1]);
    // forward maps -> backward maps
    *half2whole_map0 = CreateInverseMapping(forward_maps[0]);
    *half2whole_map1 = CreateInverseMapping(forward_maps[1]);
    delete forward_maps[0];
    delete forward_maps[1];
    return side_function;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	R3Mesh* whole_mesh = ReadMesh(g_input_mesh_name);
    RNArray<R3MeshVertex*>* axis = ReadPoints(whole_mesh, g_input_axis_name);
    // step 1: cut mesh into two halves along axis
    R3Mesh* half_meshes[2] = {NULL, NULL};
    DenseVertexCorrespondence* half2whole_maps[2] = {NULL, NULL};
    int* side_function = CutMeshAlongAxes(whole_mesh, axis,
                                          &half_meshes[0], &half_meshes[1],
                                          &half2whole_maps[0],
                                          &half2whole_maps[1]);
    // step 2: search tree for two half-meshes
    R3MeshSearchTree* search_trees[2];
    search_trees[0] = new R3MeshSearchTree(half_meshes[0]);
    search_trees[1] = new R3MeshSearchTree(half_meshes[1]);
	if (!search_trees[0] || !search_trees[1])
	{
		fprintf(stderr, "Unable to create search tree!\n");
		return 0;
	}
    // step 3: refining map
    DenseVertexCorrespondence* unrefined_map = ReadDenseVertexCorrespondence(whole_mesh, whole_mesh, g_input_map_name);
    DenseVertexCorrespondence* refined_map = new DenseVertexCorrespondence(whole_mesh, whole_mesh);
    RNTime start_time;
	start_time.Read();
    for (int i=0; i<whole_mesh->NVertices(); i++)
    {
//        printf ("Vertex ID = %d side_function[%d] = %d\n", i, i, side_function[i]);
//        fflush(stdout);
        // flip the flag
        int desired_flag = 0;
        if (side_function[i] != 2)
        {
            desired_flag = 1 - side_function[i];
        }
        // unrefned mapped vertex
        R3MeshVertex* unrefined_mapped_vertex = unrefined_map->vertices[i];
        R3MeshVertex* refined_mapped_vertex = NearestNeighbor(search_trees[desired_flag], whole_mesh->VertexPosition(unrefined_mapped_vertex));
        // mapped to whole mesh vertex
        int refined_mapped_id = half_meshes[desired_flag]->VertexID(refined_mapped_vertex);
        refined_mapped_vertex = half2whole_maps[desired_flag]->vertices[refined_mapped_id];
        refined_map->vertices[i] = refined_mapped_vertex;
    }
	cout<<"***** Time for refinement : "<<start_time.Elapsed()<<" seconds"<<endl;
    WriteDenseVertexCorrespondence(refined_map, g_output_map_name);
	return 0;
}
