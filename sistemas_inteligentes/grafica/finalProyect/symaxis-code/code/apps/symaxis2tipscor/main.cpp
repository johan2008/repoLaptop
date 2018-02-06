#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////
// Program     : symaxis2tipscor
// Usage       : symaxis2tipscor mesh tips axis symcorr  [-v, -debug] 
// Description : Generate symmetric pairs in tips from symmetry axis only (mutually closest)
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * input_mesh_name = NULL;
char * input_tips_name = NULL;
char * input_axis_name = NULL;
char * output_symcorr_name = NULL;

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * Mesh = NULL;
RNArray<R3MeshVertex*>* Tips = NULL;
RNArray<R3MeshVertex*>* Axis = NULL;
// Sides on the surface based on the symmetry axis
int * Sides = NULL;
struct Elem
{
	// Flag: on if it is a pair
	int isPair;
	R3MeshVertex* vertex[2];
};
RNArray<Elem*>* PointsRecorder = NULL;

////////////////////////////////////////////////////////////////////////
// Flags
////////////////////////////////////////////////////////////////////////

int print_verbose = 0;
int print_debug = 0;

////////////////////////////////////////////////////////////////////////
// Variables
////////////////////////////////////////////////////////////////////////
// Threshold to tell a point if it is on the axis
RNScalar thres = 0.1;

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

// Compute sides based on the axis
// Val = 1: ccw (or right on the axis)
// Val = -1: cw
// Create Sides on the mesh, from axis
// Left mesh: sides = 1
// Right mesh: sides = -1
// On the axis: sides = 0
static int*
CreateSides(R3Mesh* mesh, RNArray<R3MeshVertex*>* axis)
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
//		cout<<visited_counter<<" vertices visited during ccw."<<endl;
//		exit(-1);
//		getchar();
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
			sides[i] = 1;
		}
		else if (mesh->VertexMark(vertex) == cw_mark)
		{
			sides[i] = -1;
		}
		else if (mesh->VertexMark(vertex) == axis_mark)
		{
			sides[i] = 0;
		}
		else if (mesh->VertexMark(vertex) == init_mark)
		{
			cout<<"still init_mark!"<<endl;
			cout<<i<<endl;
			getchar();
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

/*static int**/
/*CreateSides(R3Mesh* mesh, RNArray<R3MeshVertex*>* axis)*/
/*{*/
	/*int * sides = new int [mesh->NVertices()];*/
	/*for (int i=0; i<mesh->NVertices(); i++) sides[i] = -1;*/

	/*RNArray<R3MeshEdge*> axis_edges;*/
	/*for (int i=0; i<axis->NEntries(); i++)*/
	/*{*/
		/*R3MeshVertex* v[2];*/
		/*v[0] = axis->Kth(i);*/
		/*v[1] = axis->Kth((i+1)%axis->NEntries());*/

		/*R3MeshEdge* edge = mesh->EdgeBetweenVertices(v[0], v[1]);*/
		/*if (edge)*/
			/*axis_edges.Insert(edge);*/
		/*else */
		/*{*/
			/*cout<<"Axis is disconnected!"<<endl;*/
			/*exit(-1);*/
		/*}*/
	/*}*/

	/*// Prepare : set all vertices to be intial*/
	/*R3mesh_mark ++;*/
	/*int init_mark = R3mesh_mark;*/
	/*for (int i=0; i<mesh->NVertices(); i++)*/
		/*mesh->SetVertexMark(mesh->Vertex(i), init_mark);*/
	/*// Preparation: Mark the vertices on axis*/
	/*R3mesh_mark ++;*/
	/*int axis_mark = R3mesh_mark;*/
	/*for (int i=0; i<axis->NEntries(); i++)*/
	/*{*/
		/*mesh->SetVertexMark(axis->Kth(i), axis_mark);*/
	/*}*/
	/*R3mesh_mark ++;*/
	/*int cw_mark = R3mesh_mark;*/
	/*R3mesh_mark ++;*/
	/*int ccw_mark = R3mesh_mark;*/
	/*int idx = 0;*/
	/*while (idx < axis->NEntries())*/
	/*{*/
		/*R3MeshVertex* seed = NULL;*/
		/*R3MeshVertex* vertex = axis->Kth(idx);*/
		/*R3MeshEdge* edge = axis_edges.Kth(idx);*/
		/*edge = mesh->EdgeOnVertex(vertex, edge, RN_CCW);*/
		/*if (!edge) {idx++;continue;}*/
		/*seed = mesh->VertexAcrossEdge(edge, vertex);*/
		/*if (!seed) {idx++;continue;}*/
		/*if (mesh->VertexMark(seed) == axis_mark) // another axis edge*/
		/*{ idx++;continue; }*/
		/*if (mesh->VertexMark(seed) == ccw_mark) // alread marked*/
		/*{ idx++;continue;}*/
		/*else if (mesh->VertexMark(seed) == cw_mark)*/
		/*{*/
			/*cout<<"ccw to cw conflict (init)!"<<endl;*/
			/*cout<<idx<<endl;*/
			/*if (print_debug)*/
				/*getchar();*/
			/*idx++;*/
			/*continue;*/
		/*}*/
		/*RNArray<R3MeshVertex*> stack;*/
		/*stack.Insert(seed);*/
		/*mesh->SetVertexMark(seed, ccw_mark);*/
		/*while (!stack.IsEmpty())*/
		/*{*/
			/*R3MeshVertex* vertex = stack.Tail();*/
			/*stack.RemoveTail();*/
			/*for (int i=0; i<mesh->VertexValence(vertex); i++)*/
			/*{*/
				/*R3MeshEdge* edge = mesh->EdgeOnVertex(vertex, i);*/
				/*R3MeshVertex* neighbor = mesh->VertexAcrossEdge(edge, vertex);*/
				/*if (!neighbor) continue;*/
				/*if (mesh->VertexMark(neighbor) == ccw_mark)*/
					/*continue;*/
				/*if (mesh->VertexMark(neighbor) == cw_mark)*/
				/*{*/
					/*cout<<"ccw to cw conflict!"<<endl;*/
					/*cout<<idx<<endl;*/
					/*if (print_debug)*/
						/*getchar();*/
				/*}*/
				/*// set new mark*/
				/*if (mesh->VertexMark(neighbor) == init_mark)*/
				/*{*/
					/*stack.Insert(neighbor);*/
					/*mesh->SetVertexMark(neighbor, ccw_mark);*/
				/*}*/
			/*}*/
		/*}*/
		/*// cw*/
		/*vertex = axis->Kth(idx);*/
		/*edge = axis_edges.Kth(idx);*/
		/*edge = mesh->EdgeOnVertex(vertex, edge, RN_CW);*/
		/*if (!edge) {idx++;continue;}*/
		/*seed = mesh->VertexAcrossEdge(edge, vertex);*/
		/*if (!seed) {idx++;continue;}*/
		/*if (mesh->VertexMark(seed) == axis_mark) // another axis edge*/
		/*{idx++;continue;}*/
		/*if (mesh->VertexMark(seed) == cw_mark) // alread marked*/
		/*{idx++;continue;}*/
		/*else if (mesh->VertexMark(seed) == ccw_mark)*/
		/*{*/
			/*cout<<"cw to ccw conflict (init)!"<<endl;*/
			/*cout<<idx<<endl;*/
			/*if (print_debug)*/
				/*getchar();*/
			/*idx++;*/
			/*continue;*/
		/*}*/
		/*stack.Empty();*/
		/*stack.Insert(seed);*/
		/*mesh->SetVertexMark(seed, cw_mark);*/
		/*while (!stack.IsEmpty())*/
		/*{*/
			/*R3MeshVertex* vertex = stack.Tail();*/
			/*stack.RemoveTail();*/
			/*for (int i=0; i<mesh->VertexValence(vertex); i++)*/
			/*{*/
				/*R3MeshEdge* edge = mesh->EdgeOnVertex(vertex, i);*/
				/*R3MeshVertex* neighbor = mesh->VertexAcrossEdge(edge, vertex);*/
				/*if (!neighbor) continue;*/
				/*if (mesh->VertexMark(neighbor) == cw_mark)*/
					/*continue;*/
				/*if (mesh->VertexMark(neighbor) == ccw_mark)*/
				/*{*/
					/*cout<<"cw to ccw conflict!"<<endl;*/
					/*cout<<idx<<endl;*/
					/*if (print_debug)*/
						/*getchar();*/
				/*}*/
				/*// set new mark*/
				/*if (mesh->VertexMark(neighbor) == init_mark)*/
				/*{*/
					/*stack.Insert(neighbor);*/
					/*mesh->SetVertexMark(neighbor, cw_mark);*/
				/*}*/
			/*}*/
		/*}*/
		/*idx ++;*/
	/*}	*/

	/*for (int i=0; i<mesh->NVertices(); i++)*/
	/*{*/
		/*R3MeshVertex* vertex = mesh->Vertex(i);*/
		/*if (mesh->VertexMark(vertex) == ccw_mark)*/
		/*{*/
			/*sides[i] = 1;*/
		/*}*/
		/*else if (mesh->VertexMark(vertex) == cw_mark)*/
		/*{*/
			/*sides[i] = -1;*/
		/*}*/
		/*else if (mesh->VertexMark(vertex) == axis_mark)*/
		/*{*/
			/*sides[i] = 1;*/
		/*}*/
		/*else if (mesh->VertexMark(vertex) == init_mark)*/
		/*{*/
			/*cout<<"still init_mark!"<<endl;*/
			/*sides[i] = 1;*/
/*//			cout<<i<<endl;*/
/*//			getchar();*/
		/*}*/
		/*else*/
		/*{*/
			/*cout<<"Should not be such a mark!"<<endl;*/
			/*cout<<i<<endl;*/
			/*getchar();*/
		/*}*/
	/*}*/
	/*return sides;*/
/*}*/

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
			else if (!strcmp(*argv, "-thres"))
			{
				argc--; argv++;
				thres = atof(*argv);
			}
        	else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
    	else {
        	if (!input_mesh_name) input_mesh_name = *argv;
			else if (!input_tips_name) input_tips_name = *argv;
			else if (!input_axis_name) input_axis_name = *argv;
			else if (!output_symcorr_name) output_symcorr_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!input_mesh_name || !input_tips_name || !input_axis_name
			|| !output_symcorr_name) {
    	fprintf(stderr, "Usage: symaxis2tipscor mesh tips axis symcorr  [-v, -debug, -thres]\n");
    	return 0;
  	}

	return 1;
}

static RNScalar
DiffVector(RNScalar* d0, RNScalar* d1, int dim)
{
	RNScalar sum = 0;
	for (int i=0; i<dim; i++)
	{
		RNScalar d = fabs(d0[i] - d1[i]);
		if (d0[i] > 0.0)
			sum += d/d0[i];
		else if (d > 0.0)
			sum += FLT_MAX;
	}
	return sum;
}

static RNArray<Elem*>*
FindCorrespondences(R3Mesh* mesh, int* sides, RNArray<R3MeshVertex*>* axis, RNArray<R3MeshVertex*>* tips)
{
	RNArray<Elem*>* elems = new RNArray<Elem*>;
	// build up the feature descriptors: distances to the axis
	RNScalar** tips_distances = new RNScalar*[tips->NEntries()];
	for (int i=0; i<tips->NEntries(); i++)
	{
		tips_distances[i] = new RNScalar[axis->NEntries()];
		RNLength* d = mesh->DijkstraDistances(tips->Kth(i));
		for (int j=0; j<axis->NEntries(); j++)
		{
			tips_distances[i][j] = d[mesh->VertexID(axis->Kth(j))];
		}
		delete[] d;
	}

	// find mutually closest points as for the distance between tips_distances
	int ntips = tips->NEntries();
	int naxis = axis->NEntries();
	int* tips_flag = new int[ntips];
	for (int i=0; i<ntips; i++)
	{
		tips_flag[i] = 0;
	}

	for (int i=0; i<ntips; i++)
	{
		RNScalar* d0 = tips_distances[i];
		// only process the tips with sides = -1
		if (sides[mesh->VertexID(tips->Kth(i))] == 1)
			continue;
		// find the best d1
		int best_idx = -1;
		RNScalar best_val = FLT_MAX;
		for (int j=0; j<ntips; j++)
		{
			if (sides[mesh->VertexID(tips->Kth(j))] == -1)
				continue;
			RNScalar* d1 = tips_distances[j];
			RNScalar diff = DiffVector(d0, d1, naxis);
			if (diff < best_val)
			{
				best_val = diff;
				best_idx = j;
			}
		}
		// if no one left
		if (best_idx == -1)
		{
			continue;
		}
		// if the one has been taken
		if (tips_flag[best_idx])
		{
			continue;
		}
		// find the best d2
		int best_idx2 = -1;
		RNScalar best_val2 = FLT_MAX;
		RNScalar* d1 = tips_distances[best_idx];
		for (int j=0; j<ntips; j++)
		{
			if (sides[mesh->VertexID(tips->Kth(j))] == 1)
				continue;
			RNScalar* d2 = tips_distances[j];
			RNScalar diff = DiffVector(d1, d2, naxis);
			if (diff < best_val2)
			{
				best_val2 = diff;
				best_idx2 = j;
			}
		}
		// mutually closest
		if (best_idx2 == i)
		{
			tips_flag[i] = 1;
			tips_flag[best_idx] = 1;
			Elem* elem = new Elem;
			elem->isPair = 1;
			elem->vertex[0] = tips->Kth(i);
			elem->vertex[1] = tips->Kth(best_idx);
			elems->Insert(elem);
		}
	}
	// revisit every one not visited
	RNLength* dist2axis = mesh->DijkstraDistances(*axis);
	for (int i=0; i<ntips; i++)
	{
		if (tips_flag[i] == 0 && dist2axis[mesh->VertexID(tips->Kth(i))]<thres)
		{
			Elem* elem = new Elem;
			elem->isPair = 0;
			elem->vertex[0] = tips->Kth(i);
			elem->vertex[1] = tips->Kth(i);
			elems->Insert(elem);
		}
	}
	delete[] dist2axis;
	delete[] tips_flag;
	for (int i=0; i<tips->NEntries(); i++)
	{
		delete[] tips_distances[i];
	}
	delete[] tips_distances;
	return elems;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	// Mesh
	Mesh = ReadMesh(input_mesh_name);
	// Tips
	Tips = ReadPoints(Mesh, input_tips_name);
	// Symmetry axis
	Axis = ReadPoints(Mesh, input_axis_name);	
	// Compute sides
	Sides = CreateSides(Mesh, Axis);
	if (print_debug)
	{
		R3MeshProperty * side_prp = new R3MeshProperty(Mesh);
		for (int i=0; i<Mesh->NVertices(); i++)
		{
			if (Sides[i] == 1)
				side_prp->SetVertexValue(i, 1.0);
			else if (Sides[i] == -1)
				side_prp->SetVertexValue(i, -1.0);
		}
		side_prp->Write("sides.val");
	}

	// Search for pairs and singletons
	PointsRecorder = FindCorrespondences(Mesh, Sides, Axis, Tips);

	// Write to file
	ofstream fout(output_symcorr_name);
	for (int i=0; i<PointsRecorder->NEntries(); i++)
	{
		Elem* elem = PointsRecorder->Kth(i);
		if (elem->isPair)
			fout<<Mesh->VertexID(elem->vertex[0])<<" "<<Mesh->VertexID(elem->vertex[1])<<endl;
	}
	for (int i=0; i<PointsRecorder->NEntries(); i++)
	{
		Elem* elem = PointsRecorder->Kth(i);
		if (!elem->isPair)
			fout<<Mesh->VertexID(elem->vertex[0])<<" "<<Mesh->VertexID(elem->vertex[1])<<endl;
	}
	fout.close();
	return 0;
}
