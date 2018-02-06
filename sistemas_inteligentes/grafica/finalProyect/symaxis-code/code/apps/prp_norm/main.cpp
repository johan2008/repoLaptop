#include "R3Shapes/R3Shapes.h"
#include "RNExt/RNExt.h"
////////////////////////////////////////////////////////////////////////
// Program     : prp_norm
// Usage       : prp_norm mesh in_prp out_prp
// Description : Normalize properties
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name = NULL;
char * g_input_prp_name = NULL;
char * g_output_prp_name = NULL;

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * g_mesh = NULL;
R3MeshPropertySet* g_properties = NULL;
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
        	if (!g_input_mesh_name) g_input_mesh_name = *argv;
			else if (!g_input_prp_name) g_input_prp_name = *argv;
			else if (!g_output_prp_name) g_output_prp_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_mesh_name || !g_input_prp_name || !g_output_prp_name) {
    	fprintf(stderr, "Usage: prp_norm meshname inprp outprp\n");
    	return 0;
  	}

	return 1;
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

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	g_mesh = ReadMesh(g_input_mesh_name);
	g_properties = ReadProperties(g_mesh, g_input_prp_name);
	for (int i=0; i<g_properties->NProperties(); i++)
	{
		NormalizeProperty(g_properties->Property(i));
	}
	g_properties->Write(g_output_prp_name);
	return 0;
}
