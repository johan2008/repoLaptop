#include "R3Shapes/R3Shapes.h"
#include "RNExt/RNExt.h"

////////////////////////////////////////////////////////////////////////
// Program     : myprp2prp
// Usage	   : myprp2prp mesh inprp outprp [-add prp weights, -clamp min max]
// Description : combine prps
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_mesh_name = NULL;
char * g_input_prp_name = NULL;
char * g_output_prp_name = NULL;
vector<char*> g_extra_input_prp_name;
vector<RNScalar> g_weights;

////////////////////////////////////////////////////////////////////////
// Data structure
////////////////////////////////////////////////////////////////////////

R3Mesh * g_mesh = NULL;
R3MeshPropertySet* g_properties = NULL;

RNScalar g_min = -FLT_MAX;
RNScalar g_max = FLT_MAX;
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
			else if (!strcmp(*argv, "-add"))
			{
				argc--; argv++;
				g_extra_input_prp_name.push_back(*argv);
				argc--; argv++;
				g_weights.push_back(atof(*argv));
			}
			else if (!strcmp(*argv, "-clamp"))
			{
				argc--; argv++;
				g_min = atof(*argv);
				argc--; argv++;
				g_max = atof(*argv);
			}
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
    	fprintf(stderr, "Usage: myprp2prp mesh inprp outprp [-add prp weights, -clamp min max]\n");
    	return 0;
  	}

	return 1;
}

int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	g_mesh = ReadMesh(g_input_mesh_name);
	g_properties = ReadProperties(g_mesh, g_input_prp_name);
	// -clamp
	for (int i=0; i<g_properties->NProperties(); i++)
	{
		R3MeshProperty* _p = g_properties->Property(i);
		for (int j=0; j<g_mesh->NVertices(); j++)
		{
			RNScalar val = _p->VertexValue(j);
			if (val > g_max) val = g_max;
			if (val < g_min) val = g_min;
			_p->SetVertexValue(j, val);
		}
	}
	// -add
	for (int i=0; i<int(g_extra_input_prp_name.size()); i++)
	{
		char * filename = g_extra_input_prp_name[i];
		R3MeshPropertySet* prp = ReadProperties(g_mesh, filename);
		printf("Adding a new set of properties ... (%d)\n", prp->NProperties());
		for (int k=0; k<prp->NProperties(); k++)
		{
			R3MeshProperty* _p = prp->Property(k);
			_p->Multiply(g_weights[i]);
			g_properties->Insert(_p);
		}
	}

	printf("Finally we get %d properties.\n", g_properties->NProperties());
	
	g_properties->Write(g_output_prp_name);	
	return 0;
}
