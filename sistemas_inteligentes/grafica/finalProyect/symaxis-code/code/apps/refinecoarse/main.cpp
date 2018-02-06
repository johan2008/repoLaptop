#include "R3Shapes/R3Shapes.h"
#include "RNExt/RNExt.h"
#include <tr1/unordered_map>

////////////////////////////////////////////////////////////////////////
// Program     : refinecoarse
// Usage       : refinecoarse coarse_init coarse_final
// Description : Refine a set of coarse_correspondences (not repetition)
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Input file name
////////////////////////////////////////////////////////////////////////

char * g_input_cor_name = NULL;
char * g_output_cor_name = NULL;

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
			if (!g_input_cor_name) g_input_cor_name = *argv;
			else if (!g_output_cor_name) g_output_cor_name = *argv;
			else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
        	argv++; argc--;
    	}
  	}
	// Check input filename
  	if (!g_input_cor_name || !g_output_cor_name) {
    	fprintf(stderr, "Usage       : refinecoarse coarse_init coarse_final\n");
    	return 0;
  	}

	return 1;
}

static void
ReadCor(char* filename, vector<int> c[])
{
	ifstream fin(filename);
	if (!fin)
	{
		cout<<"Unable to open "<<filename<<endl;
		exit(-1);
	}
	int temp[2];
	while (fin >> temp[0] >> temp[1])
	{
		c[0].push_back(temp[0]);
		c[1].push_back(temp[1]);
	}
	fin.close();
}


int main(int argc, char ** argv)
{
	if (!ParseArgs(argc, argv)) exit(1);
	vector<int> input_cor[2];
	ReadCor(g_input_cor_name, input_cor);

	vector<int> output_cor[2];
	tr1::unordered_map<int, int> visited[2];
	for (int i=0; i<input_cor[0].size(); i++)
	{
		int a = input_cor[0][i];
		int b = input_cor[1][i];
		if (visited[0].find(a) == visited[0].end() && visited[1].find(b) == visited[1].end())
		{
			output_cor[0].push_back(a);
			output_cor[1].push_back(b);
			visited[0][a] = 1;
			visited[1][b] = 1;
		}
	}
	ofstream fout(g_output_cor_name);
	for (int i=0; i<output_cor[0].size(); i++)
	{
		fout<<output_cor[0][i]<<" "<<output_cor[1][i]<<endl;
	}
	fout.close();
	return 0;
}
