/* read mesh from file and 
 - if one vertex is specified, for all vertices of the mesh print their distances to this vertex
 - if two vertices are specified, print the shortest path between these vertices 

	Danil Kirsanov, 01/2008 
*/

/*
 * Changed to read from off file, write to text file, one distance per line
 *
 * Raif
*/

#include <iostream>
#include <fstream>

#include "geodesic_algorithm_exact.h"


int main(int argc, char **argv) 
{
	if(argc < 3)
	{
		std::cout << "usage: infile source_vertex outfile" << std::endl; //try: "hedgehog_mesh.txt 3 14" or "flat_triangular_mesh.txt 1"
		return 0;
	}

	std::vector<double> points;	
	std::vector<unsigned> faces;

	bool success = geodesic::off_read_mesh_from_file(argv[1],points,faces);
	if(!success)
	{
		std::cout << "something is wrong with the input file" << std::endl;
		return 0;
	}

	geodesic::Mesh mesh;
	mesh.initialize_mesh_data(points, faces);		//create internal mesh data structure including edges

	geodesic::GeodesicAlgorithmExact algorithm(&mesh);	//create exact algorithm for the mesh

	unsigned source_vertex_index = atol(argv[2]);

	geodesic::SurfacePoint source(&mesh.vertices()[source_vertex_index]);		//create source 
	std::vector<geodesic::SurfacePoint> all_sources(1,source);					//in general, there could be multiple sources, but now we have only one

	
	algorithm.propagate(all_sources);	//cover the whole mesh

	std::ofstream myfile;
	myfile.open (argv[3]);


	for(unsigned i=0; i<mesh.vertices().size(); ++i)
	{
			geodesic::SurfacePoint p(&mesh.vertices()[i]);		

			double distance;
			unsigned best_source = algorithm.best_source(p,distance);		//for a given surface point, find closets source and distance to this source

			myfile << distance << std::endl;		//print geodesic distance for every vertex
	}
	
	myfile.close();	
	return 0;
}	
