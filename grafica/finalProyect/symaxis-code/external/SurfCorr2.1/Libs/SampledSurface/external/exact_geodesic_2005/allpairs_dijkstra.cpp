/* different algorithms; multiple sources and targets
	Danil Kirsanov, 01/2008 
*/
#include <iostream>
#include <fstream>
#include <string>
#include <ctime> 

#include "geodesic_algorithm_dijkstra.h"
#include "geodesic_algorithm_subdivision.h"
#include "geodesic_algorithm_exact.h"

int main(int argc, char **argv) 
{
	if(argc < 1)
	{
		std::cout << "usage: infile" << std::endl; //try: "hedgehog_mesh.txt 3 14" or "flat_triangular_mesh.txt 1"
		std::cout << "Timing Dijkstra algorithm provided by Kirsanov. All pairwise distances are intended. Only propagation time is shown." << std::endl;
		return 0;
	}
	
	std::vector<double> points;	
	std::vector<unsigned> faces;
	geodesic::off_read_mesh_from_file(argv[1],points,faces); 

	geodesic::Mesh mesh;
	mesh.initialize_mesh_data(points, faces);		//create internal mesh data structure including edges

	geodesic::GeodesicAlgorithmExact exact_algorithm(&mesh);		//exact algorithm
	geodesic::GeodesicAlgorithmDijkstra dijkstra_algorithm(&mesh);		//simplest approximate algorithm: path only allowed on the edges of the mesh
	unsigned const subdivision_level = 3;										//three additional vertices per every edge in subdivision algorithm
	geodesic::GeodesicAlgorithmSubdivision subdivision_algorithm(&mesh,2);	//with subdivision_level=0 this algorithm becomes Dijkstra, with subdivision_level->infinity it becomes exact

	std::vector<geodesic::GeodesicAlgorithmBase* > all_algorithms;		//for simplicity, store all possible geodesic algorithms here
	all_algorithms.push_back(&dijkstra_algorithm);
	//all_algorithms.push_back(&subdivision_algorithm);
	//all_algorithms.push_back(&exact_algorithm);

	std::vector<geodesic::SurfacePoint> sources;
	
	//sources.push_back(geodesic::SurfacePoint(&mesh.edges()[12]));		//second source is located in the middle of edge 12
	//sources.push_back(geodesic::SurfacePoint(&mesh.faces()[20]));		//third source is located in the middle of face 20

	std::vector<geodesic::SurfacePoint> targets;		//same thing with targets
	targets.push_back(geodesic::SurfacePoint(&mesh.vertices().back()));		
	targets.push_back(geodesic::SurfacePoint(&mesh.edges()[10]));		
	targets.push_back(geodesic::SurfacePoint(&mesh.faces()[3]));		

	for(unsigned index=0; index<all_algorithms.size(); ++index)
	{
		geodesic::GeodesicAlgorithmBase* algorithm = all_algorithms[index];		//all algorithms are derived from GeodesicAlgorithmBase
		std::cout << std::endl << "results for algorithm " << algorithm->name() << std::endl;
	
		clock_t start = clock();
		
		for(int a=0; a<mesh.vertices().size();a++){
			sources.push_back(geodesic::SurfacePoint(&mesh.vertices()[a]));		//one source is located at vertex zero
			algorithm->propagate(sources);		//cover the whole mesh
			sources.pop_back();
		}
	
		clock_t stop = clock();
		double m_time_consumed = (static_cast<double>(stop)-static_cast<double>(start))/CLOCKS_PER_SEC;
		std::cout << std::endl << "time consumed"<<  m_time_consumed << std::endl;
		
	
	}
	return 0;
}	
