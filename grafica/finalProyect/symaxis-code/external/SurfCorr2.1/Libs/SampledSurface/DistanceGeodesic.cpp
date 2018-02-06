#include "DistanceGeodesic.h"
#include "SampledSurface.h"

#include <math.h>
#include <iostream>
#include <memory>

#include "./external/exact_geodesic_2005/geodesic_algorithm_exact.h"
#include "./external/exact_geodesic_2005/geodesic_algorithm_base.h"
#include "./external/exact_geodesic_2005/geodesic_algorithm_dijkstra_alternative.h"
#include "./external/exact_geodesic_2005/geodesic_algorithm_dijkstra.h"
#include "./external/exact_geodesic_2005/geodesic_algorithm_graph_base.h"
#include "./external/exact_geodesic_2005/geodesic_algorithm_subdivision.h"
#include "./external/exact_geodesic_2005/geodesic_constants_and_simple_functions.h"
#include "./external/exact_geodesic_2005/geodesic_memory.h"
#include "./external/exact_geodesic_2005/geodesic_mesh_elements.h"
#include "./external/exact_geodesic_2005/geodesic_mesh.h"

DistanceGeodesic::DistanceGeodesic(SampledSurface * surface, CalculationType geodesicType)
: SurfaceDistance(surface, true)
{
	m_GeodesicAlgorithm = geodesicType;
	m_bVerbous = true;
	m_iPrecompWidth = surface->GetMesh()->NVertices();
	m_iPrecompHeight = surface->GetMesh()->NVertices();	
}

void DistanceGeodesic::PrecomputeDistances()
{
	m_aPrecomputedDistances = new double[m_iPrecompWidth * m_iPrecompHeight];
	
	switch(m_GeodesicAlgorithm)
	{
		case GEODIST_DIJKSTRA_FUNKHOUSER:
			PrecomputeFunkhouser();
			break;
		case GEODIST_DIJKSTRA_KIRSANOV1:
		case GEODIST_DIJKSTRA_KIRSANOV2:
		case GEODIST_EXACT_KIRSANOV:
		case GEODIST_SUBDIVISION_KIRSANOV:
			PrecomputeKirsanov();
			break;
		default:
			assert(false);
	}
}

void DistanceGeodesic::FillThresholdedDistances(GeoVertexData * vertex_data, R3Mesh * mesh, int fromVertex,
											  std::map<int, double> & verticesWithValues, double maxDistance)
{
	R3MeshVertex *source_vertex = mesh->Vertex(fromVertex);

	// distance to self
	verticesWithValues[fromVertex] = 0;
	
	// Re-set vertex data
	for (int j = 0; j < mesh->NVertices(); j++) 
	{
		vertex_data[j].distance = FLT_MAX;
		vertex_data[j].heappointer = NULL;
	}
	
	// Initialize priority queue
	GeoVertexData *data = (GeoVertexData *) mesh->VertexData(source_vertex);
	assert(data!=NULL);
	data->distance = 0;
	RNHeap<GeoVertexData *> heap(data, &(data->distance), &(data->heappointer));
	heap.Push(data);		
	
	// Visit other nodes computing shortest distance
	while (!heap.IsEmpty()) 
	{
		GeoVertexData *data = heap.Pop();
		R3MeshVertex *vertex = data->vertex;
		
		for (int j = 0; j < mesh->VertexValence(vertex); j++) 
		{
			R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, j);
			R3MeshVertex *neighbor_vertex = mesh->VertexAcrossEdge(edge, vertex);
			GeoVertexData *neighbor_data = (GeoVertexData *) mesh->VertexData(neighbor_vertex);
			RNScalar old_distance = neighbor_data->distance;
			
			RNScalar new_distance = mesh->EdgeLength(edge) + data->distance;
			assert(new_distance>=0);
			if (new_distance < old_distance && new_distance < maxDistance)
			{
				neighbor_data->distance = new_distance;
				verticesWithValues[mesh->VertexID(neighbor_data->vertex)] = new_distance;
				
				if (old_distance < FLT_MAX) heap.Update(neighbor_data);
				else heap.Push(neighbor_data);
			}
		}
	}	
}

double DistanceGeodesic::GetThresholdedDistances(int toVertex, R3Mesh * mesh, 
									  std::map<int, double> & verticesWithValues, double maxDistance)
{
	std::map<int, double>::iterator iter = verticesWithValues.find(toVertex);
	if (iter==verticesWithValues.end())
		return maxDistance;
	else
		return iter->second;
}

void DistanceGeodesic::PrecomputeFillFunkhouser(GeoVertexData * vertex_data, R3Mesh * mesh,
												int fromVertex, 
												double * toDistanceArray, 
												const std::vector<int> & toVerticesRemember)
{
	R3MeshVertex *source_vertex = mesh->Vertex(fromVertex);
	
	// Re-set vertex data
	for (int j = 0; j < mesh->NVertices(); j++) 
	{
		vertex_data[j].distance = FLT_MAX;
		vertex_data[j].heappointer = NULL;
	}

	// Initialize priority queue
	GeoVertexData *data = (GeoVertexData *) mesh->VertexData(source_vertex);
	assert(data!=NULL);
	data->distance = 0;
	RNHeap<GeoVertexData *> heap(data, &(data->distance), &(data->heappointer));
	heap.Push(data);		
		
	// Visit other nodes computing shortest distance
	while (!heap.IsEmpty()) 
	{
		GeoVertexData *data = heap.Pop();
		R3MeshVertex *vertex = data->vertex;
				
		for (int j = 0; j < mesh->VertexValence(vertex); j++) 
		{
			R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, j);
			R3MeshVertex *neighbor_vertex = mesh->VertexAcrossEdge(edge, vertex);
			GeoVertexData *neighbor_data = (GeoVertexData *) mesh->VertexData(neighbor_vertex);
			RNScalar old_distance = neighbor_data->distance;
			
			RNScalar new_distance = mesh->EdgeLength(edge) + data->distance;
			assert(new_distance>=0);
			if (new_distance < old_distance) 
			{
				neighbor_data->distance = new_distance;
				if (old_distance < FLT_MAX) heap.Update(neighbor_data);
				else heap.Push(neighbor_data);
			}
		}
	}
		
	if (toVerticesRemember.size()>0)
	{
		for (int i=0; i<(int)toVerticesRemember.size(); i++)
		{
			toDistanceArray[i] = ((GeoVertexData*)mesh->VertexData(mesh->Vertex(toVerticesRemember[i])))->distance;
			assert(toDistanceArray[i]>=0);			
		}
	}
	else
	{
		for (int i=0; i<mesh->NVertices(); i++)
		{
			toDistanceArray[i] = ((GeoVertexData*)mesh->VertexData(mesh->Vertex(i)))->distance;	
			assert(toDistanceArray[i]>=0);			
		}
	}
}

DistanceGeodesic::GeoVertexData * DistanceGeodesic::InitializeFunkPrecomputation(R3Mesh * mesh, bool verb, 
																				 GeoVertexData * vertex_data)
{
	if (verb)
		std::cout<<" Geodesics (Funk, Dijkstra) "<<std::endl;
	
	if (vertex_data==NULL)
		vertex_data = new GeoVertexData [ mesh->NVertices() ];
	assert(vertex_data!=NULL);
	
	// Initialize vertex data
	for (int i = 0; i < mesh->NVertices(); i++) 
	{
		R3MeshVertex *vertex = mesh->Vertex(i);
		mesh->SetVertexData(vertex, &(vertex_data[i]));
		vertex_data[i].vertex = vertex;	
	}
	
	return vertex_data;
}

void DistanceGeodesic::PrecomputeFunkhouser()
{
	R3Mesh * mesh = m_pSurface->GetMesh();
	// allocate vertex data
	GeoVertexData * vertex_data = InitializeFunkPrecomputation(mesh, m_bVerbous);

	if (m_bVerbous)	
		std::cout<<" ("<<mesh->NVertices()<<" x "<<mesh->NVertices()<<"): "<<flush;
		
	// Run Dijkstra for each from vertex
	for (int i=0; i<mesh->NVertices(); i++)
	{
		if (m_bVerbous)
		{	
			if (i%10000==0 && i!=0)
				std::cout<<"\t"<<i<<" / "<<mesh->NVertices()<<"\n\t\t"<<std::flush;
			else if (i%500==0 && i!=0)
				std::cout<<"."<<std::flush;
		}
		
		m_aPrecomputedDistances[i * m_iPrecompWidth + i] = 0;		
		PrecomputeFillFunkhouser(vertex_data, mesh, i, &(m_aPrecomputedDistances[i * m_iPrecompWidth]));

	}
	delete [] vertex_data;
	if (m_bVerbous)
		std::cout<<" - Done!"<<std::endl;
}

void DistanceGeodesic::InitializeKirsanovPrecomputation(geodesic::Mesh ** kirsanovMesh, 
														geodesic::GeodesicAlgorithmBase ** algorithm)
{
	R3Mesh * mesh = m_pSurface->GetMesh();
	if (m_bVerbous)
		std::cout<<"Geodesics (Kirsanov, "<<std::flush;
	
	// Initialize mesh in Kirsanov et al.'s format
	std::vector<double> points;
	std::vector<unsigned int> faces;
	
	for (int v=0; v<mesh->NVertices(); v++)
	{
		R3Point p = mesh->VertexPosition(mesh->Vertex(v));
		points.push_back(p.X());
		points.push_back(p.Y());
		points.push_back(p.Z());
	}
	
	for (int f=0; f<mesh->NFaces(); f++)
	{
		int v1 = mesh->VertexID(mesh->VertexOnFace(mesh->Face(f), 0));
		int v2 = mesh->VertexID(mesh->VertexOnFace(mesh->Face(f), 1));
		int v3 = mesh->VertexID(mesh->VertexOnFace(mesh->Face(f), 2));
		faces.push_back(v1);
		faces.push_back(v2);
		faces.push_back(v3);
	}
	
	*kirsanovMesh = new geodesic::Mesh();
	(**kirsanovMesh).initialize_mesh_data(points, faces);		//create internal mesh data structure including edges
	
	// Initialize appropriate algorithm
	switch(m_GeodesicAlgorithm)
	{
		case GEODIST_DIJKSTRA_KIRSANOV1:
			*algorithm = new geodesic::GeodesicAlgorithmDijkstra(*kirsanovMesh);
			if (m_bVerbous)
				std::cout<<", Dijkstra): "<<std::flush;
			break;
		case GEODIST_DIJKSTRA_KIRSANOV2:
			*algorithm = new geodesic::GeodesicAlgorithmDijkstraAlternative(*kirsanovMesh);	
			if (m_bVerbous)		
				std::cout<<", Dijkstra Alternative): "<<std::flush;			
			break;
		case GEODIST_EXACT_KIRSANOV:
			*algorithm = new geodesic::GeodesicAlgorithmExact(*kirsanovMesh);
			if (m_bVerbous)			
				std::cout<<", Exact): "<<std::flush;			
			break;			
		case GEODIST_SUBDIVISION_KIRSANOV:
			*algorithm = new geodesic::GeodesicAlgorithmSubdivision(*kirsanovMesh);			
			if (m_bVerbous)			
				std::cout<<", Subdivision): "<<std::flush;			
			break;
		default:
			assert(false);
	}
}

void DistanceGeodesic::FreeKirsanovPrecomputation(geodesic::Mesh * mesh, 
												  geodesic::GeodesicAlgorithmBase * algorithm)
{
	delete mesh;
	delete algorithm;	
}

void DistanceGeodesic::PrecomputeKirsanov()
{
	R3Mesh * mesh = m_pSurface->GetMesh();
	geodesic::Mesh * kirsanovMesh;
	geodesic::GeodesicAlgorithmBase * algorithm;
	InitializeKirsanovPrecomputation(&kirsanovMesh, &algorithm);
	
	assert(m_aPrecomputedDistances!=NULL);
	// run Kirsanov code for each vertex
	for (int i=0; i<mesh->NVertices(); i++)
	{
		if (m_bVerbous)
		{
			if (i%10000==0 && i!=0)
				std::cout<<"\t"<<i<<" / "<<mesh->NVertices()<<"\n\t\t"<<std::flush;
			else if (i%50==0 && i!=0)
				std::cout<<"."<<std::flush;
		}

		PrecomputeFillRowKirsanov(i, algorithm, kirsanovMesh, 
								  &(m_aPrecomputedDistances[i * m_iPrecompWidth]));
	}	
	FreeKirsanovPrecomputation(kirsanovMesh, algorithm);
}

void DistanceGeodesic::PrecomputeFillRowKirsanov(int fromVertex, 
									   geodesic::GeodesicAlgorithmBase * algorithm, 
									   geodesic::Mesh * mesh, double * toDistanceArray, 
									   const std::vector<int> & toVerticesRemember)
{
		// Create sources and propagate (e.g. find geodesics)
	geodesic::SurfacePoint source(&mesh->vertices()[fromVertex]);		//create source 
			//in general, there could be multiple sources, but now we have only one
	std::vector<geodesic::SurfacePoint> all_sources(1, source);	
				
	algorithm->propagate(all_sources);	//cover the whole mesh
	
	if (toVerticesRemember.size()>0)
	{
		for (int j=0; j<(int)toVerticesRemember.size(); j++)
		{
			double dist=0;
			geodesic::SurfacePoint p(&mesh->vertices()[toVerticesRemember[j]]);		
			algorithm->best_source(p, dist);	// returns best vertex ID - ignored
			toDistanceArray[j] = dist;
		}
	}
	else
	{
		for (int j=0; j<m_pSurface->GetMesh()->NVertices(); j++)
		{
			double dist=0;
			geodesic::SurfacePoint p(&mesh->vertices()[j]);		
			algorithm->best_source(p, dist);	// returns best vertex ID - ignored
			toDistanceArray[j] = dist;
		}
	}
}


DistanceGeodesic::CalculationType DistanceGeodesic::GetTypeFromStr(const VKString & calcType)
{
	if (calcType=="GEODIST_DIJKSTRA_FUNKHOUSER")
		return DistanceGeodesic::GEODIST_DIJKSTRA_FUNKHOUSER;
	else if (calcType=="GEODIST_DIJKSTRA_KIRSANOV1")
		return DistanceGeodesic::GEODIST_DIJKSTRA_KIRSANOV1;
	else if (calcType=="GEODIST_DIJKSTRA_KIRSANOV2")
		return DistanceGeodesic::GEODIST_DIJKSTRA_KIRSANOV2;
	else if (calcType=="GEODIST_EXACT_KIRSANOV")
		return DistanceGeodesic::GEODIST_EXACT_KIRSANOV;
	else if (calcType=="GEODIST_SUBDIVISION_KIRSANOV")
		return DistanceGeodesic::GEODIST_SUBDIVISION_KIRSANOV;
	else
		assert(false);
}






