#include "SurfaceDistance.h"
#include "SurfaceSampleSet.h"
#include "SampledSurface.h"
#include "DistanceGeodesic.h"
#include "AnalysisStats.h"

#include "FeatureMGD.h"

std::map<SampledSurface*, std::map<int, std::vector<int> > > FeatureMGD::m_EvenVertices;

FeatureMGD::FeatureMGD(SampledSurface * surface, const VKString & distanceName, 
					const VKString & initSampleSet,
					bool vertexOnly, const VKString & featurename)
: SurfaceFeature(surface, vertexOnly, featurename)
{
	m_pDistance = surface->GetDistanceMetric(distanceName);
	m_pSampleSet = surface->GetSampleSet(initSampleSet);
	assert(m_pSampleSet!=NULL);
}

FeatureMGD::FeatureMGD(SampledSurface * surface, const VKString & distanceName, 
		   SurfaceSampleSet * initSampleSet,
		   bool vertexOnly, const VKString & featurename )
: SurfaceFeature(surface, vertexOnly, featurename)
{
	m_pDistance = surface->GetDistanceMetric(distanceName);
	m_pSampleSet = initSampleSet;
}

VKString FeatureMGD::GetFeatureName()
{
	return "FeatureMGD";
}

double FeatureMGD::CalculateValue(const SurfaceSample & sample) const
{
	assert(m_pDistance!=NULL);
	double minDist=-1;
	for (int i=0; i<m_pSampleSet->NumSamples(); i++)
	{
		double dist=m_pDistance->Distance(m_pSampleSet->GetSample(i), sample);

		if (minDist<0 || minDist>dist)
			minDist=dist;
	}
	return minDist;
}

SurfaceSampleSet * FeatureMGD::ConstructSymmetryInvariantSet(SampledSurface * surface, int maxIter, 
															 double tau, int maxSamples, 
															 const VKString & initSetName,
															 const VKString & distanceName,
															 std::vector<SurfaceSampleSet *> * historySets,
															 std::vector<FeatureMGD *> * historyFeatures)
{
	SurfaceSampleSet * unionSet = new SurfaceSampleSet();
	assert(surface->GetSampleSet(initSetName)!=NULL);
	unionSet->Union(*surface->GetSampleSet(initSetName));
	
	SurfaceDistance * distance = surface->GetDistanceMetric(distanceName);
	if (distance!=NULL)
	{
		assert(false);	// parameter 1 is probably in wrong units - it should be geodesic distance
//		for (int i=0; i<maxIter; i++)
//		{
//			VKString featureName = VKString("MGD") + VKString::number(i);
//			FeatureMGD * newFeature = new FeatureMGD(surface, distanceName, unionSet, true, featureName);
//
//			SurfaceSampleSet * newSet = new SurfaceSampleSet(newFeature, 1, 0, 20,
//															 SurfaceSampleSet::FEAT_MAXIMA_ONLY);
//
//			if (newSet->NumSamples()==0)
//				break;
//			
//			double minValue = tau * newFeature->Value(newSet->GetSample(0));
//			for (int s=0; s<newSet->NumSamples(); s++)
//			{
//				if (newFeature->Value(newSet->GetSample(s))>=minValue)
//					unionSet->AddSample(newSet->GetSample(s));
//			}
//			
//			if (historySets!=NULL)
//				historySets->push_back(newSet);
//			else
//				delete newSet;
//				
//			if (historyFeatures!=NULL)
//				historyFeatures->push_back(newFeature);
//			else 
//				delete newFeature;
//			
//			if (unionSet->NumSamples()>=maxSamples)
//				break;
//		}
	}
	else
	{
		std::cout<<"Find MGD Set "<<std::endl;
		double * distancesToSet = new double[surface->GetMesh()->NVertices()];
		
		DistanceGeodesic::GeoVertexData * vertexData = DistanceGeodesic::InitializeFunkPrecomputation(surface->GetMesh(), false);
		double * distancesToSample = new double[surface->GetMesh()->NVertices()];
		
		for (int i=0; i<surface->GetMesh()->NVertices(); i++)
			distancesToSet[i] = -1;
		
		for (int i=0; i<unionSet->NumSamples(); i++)	// initialize distances to current set
		{
			bool ok;
			DistanceGeodesic::PrecomputeFillFunkhouser(vertexData, surface->GetMesh(), 
													   unionSet->GetSample(i).NearestVertex(&ok), 
													   distancesToSample);
			assert(ok);
			for (int j=0; j<surface->GetMesh()->NVertices(); j++)
			{
				if (distancesToSet[j]==-1 || distancesToSet[j] > distancesToSample[j])
					distancesToSet[j] = distancesToSample[j];
			}
		}
		for (int i=0; i<maxIter && unionSet->NumSamples()<maxSamples; i++)
		{
			if (i%10==0)
				AnalysisStats::m_GlobalStats.m_Timing.WriteProgress("CreatingSurface_Samples",
																	i, maxSamples);
//			if (i%50==0)
//				std::cout<<"."<<std::flush;
			int currCountUnionSet = unionSet->NumSamples();
			double farthestDistFromSet = -1;
			
				// find threshold for samples to add
			for (int j=0; j<surface->GetMesh()->NVertices(); j++)	
				if (farthestDistFromSet < distancesToSet[j])
					farthestDistFromSet = distancesToSet[j];
			double minValue = tau * farthestDistFromSet;
			
			int added=0;
				// add samples
			for (int j=0; j<surface->GetMesh()->NVertices(); j++)	
			{
				if (distancesToSet[j] >= minValue)
				{
					added++;
					unionSet->AddSample(SurfaceSample(j, surface->GetMesh()));
				}
			}
	//		std::cout<<"Iteration["<<i<<"]. Added="<<added<<std::endl;
			
				// update distances
			for (int j=currCountUnionSet; j<unionSet->NumSamples(); j++)	
			{
				bool ok;
				DistanceGeodesic::PrecomputeFillFunkhouser(vertexData, surface->GetMesh(), 
														   unionSet->GetSample(j).NearestVertex(&ok), 
														   distancesToSample);
				assert(ok);
				
				for (int j=0; j<surface->GetMesh()->NVertices(); j++)
				{
					if (distancesToSet[j]==-1 || distancesToSet[j] > distancesToSample[j])
						distancesToSet[j] = distancesToSample[j];
				}				
			}
		}
//		std::cout<<" - done!"<<std::endl;
		std::cout<<std::endl;
		assert(unionSet->NumSamples() < surface->GetMesh()->NVertices());
		
		delete [] distancesToSet;
		delete [] vertexData;
		delete [] distancesToSample;
	}
	return unionSet;
}

void FeatureMGD::EvenlySpreadVertices(SampledSurface * surface, int numPnts, 
									  std::vector<int> & evenVertices)
{
	std::vector<int> seeds;
	
	std::map<SampledSurface*, std::map<int, std::vector<int> > >::iterator iter1;
	iter1 = m_EvenVertices.find(surface);
	if (iter1 != m_EvenVertices.end())
	{
		std::map<int, std::vector<int> >::iterator iter2 = iter1->second.find(numPnts);
		if (iter2 != iter1->second.end())
		{
			evenVertices = iter2->second;
			return;
		}
		else
		{
			int distToPrev = -1;
			int distToNext = -1;			
			std::map<int, std::vector<int> >::iterator prev;
			std::map<int, std::vector<int> >::iterator next;
			for (iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)
			{
				if (iter2->first < numPnts)
				{
					int currDist = numPnts - iter2->first;
					if (currDist < distToPrev || distToPrev == -1)
					{
						distToPrev = currDist;
						prev = iter2;
					}
				}
				else
				{
					assert(iter2->first > numPnts);
					int currDist = iter2->first - numPnts;
					if (currDist < distToNext || distToNext == -1)
					{
						distToNext = currDist;
						next = iter2;
					}					
				}
			}
			
			assert(distToNext>=0 || distToPrev>=0);
			if (distToNext>=0)	// iteratively remove from maximal iter2
			{
				std::vector<int> verticesFromNext = next->second;
				assert((int)verticesFromNext.size() > numPnts);
				while ((int)verticesFromNext.size() > numPnts)
					verticesFromNext.erase(verticesFromNext.begin()+numPnts);
				evenVertices = verticesFromNext;
				return;
			}
			else if (distToPrev>=0)			// initialize seeds to maximal iter2
			{
				seeds = prev->second;
				EvenlySpreadVertices(surface, numPnts, evenVertices, seeds);
			}
			else
				assert(false);
			return;
		}
	}
	else
	{
		// nothing for this surface - initialize random seed
		int randSeed = rand() % surface->GetMesh()->NVertices();
		seeds.push_back(randSeed < 0 ? -randSeed : randSeed);
		std::vector<int> newEvenVertexIDs;
		EvenlySpreadVertices(surface, numPnts, newEvenVertexIDs, seeds);
		evenVertices = newEvenVertexIDs;
		m_EvenVertices[surface][numPnts] = newEvenVertexIDs;
		return;
	}
}


void FeatureMGD::EvenlySpreadVertices(SampledSurface * surface, int numPnts, 
									  std::vector<int> & evenVertices, 
									  std::vector<int> & argSeeds)
{
	R3Mesh * mesh = surface->GetMesh();
	// Allocate vertex data
	VertexData *vertexData = new VertexData [ mesh->NVertices() ];
	assert (vertexData!=NULL);
	
	// Initialize vertex data
	for (int i = 0; i < mesh->NVertices(); i++) 
	{
		R3MeshVertex *vertex = mesh->Vertex(i);
		VertexData *data = &(vertexData[i]);
		mesh->SetVertexData(vertex, data);
		data->vertex = vertex;
		data->selected = FALSE;
		data->distance = FLT_MAX;
		data->heappointer = NULL;
	}
	
	std::cout<<"Evenly Spread Vertices:"<<std::endl;
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("FeatureCalc_MGD");
	std::vector<R3MeshVertex * > seeds;
//	seeds.push_back(surface->GetMesh()->Vertex(0));
	for (int i=0; i<(int)argSeeds.size(); i++)
		seeds.push_back(surface->GetMesh()->Vertex(argSeeds[i]));
	assert(seeds.size()>0);
	for (int i=0; i<numPnts; i++)
	{
		if (i%10==0)
			AnalysisStats::m_GlobalStats.m_Timing.WriteProgress("FeatureCalc_MGD", i, numPnts);		
//		if (i%5000==0 && i!=0)
//			std::cout<<" "<<i<<" / "<<numPnts<<std::endl;
//		else if (i%200==0)
//			std::cout<<"."<<std::flush;
		R3MeshVertex * v = FindFurthestVertex(mesh, seeds, vertexData);
		
		if (v==NULL)
			break;
		
		if (i==0)
			seeds.clear();
		
		VertexData *data = &vertexData[mesh->VertexID(v)];		
		seeds.push_back(v);
		data->selected = true;
		evenVertices.push_back(mesh->VertexID(v));
	}
	std::cout<<std::endl;
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("FeatureCalc_MGD");	
	
	delete [] vertexData;
}

R3MeshVertex * FeatureMGD::FindFurthestVertex(R3Mesh *mesh, const std::vector<R3MeshVertex *>& seeds, 
													 VertexData *vertex_data)
{
	// Initialize vertex data
	for (int i = 0; i < mesh->NVertices(); i++) {
		VertexData *data = &vertex_data[i];
		data->distance = FLT_MAX;
		data->heappointer = NULL;
	}
	
	// Initialize heap
	VertexData tmp;
	RNHeap<VertexData *> heap(&tmp, &(tmp.distance), &(tmp.heappointer));
	for (int i = 0; i < (int)seeds.size(); i++) 
	{
		R3MeshVertex *vertex = seeds[i];
		VertexData *data = &vertex_data[mesh->VertexID(vertex)];
		data->distance = 0;
		heap.Push(data);
	}
	
	// Iteratively pop off heap until find further vertex
	R3MeshVertex *vertex = NULL;
	while (!heap.IsEmpty()) 
	{
		VertexData *data = heap.Pop();
		vertex = data->vertex;
		for (int j = 0; j < mesh->VertexValence(vertex); j++) 
		{
			R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, j);
			R3MeshVertex *neighbor_vertex = mesh->VertexAcrossEdge(edge, vertex);
			VertexData *neighbor_data = (VertexData *) mesh->VertexData(neighbor_vertex);
			RNScalar old_distance = neighbor_data->distance;
			RNScalar new_distance = mesh->EdgeLength(edge) + data->distance;
			if (new_distance < old_distance) 
			{
				neighbor_data->distance = new_distance;
				if (old_distance < FLT_MAX) heap.Update(neighbor_data);
				else heap.Push(neighbor_data);
			}
		}
	}
	
	// Return furthest vertex
	return vertex;
}


