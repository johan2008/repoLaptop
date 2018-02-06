#include "SampledSurface.h"
#include "SurfaceFeature.h"
#include "FeatureAGD.h"
#include "FeatureMGD.h"
#include "DistanceGeodesic.h"
#include "SurfaceDistance.h"
#include "SurfaceMap.h"
#include "VkFunctions.h"
#include "SurfaceMapConfidence.h"
#include "SurfaceMapSimilarity.h"
#include "AnalysisStats.h"
#include "DistanceOnTheFly.h"
#include "MapMultiConformal.h"
#include "MapConfidenceMultiConf.h"
#include "MeshProcessor.h"
#include "AnalysisStats.h"
#include "MapConfidenceCompareToTruth.h"
#include "SurfaceTexture.h"
#include "PipelineGeneral.h"

////////////////////// BASIC ROUTINES /////////////////////////////
SampledSurface::SampledSurface(R3Mesh * mesh)
{
	m_pMeshSearchTree = NULL;
	m_fArea = -1;
	m_bCachedBBox = false;
	m_sSurfaceName = "none";
	m_iCameraID = -1;
	if (mesh==Surface2DPlane::m_pInfinitePlanePseudomesh)
		m_pMesh = mesh;
	else
		m_pMesh = CheckMeshForProblems(mesh);
	
	m_iCachedVertexFrom = -1;
	m_sCachedNameMapOnSurf="none";		
	m_CachedMainColorSet="";	
	m_pCachedSurfaceFeature=NULL;
	m_pCachedSurfaceDistance=NULL;
	m_pCachedSurfaceMap=NULL;
	m_pCachedSurfaceOtherMap=NULL;
	m_sCachedMapSimilarity="none";
	m_sCachedMapConfidence = "none";
	m_iCachedConfID = -1;
}

void SampledSurface::AddMapActingAsFrom(const VKString & mapName, SurfaceMap * surfMap )
{
	m_mNameToMapActingAsFrom[mapName] = surfMap;
}

SurfaceMap * SampledSurface::GetMapActingAsFrom(const VKString & mapName)
{
	if (m_mNameToMapActingAsFrom.find(mapName)!=m_mNameToMapActingAsFrom.end())
		return m_mNameToMapActingAsFrom[mapName];
	return NULL;
}

void SampledSurface::AddMapActingAsTo(const VKString & mapName, SurfaceMap * surfMap )
{
	m_mNameToMapActingAsTo[mapName] = surfMap;
}

SurfaceMap * SampledSurface::GetMapActingAsTo(const VKString & mapName)
{
	if (m_mNameToMapActingAsTo.find(mapName)!=m_mNameToMapActingAsTo.end())
		return m_mNameToMapActingAsTo[mapName];
	return NULL;	
}


void SampledSurface::AddSampleSet(const VKString & setName, SurfaceSampleSet * set)
{
	m_mNameToSampleSet[setName] = set;
}

SurfaceSampleSet * SampledSurface::GetSampleSet(const VKString & setName)
{
	if(m_mNameToSampleSet.find(setName)==m_mNameToSampleSet.end())
	{
		if (setName=="AllVertices")
		{
			SurfaceSampleSet * allVertices = new SurfaceSampleSet(GetMesh());
			AddSampleSet("AllVertices", allVertices);
			return allVertices;
		}
		return NULL;
	}
	else
		return m_mNameToSampleSet[setName];
}


void SampledSurface::AddFeature(const VKString & featureName, SurfaceFeature * feature)
{
	m_mNameToFeatureSet[featureName] = feature;
}

SurfaceFeature * SampledSurface::GetFeature(const VKString & featureName)
{
	if (m_mNameToFeatureSet.find(featureName)==m_mNameToFeatureSet.end())
		return NULL;
	else
		return m_mNameToFeatureSet[featureName];
}

void SampledSurface::AddDistanceMetric(const VKString & metricName, SurfaceDistance * distances)
{
	m_mNameToDistance[metricName] = distances;
}

SurfaceDistance * SampledSurface::GetDistanceMetric(const VKString & metricName)
{
	if (m_mNameToDistance.find(metricName)==m_mNameToDistance.end())
		return NULL;
	else
		return m_mNameToDistance[metricName];
}

DistanceOnTheFly * SampledSurface::GetOnTheFlyDistanceMetric(double maxDistanceUnscaled, 
															const VKString & tryMeFirstMetric, int maxCache)
{
	VKString metricName = VKString("FLY_")+tryMeFirstMetric+VKString::number(maxDistanceUnscaled)+"_"+VKString::number(maxCache);
	SurfaceDistance * retDistance = GetDistanceMetric(metricName);
	if (retDistance==NULL)
	{
		if (maxDistanceUnscaled<0)
			retDistance = new DistanceOnTheFly(this, DistanceGeodesic::GEODIST_DIJKSTRA_FUNKHOUSER,
											   GetDistanceMetric(tryMeFirstMetric), maxCache);
		else 
			retDistance = new DistanceOnTheFly(this, DistanceGeodesic::GEODIST_DIJKSTRA_FUNKHOUSER,
											   GetDistanceMetric(tryMeFirstMetric), -1, 
											   AdjustedRadius(maxDistanceUnscaled));
		
		AddDistanceMetric(metricName, retDistance);
	}
	return (DistanceOnTheFly*)retDistance;
}

R3Mesh * SampledSurface::GetMesh()
{
	return m_pMesh;
}

SurfaceSample SampledSurface::WeightedGeodesicCentroid(std::vector<SurfaceSample> & samples,
													   std::vector<double> & weights,
													   GeodesicCentroidAlgorithmType type, 
													   double outlierThreshold)
{	
	MapConfidenceMultiConf::NormalizeWeights(weights);	
	if (type==GEO_CENTROID_EVERY_VERTEX)
		return WeightedGeoCentrEveryVertex(samples, weights);
	else if (type==GEO_CENTROID_NEIGHBORHOOD)
		return WeightedGeoCentrNhdOnly(samples, weights, 0.25);
	else if (type==GEO_CENTROID_BEST_AMONG_GIVEN)
		return WeightedGeoCentrBestAmongGiven(samples, weights);
	else if (type==GEO_CENTROID_EUCLIDEAN_APPROX)
		return WeightedGeoCentrEuclideanApprox(samples, weights, outlierThreshold);
	else
		assert(false);
}

double SampledSurface::GetWeightedGeoCentroidError(const SurfaceSample & s, 
												   std::vector<SurfaceSample> & samples, 
												   std::vector<double> & weights, 
												   SurfaceDistance * distanceMetric)
{
	double err = 0;
	for (int i=0; i<(int)samples.size(); i++)
		err += weights[i] * distanceMetric->Distance(samples[i], s);
		
	return err;
}

SurfaceSample SampledSurface::WeightedGeoCentrEveryVertex(std::vector<SurfaceSample> & samples,
														  std::vector<double> & weights)
{
	double minErr = FLT_MAX;
	int minErrVertex=-1;

	DistanceOnTheFly * onTheFlyDist = GetOnTheFlyDistanceMetric(-1, "default", -1);
	
	for (int j=0; j<(int)samples.size(); j++)		
	{
		bool atVertex;
		int vertexID = samples[j].NearestVertex(&atVertex);
		if (!atVertex)
		{
			onTheFlyDist->PrecomputeRow(samples[j].VertID(0));
			onTheFlyDist->PrecomputeRow(samples[j].VertID(1));
			onTheFlyDist->PrecomputeRow(samples[j].VertID(2));
		}
		else 
			onTheFlyDist->PrecomputeRow(vertexID);
	}
	
	for (int i=0; i<GetMesh()->NVertices(); i++)
	{
		SurfaceSample vertexSamp(i, GetMesh());
		double currErr = GetWeightedGeoCentroidError(vertexSamp, samples, weights, onTheFlyDist);
		
		if (currErr < minErr)
		{
			minErr = currErr;
			minErrVertex = i;
		}
	}
	
	return SurfaceSample(minErrVertex, GetMesh());
}

SurfaceSample SampledSurface::WeightedGeoCentrNhdOnly(std::vector<SurfaceSample> & samples,
													  std::vector<double> & weights, double thresh)
{	
	double minErr = FLT_MAX;
	int minErrVertex=-1;

	DistanceOnTheFly * onTheFlyDist = GetOnTheFlyDistanceMetric(thresh, "default", -1);
	double maxDistance = AdjustedRadius(thresh);		
	
	for (int j=0; j<(int)samples.size(); j++)
		if (weights[j] < .5 / samples.size())
			weights[j] = 0;
	
	std::vector<std::map<int, double> > sampleDistanceMaps;
	std::map<int, double> vertexValues;
	
	for (int j=0; j<(int)samples.size(); j++)	// for each sample find neighborhoods
	{
		sampleDistanceMaps.push_back(std::map<int, double>());			
		if (weights[j]==0)
			continue;
		// search in neighborhood of each vertex in triangle of a sample
		int v1 = samples[j].VertID(0);
		int v2 = samples[j].VertID(1);
		int v3 = samples[j].VertID(2);
		//			AnalysisStats::m_GlobalStats.m_Timing.startedProcess("ZCalculatingCentroidDistanceRelated1");
		onTheFlyDist->PrecomputeRow(v1);
		onTheFlyDist->PrecomputeRow(v2);
		onTheFlyDist->PrecomputeRow(v3);			
		//			AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("ZCalculatingCentroidDistanceRelated1");
		
		// now calculate distances from samples to neighbohoods
		std::map<int, double> & fromVert1 = onTheFlyDist->GetNeighborhood(v1);
		//			fromVert1[v1] = 0;
		for (std::map<int, double>::iterator iter=fromVert1.begin(); 
			 iter!=fromVert1.end(); iter++)
		{
			//					// if vertex is defined (close enough) to each vertex of sample's triangle: add vertex
			//				AnalysisStats::m_GlobalStats.m_Timing.startedProcess("ZCalculatingCentroidDistanceRelated2");				
			int toVertex = iter->first;
			double d1 = DistanceGeodesic::GetThresholdedDistances(toVertex, GetMesh(), fromVert1, 
																  maxDistance);
			//				if (d1==maxDistance)
			//					continue;
			double d2 = DistanceGeodesic::GetThresholdedDistances(toVertex, GetMesh(), 
																  onTheFlyDist->GetNeighborhood(v2), 
																  maxDistance);
			//				if (d2==maxDistance)
			//					continue;
			double d3 = DistanceGeodesic::GetThresholdedDistances(toVertex, GetMesh(), 
																  onTheFlyDist->GetNeighborhood(v3), 
																  maxDistance);
			//				if (d3==maxDistance)
			//					continue;
			
			sampleDistanceMaps[j][toVertex] = samples[j].Interpolate(d1, d2, d3);
			vertexValues[toVertex] = 0;
		}
	}
	
	//		AnalysisStats::m_GlobalStats.m_Timing.startedProcess("ZCalculatingCentroidFindingmin");		
	// find minimum
	//		std::cout<<"values = "<<vertexValues.size()<<std::endl;
	for (std::map<int, double>::iterator iter = vertexValues.begin(); iter!=vertexValues.end(); iter++)
	{
		double currErr = 0;
		
		for (int j=0; j<(int)samples.size(); j++)	// for each sample add weight
		{
			if (weights[j]==0)
				continue;
			double dist = maxDistance;				
			std::map<int, double>::iterator iter2 = sampleDistanceMaps[j].find(iter->first);
			if (iter2!=sampleDistanceMaps[j].end())
				dist = iter2->second;
			
			currErr += weights[j] * dist;
		}
		
		if (currErr < minErr)
		{
			minErr = currErr;
			minErrVertex = iter->first;
		}
	}
	
	return SurfaceSample(minErrVertex, GetMesh());
}

SurfaceSample SampledSurface::WeightedGeoCentrBestAmongGiven(std::vector<SurfaceSample> & samples,
															 std::vector<double> & weights)
{
	assert(false);	// NOTE: also needs non-trivial geodesic distances
//	int bestSampleID=0;
//	double bestSampleVal = GetWeightedGeoCentroidError(samples[0], samples, weights, 
//													   distanceMetric);
//	for (int i=0; i<(int)samples.size(); i++)
}

SurfaceSample SampledSurface::WeightedGeoCentrEuclideanApprox(std::vector<SurfaceSample> & samples,
															  std::vector<double> & weights,
															  double outlierThreshold)
{
	// ignore small weights for computational efficiency
	MapConfidenceMultiConf::NormalizeWeights(weights);
	//std::cout<<"outleir: "<<outlierThreshold<<std::endl;
	double threshold = outlierThreshold / (double)samples.size();
	int bestID=0;
	double bestWeight = -1000000;
	int nzEntries=0;
	for (int i=0; i<(int)weights.size(); i++)
	{
		R3Point p = samples[i].GetPosition();
		if (bestWeight < weights[i])
		{
			bestWeight = weights[i];
			bestID = i;
		}
		if (weights[i] < threshold)
			weights[i] = 0;
		else
			nzEntries++;
	}
	if (nzEntries==0)
	{
		weights[bestID] = 1.;
		//std::cout<<"picking best: "<<bestID<<std::endl;
	}
	else
		MapConfidenceMultiConf::NormalizeWeights(weights);
	
	// find euclidean centroid
	R3Point euclCentr(0,0,0);
	for (int i=0; i<(int)samples.size(); i++)
	{
		if (weights[i] > 0)
		{
			R3Point p = samples[i].GetPosition();
			euclCentr = euclCentr + p*weights[i];
		}
	}
	double minDist=FLT_MAX;
	int minID=-1;
	for (int i=0; i<(int)samples.size(); i++)
	{
		if (weights[i]>0)
		{
			double dist = (euclCentr-samples[i].GetPosition()).Length();
			if (dist < minDist)
			{
				minDist = dist;
				minID = i;
			}
		}
	}
	if (minID==-1)
		return samples[0];
	if (minDist==0)
		return samples[minID];
	minDist*=1.001;
	
	SurfaceSample geoCentroid = GetNearestSample(euclCentr, 0, minDist);
	if (geoCentroid.Invalid())
		return samples[minID];
	
//	static bool singlePrint = true;
//	if (singlePrint)
//	{
//		for (int i=0; i<(int)samples.size(); i++)
//			std::cout<<"\tSample["<<i<<"]: vID="<<samples.VertID()<<" w="<<weights[i]<<std::endl;
//		
//		singlePrint = false;
//		std::cout<<"centroid: faceID="<<geoCentroid.TriID()<<" bary=["<<geoCentroid.B(0)<<", ";
//		std::cout<<geoCentroid.B(1)<<", "<<geoCentroid.B(2)<<"] "<<std::endl;
//	}
	
	return geoCentroid;
}

SurfaceSample SampledSurface::GetNearestSample(const R3Point & queryPnt, 
											   double minDist, double maxDist)
{
	// create search tree
	if (m_pMeshSearchTree==NULL)
		m_pMeshSearchTree = new R3MeshSearchTree(m_pMesh);
	
	// project euclidean centroid
	R3MeshIntersection closest;
	m_pMeshSearchTree->FindClosest(queryPnt, closest, minDist, maxDist);
	if (closest.type==R3_MESH_NULL_TYPE)
		return SurfaceSample();
	else
		return SurfaceSample(closest, closest.point, m_pMesh);
}


void SampledSurface::SetParentWindow(AnalysisWindow * window)
{
	m_pWindow = window;
	m_iCameraID = m_pWindow->AddMesh(GetMesh());
}

VKStringList SampledSurface::GetAllFeatureNames()
{
	VKStringList retList;
	for (std::map<VKString, SurfaceFeature *>::iterator iter = m_mNameToFeatureSet.begin();
		 iter!=m_mNameToFeatureSet.end(); iter++)
		if (iter->first!="default" || (int)m_mNameToFeatureSet.size()==1)
			retList.push_back(iter->first);
	return retList;
}

VKStringList SampledSurface::GetAllSampleSetNames()
{
	VKStringList retList;
	for (std::map<VKString, SurfaceSampleSet *>::iterator iter = m_mNameToSampleSet.begin();
		 iter!=m_mNameToSampleSet.end(); iter++)
		if (iter->first!="default" || (int)m_mNameToSampleSet.size()==1)		
			retList.push_back(iter->first);	
	return retList;
}

VKStringList SampledSurface::GetAllDistanceMetricNames()
{
	VKStringList retList;
	for (std::map<VKString, SurfaceDistance *>::iterator iter = m_mNameToDistance.begin();
		 iter!=m_mNameToDistance.end(); iter++)
		if (iter->first!="default" || (int)m_mNameToDistance.size()==1)		
			retList.push_back(iter->first);		
	return retList;
}

void SampledSurface::SelectFeature(const VKString & featureName)
{
	assert(m_mNameToFeatureSet.find(featureName)!=m_mNameToFeatureSet.end());
	m_mNameToFeatureSet["default"] = m_mNameToFeatureSet[featureName];	
}

void SampledSurface::SelectDistanceMetric(const VKString & distanceName)
{
	assert(m_mNameToDistance.find(distanceName)!=m_mNameToDistance.end());
	m_mNameToDistance["default"] = m_mNameToDistance[distanceName];
}

void SampledSurface::SelectSampleSet(const VKString & samplesName)
{
	assert(m_mNameToSampleSet.find(samplesName)!=m_mNameToSampleSet.end());
	m_mNameToSampleSet["default"] = m_mNameToSampleSet[samplesName];
}

double SampledSurface::AdjustedRadius(double fractionOfArea)
{
	if (fractionOfArea<0)
		return 0;
	return sqrt(fractionOfArea * Area() / 3.14159265);
}

double SampledSurface::Area()
{
	if (m_fArea>0)
		return m_fArea;

	if (this==Surface2DPlane::m_pInfinitePlanePseudosurface)
		m_fArea = 1.;
	else
	{
		m_fArea = 0;
		for (int i=0; i<m_pMesh->NFaces(); i++)
			m_fArea += m_pMesh->FaceArea(m_pMesh->Face(i));
	}
	return m_fArea;
}

VKString SampledSurface::GetSurfaceType()
{
	return "Surface";
}

VKString SampledSurface::GetName()
{
	return m_sSurfaceName;
}

void SampledSurface::SetName(const VKString & name)
{
	m_sSurfaceName = name;
}

void SampledSurface::CopyFromSameSurface(SampledSurface * sameSurface)
{
	for (std::map<VKString, SurfaceSampleSet *>::iterator iter=sameSurface->m_mNameToSampleSet.begin();
		 iter != sameSurface->m_mNameToSampleSet.end(); iter++)
		AddSampleSet(iter->first, iter->second->GetCopyAnotherMesh(GetMesh()));
	
	for (std::map<VKString, SurfaceFeature *>::iterator iter=sameSurface->m_mNameToFeatureSet.begin();
		 iter != sameSurface->m_mNameToFeatureSet.end(); iter++)
		AddFeature(iter->first, iter->second);
	
	for (std::map<VKString, SurfaceDistance *>::iterator iter=sameSurface->m_mNameToDistance.begin();
		 iter != sameSurface->m_mNameToDistance.end(); iter++)
		AddDistanceMetric(iter->first, iter->second);
}

R3Mesh * SampledSurface::CreateCopy(R3Mesh * mesh)
{
	R3Mesh * copyMesh = new R3Mesh();
	for (int i=0; i<mesh->NVertices(); i++)
		copyMesh->CreateVertex(mesh->VertexPosition(mesh->Vertex(i)));
	
	for (int f=0; f<mesh->NFaces(); f++)
	{
		int v1 = mesh->VertexID(mesh->VertexOnFace(mesh->Face(f), 0));
		int v2 = mesh->VertexID(mesh->VertexOnFace(mesh->Face(f), 1));
		int v3 = mesh->VertexID(mesh->VertexOnFace(mesh->Face(f), 2));		
		copyMesh->CreateFace(copyMesh->Vertex(v1), copyMesh->Vertex(v2), copyMesh->Vertex(v3));
	}
	assert(copyMesh->NFaces()==mesh->NFaces() && mesh->NVertices()==copyMesh->NVertices());
	return copyMesh;
}

////////////////////// RENDERING ROUTINES /////////////////////////////
void SampledSurface::Draw(ParamParser * params, 
						  const VKString & renderingParams, 
						  const VKString & surfaceName)
{
	DrawSurface(params, renderingParams, surfaceName);
}

double SampledSurface::GetStandardRadius(ParamParser * , 
										 const VKString & ,
										 const VKString & )
{
	return 0.005 * GetMesh()->BBox().DiagonalLength();
}

void SampledSurface::LoadCamera()
{
	m_pWindow->LoadCamera(m_iCameraID);
}

void SampledSurface::DrawSurface(ParamParser * params, 
								 const VKString & renderingParams,
								 const VKString & surfaceName)
{
	LoadCamera();

	// Initialize mesh, materials, ligths, etc...
	static GLfloat material[4];
	
	R3Mesh * mesh = m_pMesh;
	m_pWindow->SetLights(m_pMesh);
	
	// Backfacing
	if (m_pWindow->IsFlagSet(renderingParams, "MeshFlags", "Backfacing", 
							 surfaceName, "MeshFlagsFlips"))
		glDisable(GL_CULL_FACE);		
	else 
		glEnable(GL_CULL_FACE);	

	if (m_pWindow->IsFlagSet(renderingParams, "MeshFlags", "Axis", 
							 surfaceName, "MeshFlagsFlips"))
	{
		RNScalar d = mesh->BBox().DiagonalRadius();
		glDisable(GL_LIGHTING);
		glLineWidth(3);
		R3BeginLine();
		glColor3f(1, 0, 0);
		R3LoadPoint(R3zero_point + d * R3negx_vector);
		R3LoadPoint(R3zero_point + d * R3posx_vector);
		R3EndLine();
		R3BeginLine();
		glColor3f(0, 1, 0);
		R3LoadPoint(R3zero_point + d * R3negy_vector);
		R3LoadPoint(R3zero_point + d * R3posy_vector);
		R3EndLine();
		R3BeginLine();
		glColor3f(0, 0, 1);
		R3LoadPoint(R3zero_point + d * R3negz_vector);
		R3LoadPoint(R3zero_point + d * R3posz_vector);
		R3EndLine();
		glLineWidth(1);		
	}

	VKString drawType = params->GetStrValue(renderingParams, "SignalOnSurface", valid);	
	VKString textureName = params->GetStrValue(renderingParams, "RenderTexture", valid);
	bool renderTextures = textureName!="none" && (drawType=="auto" || drawType=="Texture") 
						&& params->GetStrValue(renderingParams, "RenderMap", valid)!="none";
	bool renderColorsWithLighting = params->GetStrValue(renderingParams, "EnableLighting", valid)=="true";
	
//	std::cout<<"renderTextures="<<renderTextures<<std::endl;
//	std::cout<<"renderTexture="<<renderTextures<<std::endl;
//	std::cout<<"textureName="<<textureName.c_str()<<std::endl;
//	std::cout<<"drawType="<<drawType.c_str()<<std::endl;
//	std::cout<<"renderMap="<<params->GetStrValue(renderingParams, "RenderMap", valid).c_str()<<std::endl;
	
	// render Faces and functions on faces (distances, features, etc...)
	if (m_pWindow->IsFlagSet(renderingParams, "MeshFlags", "Faces", 
										  surfaceName, "MeshFlagsFlips"))
	{
		glEnable(GL_LIGHTING);
		m_pWindow->StandardSurfaceColor(material);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material); 
		
		if (m_sCachedNameMapOnSurf!="none")
		{
			assert((int)m_CachedPerVertexColors.size() == 3*m_pMesh->NVertices());
			bool allBlack = true;			
			for (int i=0; i<(int)m_CachedPerVertexColors.size() && allBlack; i++)
			{
				if (m_CachedPerVertexColors[i]!=0)
					allBlack = false;
			}
			if (allBlack)
			{
				glEnable(GL_LIGHTING);
				m_pMesh->DrawFaces();
			}
			else
			{
				if (renderTextures)
				{
					glEnable(GL_TEXTURE_2D);
					assert(PipelineGeneral::s_pGeneralPipeline!=NULL);
					SurfaceTexture * texture = PipelineGeneral::s_pGeneralPipeline->GetTextureByName(textureName);
					if (texture==NULL)
						renderTextures = false;
					else
						texture->LoadTexture();
				}
				
				if (renderColorsWithLighting || renderTextures)
					glEnable(GL_LIGHTING);
				else
					glDisable(GL_LIGHTING);

				glColor3d(1., 1., 1.);
				material[0] = 1.;				
				material[1] = 1.;
				material[2] = 1.;
				material[3] = 1.;
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material); 
				glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material);

				
				glBegin(GL_TRIANGLES);			
				for (int i=0; i<m_pMesh->NFaces(); i++)
				{
					for (int v=0; v<3; v++)
					{
						int vertexID = m_pMesh->VertexID(m_pMesh->VertexOnFace(m_pMesh->Face(i), v));
						R3Point vertexPos = m_pMesh->VertexPosition(m_pMesh->Vertex(vertexID));
						
						if (renderTextures)
						{
							glTexCoord2d(m_CachedPerVertexColors[vertexID*3+0],
										 m_CachedPerVertexColors[vertexID*3+1]);
						}
						else if (renderColorsWithLighting)
						{
							for (int mid=0; mid<3; mid++)
								material[mid] = m_CachedPerVertexColors[vertexID*3+mid];
						}
						else
						{
							glColor3d(m_CachedPerVertexColors[vertexID*3+0], 
									  m_CachedPerVertexColors[vertexID*3+1], 
									  m_CachedPerVertexColors[vertexID*3+2]);
						}
						glVertex3d(vertexPos.X(), vertexPos.Y(), vertexPos.Z());
					}
				}
				glEnd();	
				if (renderTextures)
					glDisable(GL_TEXTURE_2D);
			}
		}
		else
		{
			m_pMesh->DrawFaces();
		}
	}
	
	// render face IDs
	if (m_pWindow->IsFlagSet(renderingParams, "MeshFlags", "FaceIDs", 
											surfaceName, "MeshFlagsFlips"))
	{
		glDisable(GL_LIGHTING);
		glColor3f(0, 0, 0);
		for (int i = 0; i < mesh->NFaces(); i++) 
		{
			R3MeshFace *face = mesh->Face(i);
			m_pWindow->RenderText(mesh->FaceCentroid(face), 
								  VKString::number(mesh->FaceID(face)), 0, 0, 0);
		}
	}
	
	// Render edges and IDs
	if (m_pWindow->IsFlagSet(renderingParams, "MeshFlags", "Edges", surfaceName, "MeshFlagsFlips"))
	{
		glLineWidth(1);
		glDisable(GL_LIGHTING);
		glColor3f(1.0, 0.0, 0.0);
		mesh->DrawEdges();		
	}
	
	if (m_pWindow->IsFlagSet(renderingParams, "MeshFlags", "EdgeIDs", surfaceName, "MeshFlagsFlips"))
	{
		glDisable(GL_LIGHTING);
		glColor3f(0.8, 0.0, 0.0);
		for (int i = 0; i < mesh->NEdges(); i++) 
		{
			R3MeshEdge *edge = mesh->Edge(i);
			m_pWindow->RenderText(mesh->EdgeMidpoint(edge), 
								  VKString::number(mesh->EdgeID(edge)), .8, 0, 0);
		}
	}
	
	// render vertices and IDs
	if (m_pWindow->IsFlagSet(renderingParams, "MeshFlags", "Vertices", 
							 surfaceName, "MeshFlagsFlips"))
	{
		glEnable(GL_LIGHTING);
		double radius = GetStandardRadius(params, renderingParams, surfaceName);
		for (int i = 0; i < mesh->NVertices(); i++) 
		{
			R3MeshVertex *vertex = mesh->Vertex(i);
			R3Point position = mesh->VertexPosition(vertex);
			m_pWindow->MapIntToColor(i, material);
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material); 
			R3Sphere(position, radius).Draw();
		}
		
	}

	if (m_pWindow->IsFlagSet(renderingParams, "MeshFlags", "VertexIDs", 
							 surfaceName, "MeshFlagsFlips"))
	{
		glDisable(GL_LIGHTING);
		glColor3f(0.5, 0.3, 0.1);
		for (int i = 0; i < mesh->NVertices(); i++) 
		{
			R3MeshVertex *vertex = mesh->Vertex(i);
			m_pWindow->RenderText(mesh->VertexPosition(vertex), 
								  VKString::number(mesh->VertexID(vertex)), .5, .3, .1);
		}		
	}
	
	int selectedVertex = GetSelectedVertex(params);
	
	if (selectedVertex>=0)
	{
		R3Point position = mesh->VertexPosition(mesh->Vertex(selectedVertex));
		material[0] = 0;		material[1] = 1;		material[2] = 0;	material[3] = 1;
		glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material); 
		glColor3d(0., 0., 1.);
		
		SurfaceFeature * feature = GetFeature("default");		
		if (feature!=NULL && m_sCachedNameMapOnSurf=="Feature")
		{
			R3Point pos = mesh->VertexPosition(mesh->Vertex(selectedVertex));
			m_pWindow->RenderText(pos, VKString::number(feature->Value(SurfaceSample(selectedVertex, mesh))), 1, 0, 0);
			R3Sphere(position, 0.5*GetStandardRadius(params, renderingParams, surfaceName)).Draw();		
		}
		else
			R3Sphere(position, 2.0*GetStandardRadius(params, renderingParams, surfaceName)).Draw();
	}
	DrawSamples(params, renderingParams, surfaceName);
}

SurfaceSampleSet * SampledSurface::GetRenderableSampleSet(ParamParser * drawParams,
														  const VKString & renderingParams,
														  const VKString & surfaceName)
{
	bool genPipeline = drawParams->GetStrValue(renderingParams, "Type", valid)=="GeneralPipeline";
	bool interactCorrs = drawParams->GetStrValue(renderingParams, "InteractionMode", valid)=="ManualCorrespondences";

	SurfaceSampleSet * sampleSet = NULL;
	bool renderSet = m_pWindow->IsFlagSet(renderingParams, "SurfFlags", "DefaultSamples", 
										  surfaceName, "SurfFlagsFlips");
	bool renderIDs = m_pWindow->IsFlagSet(renderingParams, "SurfFlags", "SampleIDs", 
										  surfaceName, "SurfFlagsFlips");
	
	if (!genPipeline && (renderSet || renderIDs))
		sampleSet = GetSampleSet("default");
	else if (genPipeline)
	{
		if (interactCorrs)
			sampleSet = GetSampleSet(AnalysisWindow::s_InteractiveSampleSetName);
		else
			sampleSet = GetSampleSet(drawParams->GetStrValue(renderingParams, "RenderSampleSet", valid));
		renderSet = !renderIDs;
	}
	return sampleSet;
}

void SampledSurface::DrawSamples(ParamParser * params,
								 const VKString & renderingParams,
								 const VKString & surfaceName)
{	
	bool genPipeline = params->GetStrValue(renderingParams, "Type", valid)=="GeneralPipeline";
	bool interactCorrs = params->GetStrValue(renderingParams, "InteractionMode", valid)=="ManualCorrespondences";
	
	static GLfloat material[4];	
	int selectedSample = params->GetIntValue(renderingParams, "SelectedSample", valid);
	if (!valid) selectedSample = -1;
	int selectedVertex = params->GetIntValue(renderingParams, "SelectedVertex", valid);	
	if (!valid) selectedVertex = -1;

	SurfaceDistance * distance=NULL;
	if (m_pWindow->IsFlagSet(renderingParams, "SurfFlags", "DefaultMetric", 
							 surfaceName, "SurfFlagsFlips"))
		distance = GetDistanceMetric("default");
	
	if (distance!=NULL && selectedVertex==-1)
		distance=NULL;
	
	SurfaceSampleSet * sampleSet = NULL;	
	bool renderSet = m_pWindow->IsFlagSet(renderingParams, "SurfFlags", "DefaultSamples", 
										  surfaceName, "SurfFlagsFlips");
	bool renderIDs = m_pWindow->IsFlagSet(renderingParams, "SurfFlags", "SampleIDs", 
										  surfaceName, "SurfFlagsFlips");

	if (!genPipeline && (renderSet || renderIDs))
		sampleSet = GetSampleSet("default");
	else if (genPipeline)
	{
		if (interactCorrs)
			sampleSet = GetSampleSet(AnalysisWindow::s_InteractiveSampleSetName);
		else
			sampleSet = GetSampleSet(params->GetStrValue(renderingParams, "RenderSampleSet", valid));
		renderSet = !renderIDs;
	}
	
	// render samples 
	if (sampleSet!=NULL)
	{
		if (renderIDs)
		{
			glDisable(GL_LIGHTING);
			for (int i=0; i<sampleSet->NumSamples(); i++)
			{
				if (distance==NULL)
					m_pWindow->RenderText(GetDrawablePosition(sampleSet->GetSample(i), 
															  params, renderingParams, surfaceName),
										  VKString(i), 1, 0, 0);
				else
					m_pWindow->RenderText(GetDrawablePosition(sampleSet->GetSample(i), params, renderingParams, surfaceName),
										  distance->Distance(sampleSet->GetSample(i), SurfaceSample(selectedVertex, m_pMesh)), 
										  1, 0, 0);
			}
		}
		else if (renderSet)
		{
			glEnable(GL_LIGHTING);

			double radius = GetStandardRadius(params, renderingParams, surfaceName);
			for (int i=0; i<sampleSet->NumSamples(); i++)
			{
				m_pWindow->MapIntToColor(i, material);
				glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material); 
				if (selectedSample==i)
					R3Sphere(GetDrawablePosition(sampleSet->GetSample(i), params, renderingParams, surfaceName), 
							 2*radius).Draw();
				else
					R3Sphere(GetDrawablePosition(sampleSet->GetSample(i), params, renderingParams, surfaceName), 
							 radius).Draw();
			}	
		}
	}
}

R3Point SampledSurface::GetDrawablePosition(const SurfaceSample & sample, 
											ParamParser * ,
											const VKString & ,
											const VKString & )
{
	assert(sample.CheckMeshIsSame(m_pMesh));
	return m_pWindow->PositionInMyCoords(m_iCameraID, 0, sample.GetPosition());
}

void SampledSurface::ClearCachedSurfaceColors()
{
	m_iCachedVertexFrom = -1;
	m_iCachedConfID = -1;
	m_sCachedNameMapOnSurf="";		
	m_CachedMainColorSet = "";	
	m_pCachedSurfaceFeature=NULL;
	m_pCachedSurfaceDistance=NULL;
	m_pCachedSurfaceMap=NULL;
	m_pCachedSurfaceOtherMap=NULL;
	m_sCachedMapSimilarity="";
	m_sCachedMapConfidence="";
}

void SampledSurface::SetDrawSurfaceColors(ParamParser * params, const VKString & mainSet, SurfaceMap * currMap,
										  SurfaceMap * otherMap, const VKString & similarityName,
										  const VKString & confidenceName, const VKString & renderingParams,
										  const VKString & surfaceName)
{
	if (m_CachedPerVertexColors.size()==0)
		for (int i=0; i<m_pMesh->NVertices()*3; i++)
			m_CachedPerVertexColors.push_back(-1);
	
	bool genPipeline = params->GetStrValue(renderingParams, "Type", valid)=="GeneralPipeline";

	// Features, Distances, SampleSets to visualize
	VKString drawType = params->GetStrValue(renderingParams, "DrawMapColorsOnSurf", valid);
	if (genPipeline)	// try to figure out which signal to render given selection
	{
		drawType = params->GetStrValue(renderingParams, "SignalOnSurface", valid);
		if (drawType=="auto")
		{
			if (params->GetStrValue(renderingParams, "RenderFeatureValue", valid)!="none")
				drawType="Feature";
			else if (params->GetStrValue(renderingParams, "RenderDistanceMetric", valid)!="none")
				drawType="Distance";
			else if (params->GetStrValue(renderingParams, "RenderMap", valid)!="none" && currMap!=NULL)
			{
				if (params->GetStrValue(renderingParams, "RenderOtherMap", valid)!="none"
					&& params->GetStrValue(renderingParams, "RenderSimilarity", valid)!="none"
					&& otherMap!=NULL)
					drawType="MapDiscrepancy";
				else if (params->GetStrValue(renderingParams, "RenderConfidence", valid)!="none")
					drawType="MapConfidence";
				else if (params->GetStrValue(renderingParams, "RenderTexture", valid)!="none")
					drawType="Texture";
				else
					drawType="XYZ";
			}
			else
				drawType="none";
		}
	}
//	std::cout<<"drawType="<<drawType.c_str()<<std::endl;
	bool renderFaces = m_pWindow->IsFlagSet(renderingParams, "MeshFlags", "Faces", surfaceName, "MeshFlagsFlips");
	
	int subMapID = -1;
	if (currMap!=NULL && currMap->GetSurfaceMapType()=="MapScoredCollection")
		subMapID = ((MapScoredCollection*)currMap)->GetSelectedMapID();

	bool cached = (m_sCachedNameMapOnSurf==drawType);
	if (drawType=="none")
	{
		m_sCachedNameMapOnSurf = "none";
		return;
	}
	
	// precompute colors if necessary
	double minValue = 0;
	double maxValue = 1;
	double r, g, b;
	if (renderFaces)
	{
		if (drawType=="Feature")	// Draw feature
		{
			SurfaceFeature * feature = GetFeature("default");
			if (genPipeline)
				feature = GetFeature(params->GetStrValue(renderingParams, "RenderFeatureValue", valid));
			
			if (feature==NULL)
			{
				m_sCachedNameMapOnSurf = "none";
				return;
			}
			
			if (m_pCachedSurfaceFeature==feature && cached)
				return;

			m_pCachedSurfaceFeature = feature;
			m_sCachedNameMapOnSurf = drawType;
			
			minValue = feature->GetMinimalValue();
			maxValue = feature->GetMaximalValue();
			std::cout<<"\tValue in ["<<minValue<<", "<<maxValue<<"]"<<std::endl;
			for (int i=0; i<m_pMesh->NVertices(); i++)
			{
				double val = feature->Value(SurfaceSample(i, m_pMesh));
				m_pWindow->MapScalarToColor(val, minValue, maxValue, m_CachedPerVertexColors[3*i+0], 
											m_CachedPerVertexColors[3*i+1], m_CachedPerVertexColors[3*i+2]);
			}
		}
		else if (drawType=="Distance"
				 || drawType=="ThresholdedDistance")	// draw distance
		{			
			SurfaceDistance * distance = GetDistanceMetric("default");
			if (genPipeline)
				distance = GetDistanceMetric(params->GetStrValue(renderingParams, "RenderDistanceMetric", valid));
			// Selected samples / vertices
			int selectedVertex = GetSelectedVertex(params);
			
			if (distance==NULL || selectedVertex==-1)
			{
				m_sCachedNameMapOnSurf = "none";
				return;
			}
			
			if (distance==m_pCachedSurfaceDistance && cached && m_iCachedVertexFrom==selectedVertex)
				return;
			
			m_pCachedSurfaceDistance = distance;
			m_sCachedNameMapOnSurf = drawType;
			m_iCachedVertexFrom = selectedVertex;
			
			maxValue = 0;
			SurfaceSample fromVert(selectedVertex, GetMesh());
			if (drawType=="Distance")
			{
				for (int i=0; i<GetMesh()->NVertices(); i++)
				{
					m_CachedPerVertexColors[3*i+0] = distance->Distance(fromVert, SurfaceSample(i, GetMesh()));
					maxValue = (maxValue < m_CachedPerVertexColors[3*i+0]) ? m_CachedPerVertexColors[3*i+0] : maxValue;
				}
			}
			else if (drawType=="ThresholdedDistance")
			{
				std::vector<double> thresholds = params->GetDoubleValues(renderingParams, "DistanceThresholds", valid);
				assert(valid && thresholds.size()>0);
				thresholds.insert(thresholds.begin(), 0.);
				maxValue = thresholds[thresholds.size()-1];
				double distanceNorm = sqrt(Area());
				for (int i=0; i<GetMesh()->NVertices(); i++)
				{
					double dist = distance->Distance(fromVert, SurfaceSample(i, GetMesh())) / distanceNorm;
					double clampedValue = maxValue;
					assert(dist >= 0);
					for (int t=0; t<(int)thresholds.size(); t++)
						if (thresholds[t] <= dist)
							clampedValue = thresholds[t];
					m_CachedPerVertexColors[3*i+0] = clampedValue;
				}				
			}
			else
				assert(false);
			
			for (int i=0; i<GetMesh()->NVertices(); i++)
			{
				m_pWindow->MapScalarToColor(m_CachedPerVertexColors[3*i+0], minValue, maxValue, 
											m_CachedPerVertexColors[3*i+0], m_CachedPerVertexColors[3*i+1], 
											m_CachedPerVertexColors[3*i+2]);
			}
		}
		else if (currMap!=NULL && otherMap==NULL 
				 && drawType=="MapConfidence")	// Render map's area preservation
		{
			
			if (currMap==m_pCachedSurfaceMap && cached && confidenceName==m_sCachedMapConfidence &&
				m_iCachedConfID==subMapID)
				return;
			
			if (currMap->GetSurfaceMapType()=="MapScoredCollection")
			{
				MapScoredCollection * collection = (MapScoredCollection*)currMap;
				double val;
				collection->GetMapByID(subMapID, &val);

				std::cout<<"Map In Collection: "<<subMapID<<" Value="<<val<<std::endl;
			}
			
			m_iCachedConfID = subMapID;
			m_pCachedSurfaceMap = currMap;
			m_sCachedNameMapOnSurf = drawType;
			m_sCachedMapConfidence = confidenceName;
			
			SurfaceMapConfidence * confidence = currMap->GetMapConfidenceCalculator(confidenceName);

			//std::cout<<"Precomputing Confidence: "<<currMap<<" cached="<<cached<<" confName=";
			//std::cout<<m_sCachedMapConfidence.c_str()<<" ptr="<<confidence<<" subMapID="<<subMapID<<std::endl;

			assert(confidence!=NULL);
			confidence->Confidence();
			std::vector<double> perVertexVals;
			confidence->GetPerVertexConfidence(this, perVertexVals);
			assert((int)perVertexVals.size()==m_pMesh->NVertices());
			for (int i=0; i<m_pMesh->NVertices(); i++)
			{
				m_CachedPerVertexColors[3*i+0] = 0;
				m_CachedPerVertexColors[3*i+1] = 0;
				m_CachedPerVertexColors[3*i+2] = 0;
				
				double val = 1. - perVertexVals[i];
				if (val>=0 && val<=1)
					m_pWindow->MapScalarToColor(val, 0, 1, m_CachedPerVertexColors[3*i+0], 
												m_CachedPerVertexColors[3*i+1], 
												m_CachedPerVertexColors[3*i+2]);
			}	
		}
		else if (currMap!=NULL && otherMap==NULL && drawType=="ErrorColors")	// render error itself
		{
			if (currMap==m_pCachedSurfaceMap && m_CachedMainColorSet==mainSet && cached)
				return;
			
			m_pCachedSurfaceMap = currMap;
			m_sCachedNameMapOnSurf = drawType;
			m_CachedMainColorSet = mainSet;

			// check that map is coarse map with stored errors
			VKString truthName="";
			std::cout<<"[WARNING] TODO: hand-setting true confidence name. SampledSurface.cpp"<<std::endl;;
			truthName = "MapConfidence_Truth";
			if (truthName=="")
			{
				std::cout<<"[WARNING] Invalid truth name"<<std::endl;
				return;
			}
			
			SurfaceMapConfidence * confidence = currMap->GetMapConfidenceCalculator(truthName);
			if (confidence==NULL)
			{
				std::cout<<"[WARNING] Could not find confidence for the truth"<<std::endl;
				return;
			}
			
			assert(confidence->GetMapConfidenceType()=="MapConfidenceCompareToTruth");
			MapConfidenceCompareToTruth * truthConfidence = (MapConfidenceCompareToTruth*)confidence;
			double threshold = truthConfidence->ErrorThreshold();
			std::vector<double> perVertexError;			
			confidence->GetPerVertexErrors(this, perVertexError, false);
			
			assert(GetMesh()->NVertices()==(int)perVertexError.size());
			
			// color errros
			for (int i=0; i<GetMesh()->NVertices(); i++)
			{
				m_CachedPerVertexColors[3*i+0] = 0; 
				m_CachedPerVertexColors[3*i+1] = 0;
				m_CachedPerVertexColors[3*i+2] = 0;

				double error = perVertexError[i] / threshold;
				if (error>1.)
					error = 1.;
				if (error>=0 && error<=1.)
					m_pWindow->MapScalarToColor(error, 0, 1, 
												m_CachedPerVertexColors[3*i+0], 
												m_CachedPerVertexColors[3*i+1], 
												m_CachedPerVertexColors[3*i+2]);
			}			
		}
		else if (currMap!=NULL && drawType=="Texture")
		{
			if (currMap==m_pCachedSurfaceMap && cached)
				return;

			m_pCachedSurfaceMap = currMap;
			m_sCachedNameMapOnSurf = drawType;
			m_CachedMainColorSet = mainSet;
			
			// save u/v coordinates
			for (int i=0; i<GetMesh()->NVertices(); i++)
			{
				SurfaceSample s(i, GetMesh());
				R3Point p = currMap->ForwardMap(s).GetPosition();
				if (p.Z()==0)
				{
					m_CachedPerVertexColors[3*i + 0] = p.X();
					m_CachedPerVertexColors[3*i + 1] = p.Y();
				}
				else
				{
					std::cout<<"[WARNING] Texture mapping - map is not flattening"<<std::endl;
					break;
				}
			}
		}
		else if (currMap!=NULL && otherMap==NULL 
				 && (drawType=="FarVertices"|| drawType=="XYZ"))	// render map itself
		{
			if (currMap==m_pCachedSurfaceMap && m_CachedMainColorSet==mainSet && cached)
				return;
			
			bool symmetry = currMap->GetSurface(0)==currMap->GetSurface(1);
			
			m_pCachedSurfaceMap = currMap;
			m_sCachedNameMapOnSurf = drawType;
			m_CachedMainColorSet = mainSet;
			
			if (!m_bCachedBBox && drawType=="XYZ")
			{
				m_fCachedMinX=FLT_MAX, m_fCachedMinY=FLT_MAX, m_fCachedMinZ=FLT_MAX;
				m_fCachedMaxX=-FLT_MAX, m_fCachedMaxY=-FLT_MAX, m_fCachedMaxZ=-FLT_MAX;
				
				R3Mesh * msh = currMap->GetSurface(1)->m_pMesh;
				for (int i=0; i<msh->NVertices(); i++)
				{
					R3MeshVertex * vrt = msh->Vertex(i);
					m_fCachedMinX = vkMin(m_fCachedMinX, msh->VertexPosition(vrt).X());
					m_fCachedMinY = vkMin(m_fCachedMinY, msh->VertexPosition(vrt).Y());
					m_fCachedMinZ = vkMin(m_fCachedMinZ, msh->VertexPosition(vrt).Z());
					m_fCachedMaxX = vkMax(m_fCachedMaxX, msh->VertexPosition(vrt).X());
					m_fCachedMaxY = vkMax(m_fCachedMaxY, msh->VertexPosition(vrt).Y());
					m_fCachedMaxZ = vkMax(m_fCachedMaxZ, msh->VertexPosition(vrt).Z());
				}
				m_bCachedBBox = true;
			}
			

			std::cout<<"Precomputing Colors "<<m_pMesh->NVertices()<<std::endl;
			//std::cout<<"[WARNING] Precomputing Colors ("<<m_pMesh->NVertices()<<"): "<<std::flush;
			TimeProfiler profiler;
			for (int vertexID=0; vertexID<m_pMesh->NVertices(); vertexID++)
			{
				if (vertexID%100==0)
					profiler.WriteProgress("PrecomputingColors", vertexID, m_pMesh->NVertices());
				
				SurfaceSample colorSample;
				SurfaceSampleSet * colorSet=NULL;
				SurfaceDistance * distanceMetric = NULL;
				double norm=0;
				
				if (currMap->GetSurface(1)==this && !symmetry)		// interpolate colors
				{
					colorSample = SurfaceSample(vertexID, m_pMesh);
					colorSet = GetSampleSet(m_CachedMainColorSet);
					distanceMetric = GetDistanceMetric("default");
				}
				else
				{
					assert(currMap->GetSurface(0)==this);
					colorSample = currMap->ForwardMap(SurfaceSample(vertexID, m_pMesh));
					colorSet = currMap->GetSurface(1)->GetSampleSet(mainSet);
					distanceMetric = currMap->GetSurface(1)->GetDistanceMetric("default");
				}	

				r=1;	g=0;	b=0;				
				if (!colorSample.Invalid())
				{
					r=0;	g=0;	b=0;
					if (drawType=="FarVertices")
					{
						assert(distanceMetric!=NULL);
						assert(colorSet!=NULL);
						for (int i=0; i<colorSet->NumSamples(); i++)
						{
							assert(colorSample.CheckMeshIsSame(currMap->GetSurface(1)->GetMesh()));
							assert(colorSet->GetSample(i).CheckMeshIsSame(currMap->GetSurface(1)->GetMesh()));								
							double weight = 1 / (0.000001 + distanceMetric->Distance(colorSample, colorSet->GetSample(i)));
							double tr, tg, tb;
							m_pWindow->MapIntToRadicalColor(i, tr, tg, tb);
							r += tr * weight;	g += tg * weight;	b += tb * weight;
							norm += weight;
						}
						r /= norm;		g /= norm;		b /= norm;
					}
					else if (drawType=="XYZ")
					{
						R3Point pnt = colorSample.GetPosition();						
						r = (pnt.X()-m_fCachedMinX) / (m_fCachedMaxX-m_fCachedMinX);
						g = (pnt.Y()-m_fCachedMinY) / (m_fCachedMaxY-m_fCachedMinY);
						b = (pnt.Z()-m_fCachedMinZ) / (m_fCachedMaxZ-m_fCachedMinZ);
						
						if (symmetry)
						{
							SurfaceSample colorSample2(vertexID, m_pMesh);
							R3Point pnt2 = colorSample2.GetPosition();
							r += (pnt2.X()-m_fCachedMinX) / (m_fCachedMaxX-m_fCachedMinX);
							g += (pnt2.Y()-m_fCachedMinY) / (m_fCachedMaxY-m_fCachedMinY);
							b += (pnt2.Z()-m_fCachedMinZ) / (m_fCachedMaxZ-m_fCachedMinZ);
							r /= 2.;
							g /= 2.;
							b /= 2. ;
						}
						
					}
					else
						assert(false);
				}
				m_CachedPerVertexColors[vertexID*3+0] = r;
				m_CachedPerVertexColors[vertexID*3+1] = g;
				m_CachedPerVertexColors[vertexID*3+2] = b;				
			}
			std::cout<<std::endl;
		}
		else if (currMap!=NULL && otherMap!=NULL && similarityName!="none"
				 && drawType=="MapDiscrepancy")		// Render discrepancy between two maps
		{
			if (currMap==m_pCachedSurfaceMap && otherMap==m_pCachedSurfaceOtherMap 
				&& similarityName==m_sCachedMapSimilarity && cached)
				return;
			
			if (currMap->GetSurface(1)==this)	// don't render similarity on the other surface
				return;
			
			m_pCachedSurfaceMap = currMap;
			m_pCachedSurfaceOtherMap = otherMap;
			m_sCachedMapSimilarity = similarityName;
			m_sCachedNameMapOnSurf = drawType;
			
			SurfaceMapSimilarity * similarity = currMap->GetSurfaceMapSimilarity(similarityName, otherMap);
			assert(similarity!=NULL);
			std::vector<double> perVertexSimVals;
			std::cout<<"[WARNING] Calculating per vertex map discrepancy "<<std::flush;
			//similarity->GetPerVertexDissimilarity(this, perVertexSimVals, false);
			similarity->GetPerVertexSimilarity(this, perVertexSimVals, false);
			for (int i=0; i<(int)perVertexSimVals.size(); i++)
				perVertexSimVals[i] = 1. - perVertexSimVals[i];
			std::cout<<" - done!"<<std::endl;			
			for (int i=0; i<(int)perVertexSimVals.size(); i++)
			{
//				std::cout<<"val["<<i<<"] = "<<perVertexSimVals[i]<<std::endl;
				m_pWindow->MapScalarToColor(perVertexSimVals[i], 0, 1, m_CachedPerVertexColors[3*i+0], 
											m_CachedPerVertexColors[3*i+1], m_CachedPerVertexColors[3*i+2]);
			}
		}
		else if (currMap!=NULL && (drawType=="MultiMapDistGen" 
								   || drawType=="MultiMapConsistency" 
								   || drawType=="MultiMapConfidence"
								   || drawType=="MultiMapFinal"))
		{			
			if (VKString(currMap->GetSurfaceMapType().c_str()) != "MapMultiConformal")
			{
				std::cout<<"[WARNING] Cannot render "<<drawType.c_str()<<" map type=";
				std::cout<<currMap->GetSurfaceMapType().c_str()<<std::endl;
				
				return;
			}
					
			MapMultiConformal * multiConfMap = (MapMultiConformal*)currMap;
			int confMapID = params->GetIntValue(renderingParams, "ConfMapInMultiMap", valid);

			if (!valid || confMapID==-1)
			{
				std::cout<<"[WARNING] Cannot render "<<drawType.c_str();
				std::cout<<". Conformal map not selected"<<std::endl;
				return;				
			}
		
			confMapID = confMapID % multiConfMap->GetNumConformalMaps();			
			if (currMap==m_pCachedSurfaceMap && cached && m_iCachedConfID==confMapID)
				return;			
			
			m_pCachedSurfaceMap = currMap;
			m_sCachedNameMapOnSurf = drawType;
			m_iCachedConfID = confMapID;

			std::vector<double> perVertexVals;	
			multiConfMap->FillPerVertexWeights(confMapID, drawType, this, perVertexVals);
			
			for (int i=0; i<(int)perVertexVals.size(); i++)
				m_pWindow->MapScalarToColor(1.-perVertexVals[i], 0, 1, m_CachedPerVertexColors[3*i+0], 
											m_CachedPerVertexColors[3*i+1], m_CachedPerVertexColors[3*i+2]);
			
		}
		else if (currMap!=NULL && (drawType=="MultiMapTop3Weights"))
		{
			if (VKString(currMap->GetSurfaceMapType().c_str()) != "MapMultiConformal")
			{
				std::cout<<"[WARNING] Cannot render "<<drawType.c_str()<<" map type=";
				std::cout<<currMap->GetSurfaceMapType().c_str()<<std::endl;
				
				return;
			}
			MapMultiConformal * multiConfMap = (MapMultiConformal*)currMap;		
			assert(multiConfMap->GetNumConformalMaps()>=3);
			
			if (currMap==m_pCachedSurfaceMap && cached)
				return;			
			
			m_pCachedSurfaceMap = currMap;
			m_sCachedNameMapOnSurf = drawType;
			
			std::vector<double> perVertexValsR;	
			std::vector<double> perVertexValsG;	
			std::vector<double> perVertexValsB;				
			multiConfMap->FillPerVertexWeights(0, "MultiMapFinal", this, perVertexValsR);
			multiConfMap->FillPerVertexWeights(1, "MultiMapFinal", this, perVertexValsG);
			multiConfMap->FillPerVertexWeights(2, "MultiMapFinal", this, perVertexValsB);			
			
			for (int i=0; i<(int)perVertexValsR.size(); i++)
			{
				double norm = perVertexValsR[i] + perVertexValsG[i] + perVertexValsB[i];
				m_CachedPerVertexColors[3*i+0] = perVertexValsR[i] / norm;
				m_CachedPerVertexColors[3*i+1] = perVertexValsG[i] / norm;
				m_CachedPerVertexColors[3*i+2] = perVertexValsB[i] / norm;				
			}			
		}	
	}
	else
		m_sCachedNameMapOnSurf = "none";
}

int SampledSurface::GetSelectedVertex(ParamParser * drawParams)
{
	std::vector<int> selectedVertexDescription = drawParams->GetIntValues("RendererDefault", "SelectedVertex", valid);
	for (int i=0; i<(int)selectedVertexDescription.size(); i+=2)
		if (selectedVertexDescription[i]==m_iCameraID)
			return selectedVertexDescription[i+1];
	return -1;
}

///////////// CHECK / FIX MESH  /////////////////
R3Mesh * SampledSurface::CheckMeshForProblems(R3Mesh * mesh)
{
	int numIsolated=MeshProcessor::NumIsolatedVertices(mesh);
	int numConnected=MeshProcessor::NumConnectedComponents(mesh);	
	if (numIsolated>0)
		std::cout<<"[ERROR] Num Isolated Vertices = "<<numIsolated<<std::endl;
	if (numConnected!=1)
		std::cout<<"[ERROR] Num Connected Components = "<<numConnected<<std::endl;
	assert(numIsolated==0 && numConnected==1);
	return mesh;
}

//////////////// DUMMY SURFACE ////////////////////
R3Mesh * Surface2DPlane::m_pInfinitePlanePseudomesh = new R3Mesh();
SampledSurface * Surface2DPlane::m_pInfinitePlanePseudosurface = new SampledSurface(Surface2DPlane::m_pInfinitePlanePseudomesh);
