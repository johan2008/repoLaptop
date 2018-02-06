#include "SurfaceMapSimilarity.h"
#include "MapConfidenceMultiConf.h"
#include "MapMultiConformal.h"	

MapConfidenceMultiConf::MapConfidenceMultiConf(SurfaceMap * currMap, 
											   const VKString & setNameFrom,
											   const VKString & setNameTo,
											   const VKStringList & confidences, 
											   const VKString & similarity,
											   WeightCalculation weightCalculation)
: SurfaceMapConfidence(currMap, false, false, -1., setNameFrom, setNameTo, true, false)

{
	m_fIfNoGoodConformalMapTheshold = 0;
	m_pSimilarityMatrix = NULL;
	m_pEigenVectors = NULL;
	m_pEigenValues = NULL;
	
	m_bWeightsReady = false;
	m_bInterpolatedWeights = true;
	m_WeightCalculation = weightCalculation;
	//m_bPresmoothSimConf = m_WeightCalculation==MULTICONF_WEIGHT_CONFIDENCE_AND_CONSISTENCY;
	//m_bPresmoothSimConf = m_bPresmoothSimConf || m_WeightCalculation==MULTICONF_WEIGHT_CONFIDENCE;
	m_bPresmoothSimConf = false;
	
	m_bCacheWeights = true;
	m_Confidences = confidences;
	m_Similarity = similarity;
	
	m_bDeleteMultiConformal = false;	
	m_pMultiConformalMap = NULL;

	if (currMap!=NULL)
	{		
		if (currMap->GetSurfaceMapType()=="MapMultiConformal")
			m_pMultiConformalMap = (MapMultiConformal*)currMap;
		else if (currMap->GetSurfaceMapType()=="MapConformal")
		{
			m_bDeleteMultiConformal = true;
			// TODO: Get self's name (instead of "none")
			m_pMultiConformalMap = new MapMultiConformal((MapConformal*)currMap, "none");
		}	
		else
			assert(false);
	}
}

MapConfidenceMultiConf::WeightCalculation 
						MapConfidenceMultiConf::StrToWeightCalculation(const VKString & str)
{	
	if (str=="MULTICONF_WEIGHT_CONFIDENCE")
		return MULTICONF_WEIGHT_CONFIDENCE;
	else if (str=="MULTICONF_WEIGHT_CONFIDENCE_AND_CONSISTENCY")
		return MULTICONF_WEIGHT_CONFIDENCE_AND_CONSISTENCY;
	else if (str=="MULTICONF_WEIGHT_PER_POINT_EIGENANALYSIS")
		return MULTICONF_WEIGHT_PER_POINT_EIGENANALYSIS;
	else if (str=="MULTICONF_WEIGHT_PER_POINT_EIGENANALYSIS_CONFIDENCE")
		return MULTICONF_WEIGHT_PER_POINT_EIGENANALYSIS_CONFIDENCE;
	else if (str=="MULTICONF_WEIGHT_GLOBAL_EIGENANALYSIS")
		return MULTICONF_WEIGHT_GLOBAL_EIGENANALYSIS;
	else 
		assert(false);
}

MapConfidenceMultiConf::~MapConfidenceMultiConf()
{
//	int N = (int)m_pMultiConformalMap->m_ConformalMaps.size();	
	
	if (m_bDeleteMultiConformal)
		delete m_pMultiConformalMap;
	
	if (m_pSimilarityMatrix!=NULL)
		delete m_pSimilarityMatrix;
	
	if (m_pEigenVectors!=NULL)
		delete m_pEigenVectors;
	
	if (m_pEigenValues!=NULL)
		delete m_pEigenValues;
	
	// TODO: memory leak 
//	for (int surf=0; surf<2; surf++)
//	{
//		if (m_SuggestedWeights.find(m_pSurfaceMap->GetSurface(surf))!=m_SuggestedWeights.end());
//		{
//			assert((int)m_SuggestedWeights[m_pSurfaceMap->GetSurface(surf)].size()==N);
//			for (int i=0; i<N; i++)
//				delete m_SuggestedWeights[m_pSurfaceMap->GetSurface(surf)][i];
//		}
//	}

}

void MapConfidenceMultiConf::SetIfNoGoodConformalMapThreshold(double threshold)
{
	m_fIfNoGoodConformalMapTheshold = threshold;
}

void MapConfidenceMultiConf::FillFinalWeights(std::vector< std::map<int, double> > & weights)
{
	ClearMultiConfWeights();
	
	SampledSurface * surf = m_pSurfaceMap->GetSurface(0);
	
	int N = (int)weights.size();
	if (m_SuggestedWeights.find(surf)==m_SuggestedWeights.end())
	{
		for (int i=0; i<N; i++)
			m_SuggestedWeights[surf].push_back(new SurfacePerVertexValue(surf, true));
	}	
	assert((int)m_SuggestedWeights[surf].size()==N);
	for (int i=0; i<N; i++)
	{
		SurfacePerVertexValue * perVertWeight = m_SuggestedWeights[surf][i];
		
		for (std::map<int, double>::iterator iter = weights[i].begin(); 
			 iter != weights[i].end(); iter++)
			perVertWeight->SetVertexValue(iter->first, iter->second);
		perVertWeight->SmoothGaussianFromKnown(-1);
	}
	
//	for (int v=0; v<surf->GetMesh()->NVertices(); v++)
//	{
//		std::vector<double> & perVertexW = m_mVertexToOptimalWeights[surf][v];
//		SurfaceSample vertSamp(v, surf->GetMesh());
//		for (int i=0; i<N; i++)
//			perVertexW.push_back(m_SuggestedWeights[surf][i]->GetInterpolatedValue(vertSamp));
//	}
	m_bWeightsReady = true;
//	m_bInterpolatedWeights = true;
}

void MapConfidenceMultiConf::ClearMultiConfWeights()
{
	m_bWeightsReady = false;
	m_bInterpolatedWeights = false;
	assert(m_pMultiConformalMap!=NULL);
	int N = (int)m_pMultiConformalMap->m_ConformalMaps.size();	
	for (int surf=0; surf<2; surf++)
	for (int i=0; i<N; i++)
	{
		if (m_SuggestedWeights.find(m_pSurfaceMap->GetSurface(surf))!=m_SuggestedWeights.end())
		{
			assert ((int)m_SuggestedWeights[m_pSurfaceMap->GetSurface(surf)].size()==N);
			m_SuggestedWeights[m_pSurfaceMap->GetSurface(surf)][i]->ClearValues();
		}
	}
	m_mVertexToOptimalWeights.clear();
}

bool MapConfidenceMultiConf::WeightsReady()
{
	return m_bWeightsReady;
}

VKString MapConfidenceMultiConf::GetMapConfidenceType()
{
	return "MapConfidenceMultiConf";
}

SurfaceMapConfidence * MapConfidenceMultiConf::Clone(SurfaceMap * currMap)
{
	return new MapConfidenceMultiConf(currMap, m_SampleSet1, m_SampleSet2,
									  m_Confidences, m_Similarity, m_WeightCalculation);
}

void MapConfidenceMultiConf::PrepareFinalWeights()
{
	if (m_bWeightsReady)
		return;
	m_bWeightsReady = true;
	int N = (int)m_pMultiConformalMap->m_ConformalMaps.size();	
	if (m_bPresmoothSimConf)
	{
		assert(false);
		assert(m_WeightCalculation==MULTICONF_WEIGHT_CONFIDENCE
			   || m_WeightCalculation==MULTICONF_WEIGHT_CONFIDENCE_AND_CONSISTENCY);
		if (m_Similarity=="none")
			return;
		
		// for each similarity - smooth
		std::cout<<"Preparing Smooth Weights: "<<std::endl;
		std::cout<<"\tSmoothing Similarities ("<<N<<"): "<<std::flush;		
		for (int i=0; i<N; i++)
		{
			std::cout<<i<<" "<<std::flush;
			MapConformal * m_i = m_pMultiConformalMap->m_ConformalMaps[i];
			for (int j=0; j<N; j++)
			{
				MapConformal * m_j = m_pMultiConformalMap->m_ConformalMaps[j];
				if (m_i!=m_j)
					m_i->GetSurfaceMapSimilarity(m_Similarity, m_j)->PrecomputeSmoothDissimilarities();
			}
		}
		std::cout<<std::endl;
		
		for (int surfaceID = 0; surfaceID < 2; surfaceID++)
		{
			SampledSurface * surface = m_pMultiConformalMap->GetSurface(surfaceID);
			
			std::cout<<"\tFinal Weights {surf"<<surfaceID<<"}";
			std::cout<<" ("<<surface->GetMesh()->NVertices()<<"): "<<std::flush;
			if (m_SuggestedWeights.find(surface)==m_SuggestedWeights.end())
			{
				for (int mapID = 0; mapID < N; mapID++)
				{
					SurfacePerVertexValue * perVertexValue = new SurfacePerVertexValue(surface, true);
					m_SuggestedWeights[surface].push_back(perVertexValue);
				}
			}
			else
				assert((int)m_SuggestedWeights[surface].size()==N);
			
			for (int vID=0; vID<surface->GetMesh()->NVertices(); vID++)
			{
				if (vID%1000==0)
					std::cout<<vID<<" "<<std::flush;
				std::vector<double> weights;
				SurfaceSample atVertexSample = SurfaceSample(vID, surface->GetMesh());
				CalculateWeightsAtSample(weights, atVertexSample);
				assert((int)weights.size() == N);
				for (int mapID=0; mapID<N; mapID++)
					m_SuggestedWeights[surface][mapID]->SetVertexValue(atVertexSample, weights[mapID]);
			}
			std::cout<<std::endl;
		}
	}
	else	// No smoothing
	{
		std::cout<<"Calculating Optimal Weights. Method="<<m_WeightCalculation<<std::endl;
		AnalysisStats::m_GlobalStats.m_Timing.startedProcess("Extrapolation_OptimalWeights");
		int MAX_OPTIM=0;
		for (int surfaceID = 0; surfaceID < 2; surfaceID++)
		{
			SampledSurface * surface = m_pMultiConformalMap->GetSurface(surfaceID);
			SurfaceSampleSet * sampleSet = surface->GetSampleSet(m_SampleSet1);
			MAX_OPTIM += sampleSet->NumSamples();
		}

		for (int surfaceID = 0; surfaceID < 2; surfaceID++)
		{
			SampledSurface * surface = m_pMultiConformalMap->GetSurface(surfaceID);
//			std::cout<<"SurfaceID="<<surfaceID<<" surf="<<surface<<std::endl;
			if (m_SuggestedWeights.find(surface)==m_SuggestedWeights.end())
			{			
				for (int mapID = 0; mapID < N; mapID++)
				{
					SurfacePerVertexValue * perVertexValue = new SurfacePerVertexValue(surface, true);
					m_SuggestedWeights[surface].push_back(perVertexValue);
				}
			}
			else
				assert((int)m_SuggestedWeights[surface].size()==N);
			
			switch(m_WeightCalculation)
			{
				case MULTICONF_WEIGHT_CONFIDENCE:
				case MULTICONF_WEIGHT_CONFIDENCE_AND_CONSISTENCY:
				case MULTICONF_WEIGHT_PER_POINT_EIGENANALYSIS:
				case MULTICONF_WEIGHT_PER_POINT_EIGENANALYSIS_CONFIDENCE:
				{
					SurfaceSampleSet * sampleSet = surface->GetSampleSet(m_SampleSet1);
					for (int sampleID=0; sampleID<sampleSet->NumSamples(); sampleID++)
					{
						if (sampleID%10==0)
						{
							int prog = surfaceID * sampleSet->NumSamples() + sampleID;
							static VKString str ="Extrapolation_OptimalWeights";
							AnalysisStats::m_GlobalStats.m_Timing.WriteProgress(str, prog, MAX_OPTIM);
						}
						const SurfaceSample & sample = sampleSet->GetSample(sampleID);
						bool atVertex;
						sample.NearestVertex(&atVertex);
						assert(atVertex);
						std::vector<double> weights;
						CalculateWeightsAtSample(weights, sample);
						// condition: ad-hoc stuff for extrapolation: see m_fIfNoGoodConformalMapTheshold
						if (weights[0]>=0)
						{
							for (int mapID=0; mapID<N; mapID++)
								m_SuggestedWeights[surface][mapID]->SetVertexValue(sample, weights[mapID]);
						}
					}
					
				}
					break;
				case MULTICONF_WEIGHT_GLOBAL_EIGENANALYSIS:
					FillWeightsUsingGlobalEigenAnalysis(surface);
					break;
				default:
					assert(false);
			}
		}
		std::cout<<std::endl;			
		AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("Extrapolation_OptimalWeights");		
	}
}

void MapConfidenceMultiConf::SmoothFinalWeights()
{
	if (m_bInterpolatedWeights)
		return;
	m_bInterpolatedWeights = true;
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("Extrapolation_SmoothWeights");
	int N = (int)m_pMultiConformalMap->m_ConformalMaps.size();		
	bool needToDoAnything = false;
	int MAX_TO_SMOOTH = N * 2;
	for (int surfaceID = 0; surfaceID < 2; surfaceID++)
	{
		SampledSurface * surface = m_pMultiConformalMap->GetSurface(surfaceID);
		if (m_SuggestedWeights.find(surface)==m_SuggestedWeights.end())	// map from surface not defined
			continue;
		for (int mapID=0; mapID < N ; mapID++)
		{
			assert((int)m_SuggestedWeights[surface].size() > mapID);
			if (!m_SuggestedWeights[surface][mapID]->AllValuesSet())
			{
				if (!needToDoAnything)
				{
					std::cout<<"Smoothing Final Weights ("<<N<<"): "<<std::endl;	
					needToDoAnything = true;
				}
				if (mapID%2==0)
					AnalysisStats::m_GlobalStats.m_Timing.WriteProgress("Extrapolation_SmoothWeights",
																		surfaceID*N + mapID, MAX_TO_SMOOTH);
				const VKString & sampleSetName = (surface == m_pSurfaceMap->GetSurface(0)) ? m_SampleSet1 : m_SampleSet2;
				//m_SuggestedWeights[surface][mapID]->SmoothGaussianFromKnown(sampleSetName, -1);	
				m_SuggestedWeights[surface][mapID]->SmoothGaussianFromKnown(-1);	
			}
		}
	}
	if (needToDoAnything)
		std::cout<<std::endl;	
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("Extrapolation_SmoothWeights");	
}

double MapConfidenceMultiConf::Weight(int mapID, const SurfaceSample & sample)
{	
	SampledSurface * surface = m_pMultiConformalMap->GetSurface(sample);
	
	assert (m_SuggestedWeights.find(surface)!=m_SuggestedWeights.end()
			&& (int)m_SuggestedWeights[surface].size()>mapID);
	return m_SuggestedWeights[surface][mapID]->GetInterpolatedValue(sample);
}

void MapConfidenceMultiConf::CalculateWeightsAtSample(std::vector<double> & weights, 
													  const SurfaceSample & sample)
{
	switch (m_WeightCalculation)
	{
		case MULTICONF_WEIGHT_CONFIDENCE:
			CalcWeightsConfidence(weights, sample);
			break;
		case MULTICONF_WEIGHT_CONFIDENCE_AND_CONSISTENCY:
			CalcWeightsConfidenceConsistency(weights, sample);
			break;
		case MULTICONF_WEIGHT_PER_POINT_EIGENANALYSIS:
		case MULTICONF_WEIGHT_PER_POINT_EIGENANALYSIS_CONFIDENCE:
			CalcWeightsLocalEigenAnalysis(weights, sample);
			break;
		default:
			assert(false);
	}
}



void MapConfidenceMultiConf::FillWeightsUsingGlobalEigenAnalysis(SampledSurface * surface)
{
	// DO NOT FORGET TO NORMALIZE
	assert(false);
}

void MapConfidenceMultiConf::CalcWeightsConfidenceTakeBestConf(std::vector<double> & weights, 
															   const SurfaceSample & sample)
{
	weights.clear();
	int bestID = -1;
	double bestVal=-FLT_MAX;
	int N = (int)m_pMultiConformalMap->m_ConformalMaps.size();
	for (int i=0; i<N; i++)
	{
		double currVal = GetTotalConfidence(i, sample);
		if (currVal > bestVal)
		{
			bestVal = currVal;
			bestID = i;
		}
		weights.push_back(0.);
	}
	assert(bestID>=0);
	weights[bestID] = 1.;
}

void MapConfidenceMultiConf::CalcWeightsConfidence(std::vector<double> & weights, 
												  const SurfaceSample & sample)
{
	int N = (int)m_pMultiConformalMap->m_ConformalMaps.size();
	double maxConfidenceValue = 0;
	for (int i=0; i<N; i++)
	{
		weights.push_back(GetTotalConfidence(i, sample));
		maxConfidenceValue = vkMax(weights[i], maxConfidenceValue);
		
	}
//	if (sample.NearestVertex()==5342
//		|| sample.NearestVertex()==7127)
//		std::cout<<"Vertex["<<sample.NearestVertex()<<"] maxConfVal="<<maxConfidenceValue<<std::endl;
	
	if (maxConfidenceValue < m_fIfNoGoodConformalMapTheshold)
	{
		weights[0] = -1;
	}
	else
		NormalizeWeights(weights);
}

void MapConfidenceMultiConf::CalcWeightsConstant(std::vector<double> & weights)
{
	int N = (int)m_pMultiConformalMap->m_ConformalMaps.size();
	for (int i=0; i<N; i++)
		weights.push_back(1./(double)N);
}

void MapConfidenceMultiConf::CalcWeightsConfidenceConsistency(std::vector<double> & weights, 
															  const SurfaceSample & sample)
{
	int N = (int)m_pMultiConformalMap->m_ConformalMaps.size();
	for (int i=0; i<N; i++)
		weights.push_back(GetTotalConfidence(i, sample) * GetConsistency(i, sample));
	NormalizeWeights(weights);
}

void MapConfidenceMultiConf::CalcWeightsLocalEigenAnalysis(std::vector<double> & weights, 
														  const SurfaceSample & sample)
{	
	weights.clear();
	// Initialize datastructures
	int N = (int)m_pMultiConformalMap->m_ConformalMaps.size();
	if (m_pSimilarityMatrix==NULL)
		m_pSimilarityMatrix = new LinAlgMatrixReal(N, N);
	if (m_pEigenVectors==NULL)
		m_pEigenVectors = new LinAlgMatrixReal(N, N);
	if (m_pEigenValues==NULL)
		m_pEigenValues = new LinAlgVectorReal(N);	
	
	// Fill matrix m_pSimilarityMatrix = S_ij c_i w_i c_j w_j
	for (int i=0; i<N; i++)
	{
		double c_i = GetTotalConfidence(i, sample);		
//		m_pSimilarityMatrix->AddVal(i, i, c_i*c_i); 		
		(*m_pSimilarityMatrix)(i, i) = c_i*c_i;
//		weights.push_back(c_i);
		for (int j=i+1; j<N; j++)
		{
			double c_j = GetTotalConfidence(j, sample);		
			double S_i_j = GetSimilarity(i, j, sample);
//			m_pSimilarityMatrix->AddVal(i, j, S_i_j*c_i*c_j);	// ENERGY:
//			m_pSimilarityMatrix->AddVal(j, i, S_i_j*c_i*c_j);	
			(*m_pSimilarityMatrix)(i, j) = S_i_j*c_i*c_j;
			(*m_pSimilarityMatrix)(j, i) = S_i_j*c_i*c_j;			
		}
	}
	
	// Eigendecomposition of A: 
//	m_pSimilarityMatrix->RowNormalize();
	m_pSimilarityMatrix->EigendecompositionSymmetric(*m_pEigenValues, *m_pEigenVectors);
	
	// Use top eigenvector
	for (int i=0; i<N; i++)
	{
		weights.push_back((*m_pEigenVectors)(i, 0));
		if (weights[i] < 0)
			weights[i] *= -1;	// top eigenvector can have all negative entries
	}

	NormalizeWeights2(weights);
	SampledSurface * surf = m_pMultiConformalMap->GetSurface(sample);;	
	bool atVertex = false;
	int vertex = sample.NearestVertex(&atVertex);
	assert(atVertex);
	for (int i=0; i<N; i++)
		m_mVertexToOptimalWeights[surf][vertex].push_back(weights[i]);

	if (m_WeightCalculation==MULTICONF_WEIGHT_PER_POINT_EIGENANALYSIS_CONFIDENCE)
	{
		for (int i=0; i<N; i++)
			weights[i] *= GetTotalConfidence(i, sample);
	}	
	NormalizeWeights(weights);
}

double MapConfidenceMultiConf::GetTotalConfidence(int mapID, const SurfaceSample & sample)
{
	// introduce global scale
	double globalScale = 1.;
	if ((int)m_pMultiConformalMap->m_MapWeight.size() > mapID)
		globalScale = m_pMultiConformalMap->m_MapWeight[mapID];
	
	double confidence = 1.;	// combined confidence
	for (int cID = 0; cID < m_Confidences.count(); cID++)
	{
		double val = GetConfidence(mapID, cID, sample);
		confidence *= val;
	}
	return confidence * globalScale;
}

double MapConfidenceMultiConf::GetConfidence(int mapID, int confID, const SurfaceSample & sample)
{	
	if (confID >= m_Confidences.count() || m_Confidences[confID]=="none")
		return 1.;
	else
	{
		MapConformal * m_j = m_pMultiConformalMap->m_ConformalMaps[mapID];
		return m_j->GetMapConfidenceCalculator(m_Confidences[confID])->ConfidenceAtSampleNoSmoothing(sample);
	}
}

double MapConfidenceMultiConf::GetConsistency(int mapID, const SurfaceSample & sample)
{	
	int N = (int)m_pMultiConformalMap->m_ConformalMaps.size();
	
	double consistency = 0;
	for (int j=0; j<N; j++)
	{
		double c_j = GetTotalConfidence(mapID, sample);		
		double S_i_j = GetSimilarity(mapID, j, sample);
		
//		consistency += sqrt(S_i_j * c_j);
		consistency += S_i_j * c_j;
	}
	return consistency;
}

double MapConfidenceMultiConf::GetSimilarity(int i, int j, const SurfaceSample & sample)
{
	if (i==j || m_Similarity=="none")
		return 1.;
	MapConformal * m_i = m_pMultiConformalMap->m_ConformalMaps[i];
	MapConformal * m_j = m_pMultiConformalMap->m_ConformalMaps[j];	
	return m_i->GetSurfaceMapSimilarity(m_Similarity, m_j)->SimilarityAtSample(sample);
}

void MapConfidenceMultiConf::NormalizeWeights(std::vector<double> & weights)
{
	double norm = 0;
	for (int i=0; i<(int)weights.size(); i++)
		norm += weights[i];
	norm = norm;
	if (norm!=0)
	for (int i=0; i<(int)weights.size(); i++)
		weights[i] /= norm;
}

void MapConfidenceMultiConf::NormalizeWeights2(std::vector<double> & weights)
{	
	double norm = 0;	
	for (int i=0; i<(int)weights.size(); i++)
		norm += weights[i]*weights[i];
	norm = sqrt(norm);
	if (norm!=0)
		for (int i=0; i<(int)weights.size(); i++)
			weights[i] /= norm;
}

double MapConfidenceMultiConf::CalculateConfidenceAtSample(const SurfaceSample & sample)
{
	double confidenceVal=0;
	int N = (int)m_pMultiConformalMap->m_ConformalMaps.size();	
	
	// approximate weights
	std::vector<double> weights;
	if (m_mVertexToOptimalWeights.size()==0)
	{
		//CalcWeightsConstant(weights);
		CalcWeightsConfidence(weights, sample);
		NormalizeWeights2(weights);		
	}
	else
	{
		bool atVertex;
		SampledSurface * surf = m_pMultiConformalMap->GetSurface(sample);
		assert(m_mVertexToOptimalWeights.find(surf)!=m_mVertexToOptimalWeights.end());
		std::map<int, std::vector<double> > & vertexToWeight = m_mVertexToOptimalWeights[surf];
		int vertID = sample.NearestVertex(&atVertex);
		std::map<int, std::vector<double> >::iterator iter = vertexToWeight.find(vertID);
		assert(iter!=vertexToWeight.end() && atVertex);
		assert((int)iter->second.size()==N);
		for (int i=0; i<N; i++)
			weights.push_back(iter->second[i]);
	}
	
	// calculate energy
	for (int i=0; i<N; i++)
	{
		double c_i = GetTotalConfidence(i, sample);		
		for (int j=0; j<N; j++)
		{
			double c_j = GetTotalConfidence(j, sample);		
			double S_i_j = GetSimilarity(i, j, sample);
	//		confidenceVal += sqrt(S_i_j*c_i*c_j * weights[i] * weights[j]);
			confidenceVal += S_i_j*c_i*c_j * weights[i] * weights[j];	// objective function:
		}
	}	
	return confidenceVal;
}

void MapConfidenceMultiConf::SaveConfidence(std::ofstream & textStream)
{	
	int N = (int)m_pMultiConformalMap->m_ConformalMaps.size();
	
	// save values
	SaveConfidenceValues(textStream);
	
	// save names
	textStream<<"Confidences "<<m_Confidences.count()<<"\n";
	for (int i=0; i<m_Confidences.count(); i++)
		textStream<<m_Confidences[i].c_str()<<" ";
	
	textStream<<"\nSimilarity "<<m_Similarity.c_str()<<"\n";
	
	// save weights
	textStream<<"CacheWeights "<<m_bCacheWeights<<"\n";
	if (m_bCacheWeights)
	{
		for (int i=0; i<2; i++)
		{
			SampledSurface * surface = m_pMultiConformalMap->GetSurface(i);
			bool exists = (m_SuggestedWeights.find(surface)!=m_SuggestedWeights.end());
			textStream<<"SurfaceExists "<<exists<<"\n";
			if (exists)
			{
				std::vector<SurfacePerVertexValue * > & precomputedWeights = m_SuggestedWeights[surface];
				if(N!=(int)precomputedWeights.size() && N!=0)
					std::cout<<"[ERROR] Precomputed = "<<precomputedWeights.size()<<" =? "<<N<<std::endl;
				assert(N==(int)precomputedWeights.size() && N!=0);	
				textStream<<"NumConfMaps "<<N<<"\n";
				for (int confMapID = 0; confMapID < N; confMapID++)
				{
					SurfacePerVertexValue * perVertexValue = precomputedWeights[confMapID];
					assert(perVertexValue!=NULL);
					perVertexValue->WritePerVertexValues(textStream);
				}
			}
		}
	}

	textStream<<"OptimalWeights "<<(m_mVertexToOptimalWeights.size()!=0)<<"\n";
	if (m_mVertexToOptimalWeights.size()!=0)
	{
		for (int surfID=0; surfID<2; surfID++)	
		{
			SampledSurface * surfPtr = m_pMultiConformalMap->GetSurface(surfID);
			assert(m_mVertexToOptimalWeights.find(surfPtr)!=m_mVertexToOptimalWeights.end());
			textStream<<"Surface "<<surfID;
			textStream<<" NumSamp "<<m_mVertexToOptimalWeights[surfPtr].size()<<"\n";
			assert(m_mVertexToOptimalWeights[surfPtr].size()!=0);
			std::map<int, std::vector<double> >::iterator iter;
			for (iter = m_mVertexToOptimalWeights[surfPtr].begin();
				 iter !=m_mVertexToOptimalWeights[surfPtr].end(); iter++)
			{
				assert((int)iter->second.size()==N);
				textStream<<"SampleVertex "<<iter->first<<" ";
				for (int i=0; i<N; i++)
					textStream<<iter->second[i]<<" ";
				textStream<<"\n";
			}
				
		}
	}
}

void MapConfidenceMultiConf::LoadConfidence(std::ifstream & textStream)
{
	int N = -1;
	std::string tempStr;
	// save values
	LoadConfidenceValues(textStream);
	
	// save names
	int numConfidences;
	textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "Confidences")==0);
	textStream>>numConfidences;		assert(numConfidences==m_Confidences.count());

	for (int i=0; i<m_Confidences.count(); i++)
	{
		textStream>>tempStr;
		assert(strcmp(tempStr.c_str(), m_Confidences[i].c_str())==0);
	}
	textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "Similarity")==0);
	textStream>>tempStr;	assert(strcmp(tempStr.c_str(), m_Similarity.c_str())==0);
	
	// save weights
	textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "CacheWeights")==0);
	bool cacheWeights;		textStream>>cacheWeights;	assert(cacheWeights==m_bCacheWeights);
	
	m_bWeightsReady = false;
	if (m_bCacheWeights)
	{
		m_bWeightsReady = true;
		for (int i=0; i<2; i++)
		{
			SampledSurface * surface = m_pMultiConformalMap->GetSurface(i);
			textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "SurfaceExists")==0);
			bool exists;			textStream>>exists;
			m_bWeightsReady = m_bWeightsReady && exists;
			if (exists)
			{
				std::vector<SurfacePerVertexValue * > & precomputedWeights = m_SuggestedWeights[surface];
				assert(precomputedWeights.size()==0);
				textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "NumConfMaps")==0);
				textStream>>N;
				for (int confMapID = 0; confMapID < N; confMapID++)
				{
					precomputedWeights.push_back(new SurfacePerVertexValue(surface, true));
					bool ok = precomputedWeights[confMapID]->LoadPerVertexValues(textStream);
					assert(ok);
				}		
			}
		}
	}
	
	textStream>>tempStr;		assert(strcmp(tempStr.c_str(), "OptimalWeights")==0);
	bool loadOptWeights;		textStream>>loadOptWeights; 
	if (loadOptWeights)
	{
		assert(N!=-1);
		for (int surfID=0; surfID<2; surfID++)	
		{
			SampledSurface * surfPtr = m_pMultiConformalMap->GetSurface(surfID);
			textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "Surface")==0);
			int loadSurfID;			textStream>>loadSurfID;		assert(loadSurfID==surfID);
			textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "NumSamp")==0);
			int numSamp;			textStream>>numSamp;

			for (int i=0; i<numSamp; i++)
			{
				textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "SampleVertex")==0);
				int sampVertex;		textStream>>sampVertex;
				m_mVertexToOptimalWeights[surfPtr][sampVertex] = std::vector<double>();
				
				std::vector<double> & currVect = m_mVertexToOptimalWeights[surfPtr][sampVertex];
				for (int mapID=0; mapID<N; mapID++)
				{
					double val;		textStream>>val;
					currVect.push_back(val);
				}
			}
		}
	}
	else
		m_bInterpolatedWeights = false;
}





