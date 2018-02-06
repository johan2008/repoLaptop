#include "PipelineExploreConfmaps.h"
#include "MapMultiConformal.h"
#include <algorithm>
#include "Sortable.h"
#include "LinAlgMatrixSparseReal.h"

//#define MAX_EIGENVECTORS 10
//// how many eigenvectors to take (prune based on eigenvalue and total number)
//#define FRAC_EIGENVECTOR_ENERGY .25		//.75
//
//// how many conformal maps to take in an eigenvector (relative to the top weights)
//#define FRAC_CONF_MAP_ENERGY .25			//.25
//#define SEED_FROM_CONF_ENERGY .75		//.75


bool OneStepSpectralConsistentConfMaps::AddMap(MapConformal* confMap, 
											   double eigenvectorValue, int id)
{
	std::vector<int>  * confMapCorrs;
	SurfaceSampleSet * confMapGenSet1;
	SurfaceSampleSet * confMapGenSet2;	
	confMap->GetGeneratingCorrespondences(&confMapCorrs, &confMapGenSet1, &confMapGenSet2);
	
	bool consistent = true;;
	for (int c=0; c<(int)confMapCorrs->size(); c+=2)
	{
		std::map<int, int>::iterator iter1 = m_Corrs12.find((*confMapCorrs)[c+0]);
		if (iter1!=m_Corrs12.end() 
			&& iter1->second!=(*confMapCorrs)[c+1])	// conflicts - do not add corr
		{
			consistent = false;
			break;
		}
		
		std::map<int, int>::iterator iter2 = m_Corrs21.find((*confMapCorrs)[c+1]);
		if (iter2!=m_Corrs21.end() 
			&& iter2->second!=(*confMapCorrs)[c+0])	// conflicts - do not add corr
		{
			consistent = false;
			break;
		}
	}
	
	if (consistent)
	{
		for (int c=0; c<(int)confMapCorrs->size(); c+=2)
		{
			m_Corrs12[(*confMapCorrs)[c+0]] = (*confMapCorrs)[c+1];
			m_Corrs21[(*confMapCorrs)[c+1]] = (*confMapCorrs)[c+0];
		}
	}
	else 
		return false;
	m_ConfMaps.push_back(confMap);
	m_OriginalIDs.push_back(id);
	m_EigenvectorValue.push_back(eigenvectorValue);
	//m_EigenvectorValue.push_back(1.);
	return true;
}

MapMultiConformal * OneStepSpectralConsistentConfMaps::CreateMultiConformalMap(std::map<VKString, bool> & existingCorrs,
																			   LinAlgMatrixReal & cachedConfVals,
																			   SurfaceDistance * dist1, double distNorm,
																			   SurfaceSampleSet * integrationSet,
																			   const VKString & multiConfWeights,
																			   const VKString & interpStr)
{
	int NEvenSamp = integrationSet->NumSamples();
	VKString encodedCorrs("");
	for (std::map<int, int>::iterator iter = m_Corrs12.begin(); 
		 iter!=m_Corrs12.end(); iter++)
		encodedCorrs=encodedCorrs+"_X_"+VKString::number(iter->first)+"->"+VKString::number(iter->second);

	if (existingCorrs.find(encodedCorrs)!=existingCorrs.end())
		return NULL;
	
	existingCorrs[encodedCorrs] = true;
	
	std::vector< std::map<int, double> > perPointWeights;
	for (int mapID=0; mapID<(int)m_EigenvectorValue.size(); mapID++)
	{
		perPointWeights.push_back(std::map<int, double>());
		double totalConf = 0;
		for (int p=0; p<NEvenSamp; p++)
		{
			int vID = integrationSet->GetSample(p).NearestVertex();
			double val = cachedConfVals(m_OriginalIDs[mapID], p);

			perPointWeights[mapID][vID] = val;
			totalConf += val;
		}
		std::cout<<"TotalConf Map = ["<<mapID<<"] = "<<totalConf<<std::endl;
	}
	
	std::cout<<"Creating Multi-Conformal Map:"<<std::endl;
	std::cout<<"\tWeights:"<<perPointWeights.size()<<"     Maps="<<m_ConfMaps.size()<<std::flush;
	std::cout<<"     Points="<<NEvenSamp<<std::endl;
	
	return new MapMultiConformal(m_ConfMaps, perPointWeights, multiConfWeights, interpStr);
}

void PipelineExploreConfmaps::OneStepSpectralLoadConformalMaps(const VKString & algorithmName)
{
	// fill m_pCollectionOfConformalMaps and m_pCachedConfidenceVals
	
	VKString genSetName = m_Params.GetStrValue(algorithmName, "SamplesMobVoteCast", valid);
	assert(valid);
	SurfaceSampleSet * genSet1 = GetSurface(0)->GetSampleSet(genSetName);
	assert(genSet1!=NULL);
	SurfaceSampleSet * genSet2 = GetSurface(1)->GetSampleSet(genSetName);	
	assert(genSet2!=NULL);	
//	int NGen1 = genSet1->NumSamples();
//	int NGen2 = genSet2->NumSamples();	
	
	VKString confidenceName = m_Params.GetStrValue(algorithmName, "MapEvaluator", valid);
	assert(valid);
	VKString evenSetName = m_Params.GetStrValue(algorithmName, "SamplesFineCorr", valid);
	assert(valid);
	
	SurfaceSampleSet * evenSet1 = GetSurface(0)->GetSampleSet(evenSetName);
	assert(evenSet1!=NULL);	
	SurfaceSampleSet * evenSet2 = GetSurface(1)->GetSampleSet(evenSetName);
	assert(evenSet2!=NULL);	
	int NEvenSamp = evenSet1->NumSamples();

	SurfaceMidEdgeConf * MConf1 = (SurfaceMidEdgeConf*)GetSurface(0);
	SurfaceMidEdgeConf * MConf2 = (SurfaceMidEdgeConf*)GetSurface(1);	

	VKString loadTriplets = m_Params.GetStrValue(algorithmName, "LoadConfMaps", valid);
	assert(valid);
	
	std::vector<int> triplets;
	std::vector<double> weight;
	std::ifstream tripletsStream(loadTriplets.c_str());
	assert(tripletsStream.is_open());
	//std::cout<<"Loading triplets:"<<std::endl;
	while (!tripletsStream.eof())
	{
		int a1, a2, a3;
		int b1, b2, b3;
		double w;
		tripletsStream>>a1>>b1>>a2>>b2>>a3>>b3>>w;
		if (a1==a2 || a1==a3 || a2 == a3)
			continue;
		if (w > 0)
		{
			//std::cout<<"\t["<<x<<" "<<y<<" "<<z<<" "<<w<<"]"<<std::endl;
			triplets.push_back(a1);
			triplets.push_back(b1);
			triplets.push_back(a2);
			triplets.push_back(b2);
			triplets.push_back(a3);
			triplets.push_back(b3);
			weight.push_back(w);
		}
	}
	
	VKString loadConfidences = m_Params.GetStrValue(algorithmName, "LoadConfidenceValues", valid);
	assert(!valid || loadConfidences=="none");	// cannot load for now - need for DB analysis project
	
	m_pCollectionOfConformalMaps = new MapScoredCollection(GetSurface(0), GetSurface(1));
	m_pCachedConfidenceVals = new LinAlgMatrixReal((int)weight.size(), NEvenSamp);
	TimeProfiler progressBar;
	std::cout<<"Loading Conformal Maps: "<<weight.size()<<std::endl;
	for (int i=0; i<(int)weight.size(); i++)	
	{
		double w = weight[i];
		std::vector<int> corrs;
		corrs.push_back(triplets[6*i + 0]);
		corrs.push_back(triplets[6*i + 1]);
		corrs.push_back(triplets[6*i + 2]);
		corrs.push_back(triplets[6*i + 3]);
		corrs.push_back(triplets[6*i + 4]);
		corrs.push_back(triplets[6*i + 5]);
		
		MapConformal * confMap = new MapConformal(MConf1, MConf2, corrs, genSetName);
		double ci=0;
		for (int p=0; p<NEvenSamp; p++)
		{
			SurfaceMapConfidence * confCalc = confMap->GetMapConfidenceCalculator(confidenceName);
			const SurfaceSample & s = evenSet1->GetSample(p);
			(*m_pCachedConfidenceVals)(i, p) = confCalc->CalculateConfidenceAtSample(s);
			ci += (*m_pCachedConfidenceVals)(i, p) * w;
		}
		m_pCollectionOfConformalMaps->AddMap(confMap, ci/NEvenSamp);
		
		if (i%10==0)
			progressBar.WriteProgress("CalcConf", i, (int)(weight.size()));
	}
}

void PipelineExploreConfmaps::OneStepSpectralGenerateConformalMaps(const VKString & algorithmName)
{
	//	srand(0);
	int KEEP1_AFTER_MCN = INT_MAX;
	int KEEP2_RANDOM	= 10000;	
	bool SUGGESTED_KEEP2_PER_CORR = -1;	// TODO
	int KEEP3_AFTER_CONFIDENCE = INT_MAX;
	
	VKString genSetName = m_Params.GetStrValue(algorithmName, "SamplesMobVoteCast", valid);
	assert(valid);
	SurfaceSampleSet * genSet1 = GetSurface(0)->GetSampleSet(genSetName);
	assert(genSet1!=NULL);
	SurfaceSampleSet * genSet2 = GetSurface(1)->GetSampleSet(genSetName);	
	assert(genSet2!=NULL);	
	int NGen1 = genSet1->NumSamples();
	int NGen2 = genSet2->NumSamples();	
	
	VKString confidenceName = m_Params.GetStrValue(algorithmName, "MapEvaluator", valid);
	assert(valid);
	VKString evenSetName = m_Params.GetStrValue(algorithmName, "SamplesFineCorr", valid);
	assert(valid);
	
	SurfaceSampleSet * evenSet1 = GetSurface(0)->GetSampleSet(evenSetName);
	assert(evenSet1!=NULL);	
	SurfaceSampleSet * evenSet2 = GetSurface(1)->GetSampleSet(evenSetName);
	assert(evenSet2!=NULL);	
	int NEvenSamp = evenSet1->NumSamples();
	SurfaceMidEdgeConf * MConf1 = (SurfaceMidEdgeConf*)GetSurface(0);
	SurfaceMidEdgeConf * MConf2 = (SurfaceMidEdgeConf*)GetSurface(1);	
	m_pCollectionOfConformalMaps = new MapScoredCollection(GetSurface(0), GetSurface(1));
	LinAlgComplex * mcn_orig1;
	LinAlgComplex * mcn_orig2;	
	LinAlgComplex * mcn_trans1;
	LinAlgComplex * mcn_trans2;	
	MapConfidenceMutuallyClosest::ConstructCoordsHolderForFastMCN(MConf1, MConf2, evenSet1, evenSet2,
																  &mcn_orig1, &mcn_orig2, 
																  &mcn_trans1, &mcn_trans2);
	std::vector<Sortable> conformalMaps;
	int MAX_CONF_MAPS = NGen1*(NGen1-1)*(NGen1-2) * NGen2*(NGen2-1)*(NGen2-2) / 6;
	int NUM_CORRS = NGen1 * NGen2;
	int KEEP2_PER_CORR = SUGGESTED_KEEP2_PER_CORR;
	if (KEEP2_PER_CORR<0)
	{
		KEEP2_PER_CORR = KEEP2_RANDOM / NUM_CORRS;
		KEEP2_RANDOM = KEEP2_PER_CORR * NUM_CORRS;
	}
	std::vector< std::vector<int> > corrToMapID;
	std::vector< std::vector<int> > mapIDToCorr;
	for (int i=0; i<MAX_CONF_MAPS; i++)
		mapIDToCorr.push_back(std::vector<int>());
	for (int i=0; i<NGen1*NGen2; i++)
		corrToMapID.push_back(std::vector<int>());
	
	TimeProfiler progressBar;
	std::cout<<"Generating Conformal Maps ("<<MAX_CONF_MAPS<<")"<<std::endl;
	int progress=0;
	for (int c1_1=0; c1_1<NGen1; c1_1++)
		for (int c2_1=c1_1+1; c2_1<NGen1; c2_1++)
			for (int c3_1=c2_1+1; c3_1<NGen1; c3_1++)
				for (int c1_2=0; c1_2<NGen2; c1_2++)
					for (int c2_2=0; c2_2<NGen2; c2_2++)
						for (int c3_2=0; c3_2<NGen2; c3_2++)
						{
							if (c1_2!=c2_2 && c1_2!=c3_2 && c2_2!=c3_2)
							{
								if (progress%500==0)
									progressBar.WriteProgress("CreateConfMaps", progress, MAX_CONF_MAPS);
								progress++;
								std::vector<int> corrs;
								corrs.push_back(c1_1);			corrs.push_back(c1_2);
								corrs.push_back(c2_1);			corrs.push_back(c2_2);
								corrs.push_back(c3_1);			corrs.push_back(c3_2);
								
								MapConformal * confMap = new MapConformal(MConf1, MConf2, corrs, genSetName);
								double score = progress;
								if (KEEP1_AFTER_MCN < MAX_CONF_MAPS)
								{
									score = MapConfidenceMutuallyClosest::GetFractionOfMutuallyClosest(MConf1, MConf2, 
																									   evenSet1, evenSet2, 
																									   mcn_orig1, mcn_orig2, 
																									   mcn_trans1, mcn_trans2,
																									   confMap);
								}
								int currID = (int)conformalMaps.size();
								int corr1ID = c1_1 * NGen2 + c1_2;
								int corr2ID = c2_1 * NGen2 + c2_2;
								int corr3ID = c3_1 * NGen2 + c3_2;
								assert(corr1ID < (int)corrToMapID.size());
								assert(corr2ID < (int)corrToMapID.size());
								assert(corr3ID < (int)corrToMapID.size());			
								corrToMapID[corr1ID].push_back(currID);
								corrToMapID[corr2ID].push_back(currID);
								corrToMapID[corr3ID].push_back(currID);	
								mapIDToCorr[currID].push_back(corr1ID);
								mapIDToCorr[currID].push_back(corr2ID);
								mapIDToCorr[currID].push_back(corr3ID);
								
								conformalMaps.push_back(Sortable(score, confMap, currID));
							}
						}
	std::cout<<std::endl;
	std::sort(conformalMaps.begin(), conformalMaps.end());
	std::cout<<"Before Prunning: "<<conformalMaps.size()<<std::endl;
	// First Pruning: based on fraction of MCN
	while ((int)conformalMaps.size() > KEEP1_AFTER_MCN)
	{
		// TODO check KEEP_PER_SAMPLE
		conformalMaps.erase(conformalMaps.begin());
	}
	
	// Second pruning: random
	while ((int)conformalMaps.size()>KEEP2_RANDOM)
	{
		double random = (double)rand() / (double)RAND_MAX;
		int delMapID = (int)(random * (double)conformalMaps.size()); 
		if (delMapID<0 || delMapID>=(int)conformalMaps.size())
			continue;
		
		conformalMaps.erase(conformalMaps.begin()+delMapID);
	}
	std::cout<<"After Prunning: "<<conformalMaps.size()<<std::endl;
	
	// Third Pruning: based on integrated confidences
	VKString useFuzzy = m_Params.GetStrValue("Pipeline", "UseFuzzy", valid);
	if (!valid)
		useFuzzy = "none";
	
	std::cout<<"Calculating Confidences"<<std::endl;
	LinAlgMatrixReal cachedConfidence1((int)conformalMaps.size(), NEvenSamp);
	for (int i=0; i<(int)conformalMaps.size(); i++)
	{
		MapConformal * confMap = (MapConformal*)conformalMaps[i].ptr;		
		double ci=0;
		for (int p=0; p<NEvenSamp; p++)
		{
			SurfaceMapConfidence * confCalc = confMap->GetMapConfidenceCalculator(confidenceName);
			const SurfaceSample & s = evenSet1->GetSample(p);
			if (useFuzzy!="FuzzyScaleConfidence")
				cachedConfidence1(i, p) = confCalc->CalculateConfidenceAtSample(s);
			else
				cachedConfidence1(i, p) = FuzzyConsistentConfidenceAtSample(confMap, s, confCalc);
			ci += cachedConfidence1(i, p);
			
		}
		conformalMaps[i].id = i;		
		conformalMaps[i].value = ci/NEvenSamp;
		
		if (i%10==0)
			progressBar.WriteProgress("CalcConf", i, (int)conformalMaps.size());
	}
	std::cout<<std::endl;
	
	std::sort(conformalMaps.begin(), conformalMaps.end());
	std::cout<<"Best Confidence = "<<conformalMaps[conformalMaps.size()-1].value<<std::endl;		
	while((int)conformalMaps.size() > KEEP3_AFTER_CONFIDENCE)
	{
		// TODO: check KEEP_PER_SAMPLE
		conformalMaps.erase(conformalMaps.begin());
	}
	
	// Construct the final map collection
	m_pCachedConfidenceVals = new LinAlgMatrixReal((int)conformalMaps.size(), NEvenSamp);
	for (int i=0; i<(int)conformalMaps.size(); i++)
	{
		m_pCollectionOfConformalMaps->AddMap((SurfaceMap*)conformalMaps[i].ptr, 0);
		int id = conformalMaps[i].id;
		for (int p=0; p<NEvenSamp; p++)
		{
			assert(cachedConfidence1(id, p)>=0);
			(*m_pCachedConfidenceVals)(i, p) = cachedConfidence1(id, p);
		}
	}
	
}

void PipelineExploreConfmaps::OneStepSpectralFillIntegratedBlendedMatrix(const VKString & algorithmName)
{
	TimeProfiler progressBar;
	VKString simName = m_Params.GetStrValue(algorithmName, "MapSimilarity", valid);
	assert(valid);
	
	VKString evenSetName = m_Params.GetStrValue(algorithmName, "SamplesFineCorr", valid);
	assert(valid);
	SurfaceSampleSet * evenSet1 = GetSurface(0)->GetSampleSet(evenSetName);
	assert(evenSet1!=NULL);
	int NEvenSamp = evenSet1->NumSamples();
	
	assert(m_pCachedConfidenceVals!=NULL);

	assert(m_pCollectionOfConformalMaps!=NULL);
	int NConfMaps = m_pCollectionOfConformalMaps->GetNumMaps();	
	double sigma2=pow(.5, 2.);
	SurfaceDistance * dist1 = GetSurface(0)->GetDistanceMetric("default");
	double distNorm = sqrt(GetSurface(0)->Area());
	int minNumShared=m_Params.GetIntValue(algorithmName, "MinNumShared", valid);
	if (!valid)
		minNumShared = 2;

	int numNZelements = 0;
	std::vector<int> pairsToExamine;
	for (int i=0; i<NConfMaps; i++)	
	{
		MapConformal * mi = (MapConformal*)m_pCollectionOfConformalMaps->GetMapByID(i);
		for (int j=0; j<NConfMaps; j++)
		{
			if (i%20==0 && j==0)
				progressBar.WriteProgress("NonzeroEntriesSearch", i * NConfMaps + j, NConfMaps*NConfMaps);

			MapConformal * mj = (MapConformal*)m_pCollectionOfConformalMaps->GetMapByID(j);
			bool consistent = MapSimilarityDistance::NumSharedCorrespondences(mi, mj)>=minNumShared;
//			std::cout<<"i: "<<i<<" j:"<<j<<" consistent: "<<consistent<<std::endl;
			if (consistent)
			{
				numNZelements++;
				pairsToExamine.push_back(i * NConfMaps + j);
			}
		}
	}
	std::cout<<"\n\tFound non-zero elements:"<<numNZelements<<" / "<<(NConfMaps*NConfMaps)<<" = ";
	std::cout<<((double)numNZelements / ((double)NConfMaps*NConfMaps))<<std::endl;
	m_pIntegratedBlendingMatrix = new LinAlgMatrixSparseReal(NConfMaps, NConfMaps, 
															 LinAlgMatrixSparseReal::SUPPORT_ALL_DYNAMIC_CAST,
															 numNZelements);
	
	std::cout<<"Filling non-zero entries:"<<std::endl;
	std::map<int, std::map<int, double> > nonZeroValues;
	int numAdded=0;
	for (int pairID=0; pairID<(int)pairsToExamine.size(); pairID++)
	{
		if (pairID%1000==0)
			progressBar.WriteProgress("FillMatrix", pairID, numNZelements);
		int i = pairsToExamine[pairID] / NConfMaps;
		int j = pairsToExamine[pairID] % NConfMaps;	
		
		if (i > j)	// must be precomputed
		{
			std::map<int, std::map<int, double> >::iterator iter1 = nonZeroValues.find(i);
			assert(iter1!=nonZeroValues.end());
			std::map<int, double>::iterator iter2 = iter1->second.find(j);
			assert(iter2!=iter1->second.end());
			m_pIntegratedBlendingMatrix->AddVal(j, i, iter2->second);
		}
		else	// need to compute
		{
			assert(i >=0 && i < m_pCollectionOfConformalMaps->GetNumMaps());
			assert(j >=0 && j < m_pCollectionOfConformalMaps->GetNumMaps());		
			
			MapConformal * mi = (MapConformal*)m_pCollectionOfConformalMaps->GetMapByID(i);
			MapConformal * mj = (MapConformal*)m_pCollectionOfConformalMaps->GetMapByID(j);
			
			double integral = 0;
			for (int p=0; p<NEvenSamp; p++)
			{
				double Sij = MapSimilarityDistance::FastSimilarity(sigma2, mi, mj,
																   evenSet1->GetSample(p), 
																   dist1, distNorm);
				double ci = (*m_pCachedConfidenceVals)(i, p);
				double cj = (*m_pCachedConfidenceVals)(j, p);
				
				assert(Sij*ci*cj>=0);
				integral += Sij*ci*cj;
			}
			//std::cout<<"B["<<j<<" "<<i<<"] = "<<integral<<std::endl;
			m_pIntegratedBlendingMatrix->AddVal(j, i, integral);
			if (i!=j)	// minor detail: don't cache diagonal entries
				nonZeroValues[j][i] = integral;
			if (i==j)
				numAdded++;
			else
				numAdded+=2;
		}
	}
	assert(numAdded==numNZelements);
	std::cout<<std::endl;	
	
	OneStepSpectralFindEigenvalues(algorithmName);		
}

void PipelineExploreConfmaps::OneStepSpectralFindEigenvalues(const VKString & algorithmName)
{
	int MAX_EIGENVECTORS = m_Params.GetIntValue(algorithmName, "MAX_EIGENVECTORS", valid);
	assert(valid);
	
	int NConfMaps = m_pCollectionOfConformalMaps->GetNumMaps();	
	MAX_EIGENVECTORS = vkMin(MAX_EIGENVECTORS, NConfMaps);
	std::cout<<"NEigen = "<<MAX_EIGENVECTORS<<std::endl;;
	// perform spectral analysis of blending matrix
	std::cout<<"Eigendecomposition of "<<NConfMaps<<" x "<<NConfMaps<<" matrix"<<std::endl;
	m_pEigenvaluesIBM = new LinAlgVectorReal(MAX_EIGENVECTORS);
	m_pEigenvectorsIBM = new LinAlgMatrixReal(NConfMaps, MAX_EIGENVECTORS);
	
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("EigenDecomposition");
	
	VKString svdType = m_Params.GetStrValue(algorithmName, "SVDSparse", valid);
	std::cout<<"SvdType="<<svdType.c_str()<<std::endl;
	if (svdType=="Matlab" || svdType=="matlab") 
	{
		VKString surfName = m_Params.GetStrValue(algorithmName, "Surface", valid);
		assert(valid);
		VKStringList meshes = m_Params.GetStrValues(surfName, "MeshName", valid);
		assert(valid);
		VKString tmpStr = ((meshes.count()==2) ? (meshes[0]+"_"+meshes[1]) : meshes[0]);
		m_pIntegratedBlendingMatrix->EigendecompositionSymmetricMatlab(*m_pEigenvaluesIBM, 
																	   *m_pEigenvectorsIBM, 
																	   MAX_EIGENVECTORS, 
																	   (tmpStr).c_str());
	}
	else if (svdType=="SVDLib")
	{
		std::cout<<"evals: "<<m_pEigenvaluesIBM->Rows()<<" evecs: "<<m_pEigenvectorsIBM->Rows();
		std::cout<<" x "<<m_pEigenvectorsIBM->Cols()<<" NEigs="<<MAX_EIGENVECTORS<<std::endl;
		m_pIntegratedBlendingMatrix->EigendecompositionSymmetricSVDLib(*m_pEigenvaluesIBM, 
																	   *m_pEigenvectorsIBM, 
																	   MAX_EIGENVECTORS);
	}
	else
	{
		std::cout<<"evals: "<<m_pEigenvaluesIBM->Rows()<<" evecs: "<<m_pEigenvectorsIBM->Rows();
		std::cout<<" x "<<m_pEigenvectorsIBM->Cols()<<" NEigs="<<MAX_EIGENVECTORS<<std::endl;

		m_pIntegratedBlendingMatrix->EigendecompositionSymmetric(*m_pEigenvaluesIBM, 
																 *m_pEigenvectorsIBM, 
																 MAX_EIGENVECTORS);
	}
	
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("EigenDecomposition", true);
}

void PipelineExploreConfmaps::OneStepSpectralCreateMultiConformalFromIBM(const VKString & algorithmName)
{
	int MAX_EIGENVECTORS = m_Params.GetIntValue(algorithmName, "MAX_EIGENVECTORS", valid);
	assert(valid);
	int NConfMaps = m_pCollectionOfConformalMaps->GetNumMaps();	
	MAX_EIGENVECTORS = vkMin(MAX_EIGENVECTORS, NConfMaps);
	std::cout<<"NEigen = "<<MAX_EIGENVECTORS<<std::endl;

	double FRAC_EIGENVECTOR_ENERGY = m_Params.GetDoubleValue(algorithmName, "FRAC_EIGENVECTOR_ENERGY", valid);
	assert(valid);
	double FRAC_CONF_MAP_ENERGY = m_Params.GetDoubleValue(algorithmName, "FRAC_CONF_MAP_ENERGY", valid);
	assert(valid);
	double SEED_FROM_CONF_ENERGY = m_Params.GetDoubleValue(algorithmName, "SEED_FROM_CONF_ENERGY", valid);
	assert(valid);
	
	VKString genSetName = m_Params.GetStrValue(algorithmName, "SamplesMobVoteCast", valid);
	assert(valid);
	SurfaceSampleSet * genSet1 = GetSurface(0)->GetSampleSet(genSetName);
	assert(genSet1!=NULL);
	SurfaceSampleSet * genSet2 = GetSurface(1)->GetSampleSet(genSetName);
	assert(genSet2!=NULL);
//	int NGen1 = genSet1->NumSamples();
	
	VKString evenSetName = m_Params.GetStrValue(algorithmName, "SamplesFineCorr", valid);
	assert(valid);
	SurfaceSampleSet * evenSet1 = GetSurface(0)->GetSampleSet(evenSetName);
	assert(evenSet1!=NULL);	
//	int NEvenSamp = evenSet1->NumSamples();
	
	assert(m_pEigenvaluesIBM!=NULL);
	assert(m_pEigenvectorsIBM!=NULL);	
		
	int eigV = 0;
	
	VKString finalChoiceFunction = m_Params.GetStrValue(algorithmName, "RefinementMapEval", valid);
	assert(valid);
	VKString interpStr =  m_Params.GetStrValue(algorithmName, "MultiConfInterp", valid);
	assert(valid);		
	VKString multiConfWeights =  m_Params.GetStrValue(algorithmName, "MultiConfWeights", valid);
	assert(valid);		
	
	m_pCollectionOfMaps = new MapScoredCollection(GetSurface(0), GetSurface(1));
	
	std::vector<int> * confMapCorrs;
	SurfaceSampleSet * confMapGenSet1;
	SurfaceSampleSet * confMapGenSet2;	
	std::map<VKString, bool> encodedCorrsExist;	
	//std::cout<<"Eigenvalue[0]="<<(*m_pEigenvaluesIBM)(0)*FRAC_EIGENVECTOR_ENERGY<<std::endl;
	bool ADD_CONFORMAL_MAPS = true;
	
	std::vector<OneStepSpectralConsistentConfMaps*> consistentMaps;
	// temp stuff:
	std::vector<int> orderOfConfMaps;
	std::vector<bool> addedConfMap;
	for (int i=0; i<NConfMaps; i++)
		addedConfMap.push_back(false);

	VKString useFuzzy = m_Params.GetStrValue("Pipeline", "UseFuzzy", valid);
	if (!valid)
		useFuzzy = "none";
	
	
	// analyze top eigenvectors
	while (eigV < MAX_EIGENVECTORS && eigV < (*m_pEigenvaluesIBM).Rows()
		   && (*m_pEigenvaluesIBM)(eigV) >= (*m_pEigenvaluesIBM)(0)*FRAC_EIGENVECTOR_ENERGY)
	{		
		std::cout<<"Creating eigenvector: "<<eigV<<" val="<<(*m_pEigenvaluesIBM)(eigV)<<std::flush;
		std::cout<<" MaxVal="<<((*m_pEigenvaluesIBM)(0)*FRAC_EIGENVECTOR_ENERGY)<<std::endl;
		std::vector<Sortable> sortedEigenvector;
		for (int i=0; i<NConfMaps; i++)
		{
			double val = (*m_pEigenvectorsIBM)(i, eigV);
			if (val < 0)
				val *= -1;
			sortedEigenvector.push_back(Sortable(val, NULL, i));
		}
		std::sort(sortedEigenvector.begin(), sortedEigenvector.end());
		double maxValue = sortedEigenvector[(int)sortedEigenvector.size()-1].value*FRAC_CONF_MAP_ENERGY;
		double maxSeedVal = sortedEigenvector[(int)sortedEigenvector.size()-1].value*SEED_FROM_CONF_ENERGY;
		
		// analyze eigenvector
		for (int i=0; i<(int)sortedEigenvector.size(); i++)
		{
			int id = sortedEigenvector[sortedEigenvector.size()-i-1].id;
			double val = sortedEigenvector[sortedEigenvector.size()-i-1].value;

			if (val < maxValue)
				break;
			//val = 1.;
			MapConformal * confMap = (MapConformal*) m_pCollectionOfConformalMaps->GetMapByID(id);
			confMap->GetGeneratingCorrespondences(&confMapCorrs, &confMapGenSet1, &confMapGenSet2);
			assert(confMapGenSet1==genSet1);
			assert(confMapGenSet2==genSet2);
			//bool conflictingCorrs = false;
			if (!ADD_CONFORMAL_MAPS)	// ADD_CORRESPONDENCES
			{
				assert(false);
			}
			else 
			{				
				bool couldAdd = false;
				int tryingToAddID=0;
				while(!couldAdd && tryingToAddID < (int)consistentMaps.size())
				{
					couldAdd = consistentMaps[tryingToAddID++]->AddMap(confMap, val, id);
				}
				if (!couldAdd && val > maxSeedVal)
				{
					consistentMaps.push_back(new OneStepSpectralConsistentConfMaps());
					consistentMaps[(int)consistentMaps.size()-1]->AddMap(confMap, val, id);
				}
			}
		}
			
		if (!ADD_CONFORMAL_MAPS)	// ADD_CORRESPONDENCES
		{
			assert(false);
		}
		
		std::cout<<"\tEncountered "<<consistentMaps.size()<<" maps in eigenvector "<<eigV<<std::endl;
		SurfaceDistance * dist1 = GetSurface(0)->GetDistanceMetric("default");
		double distNorm = sqrt(GetSurface(0)->Area());
		for (int i=0; i<(int)consistentMaps.size(); i++)
		{
			// order: temp stuff
			for (int j=0; j<(int)consistentMaps[i]->m_OriginalIDs.size(); j++)
			{
				int tempAddingID = consistentMaps[i]->m_OriginalIDs[j];
				if (!addedConfMap[tempAddingID])
				{
					orderOfConfMaps.push_back(tempAddingID);
					addedConfMap[tempAddingID] = true;
				}
			}
			
			MapMultiConformal * multiConf = consistentMaps[i]->CreateMultiConformalMap(encodedCorrsExist,
																					   *m_pCachedConfidenceVals,
																					   dist1, distNorm,
																					   evenSet1,
																					   multiConfWeights,
																					   interpStr);
			std::cout<<"\tMap "<<i<<" - Added "<<multiConf<<std::endl;
			std::cout<<"\tCorrs: "<<std::endl;
			for (std::map<int, int>::iterator iter = consistentMaps[i]->m_Corrs12.begin();
				 iter != consistentMaps[i]->m_Corrs12.end(); iter++)
				std::cout<<"\t\t"<<iter->first<<" -> "<<iter->second<<std::endl;
			if (multiConf!=NULL)
			{
				int newMapID = -1;
				double assignedScore;
				if (useFuzzy=="none")
				{
					newMapID = m_pCollectionOfMaps->AddMap(multiConf, finalChoiceFunction);
					m_pCollectionOfMaps->GetMapByID(newMapID, &assignedScore);
				}
				else if (useFuzzy=="FuzzyFinalChoice"
						 || useFuzzy=="FuzzyScaleConfidence")
				{
					assignedScore = FuzzyConsistentConfidence(multiConf, evenSetName,
															  finalChoiceFunction);
					m_pCollectionOfMaps->AddMap(multiConf, assignedScore);
				}
				else
					assert(false);
				
				m_MapIDToEigenvalue.push_back(Sortable(-assignedScore, multiConf, eigV));
			}
		}
		consistentMaps.clear();
//		std::cout<<"\tCollection Size = "<<m_pCollectionOfMaps->GetNumMaps()<<std::endl;
		
		eigV++;
	}

	std::sort(m_MapIDToEigenvalue.begin(), m_MapIDToEigenvalue.end());
	m_pCollectionOfMaps->SortMaps();
	m_pTheFinalMap = m_pCollectionOfMaps->GetBestMap();
	std::cout<<"Created "<<m_pCollectionOfMaps->GetNumMaps()<<" maps"<<std::endl;
	
	// temp stuff, removeme
//	std::ofstream mapOrderFile("MapOrderFile.txt");
//	for (int i=0; i<(int)orderOfConfMaps.size(); i++)
//		mapOrderFile<<orderOfConfMaps[i]<<" ";
//	mapOrderFile.close();
}

double PipelineExploreConfmaps::GetFuzzyValue(int vertex1, int vertex2)
{
	// just a sanity check in case fuzzy corrs were loaded, but do not exist
	assert(m_VertexToFuzzySetID1!=NULL && m_VertexToFuzzyVertexID2!=NULL);
	
	if (m_aFuzzyValues==NULL)	// were loaded / are not defined
		return .1;
	assert(m_aFuzzyValues!=NULL);
	int setID1 = m_VertexToFuzzySetID1[vertex1];
	int vtxID2 = m_VertexToFuzzyVertexID2[vertex2];
	return .1 + .9 * m_aFuzzyValues[setID1][vtxID2];
}

void PipelineExploreConfmaps::ReadFuzzyCorrespondences(const VKString & filename)
{
	FILE * fuzzyFile = fopen(filename.c_str(), "rb");
	if (fuzzyFile==NULL)
	{
		std::cout<<"[ERROR] Could not read fuzzy corrs: "<<filename.c_str()<<std::endl;
		assert(fuzzyFile!=NULL);
	}
	
	std::map<int, std::map<int, double> > fuzzyMap;
	std::map<int, bool> toSetAdded;	// stores set of points on model 2
	std::vector<int> toSet;
	std::vector<int> fromSet;
		
	// read fuzzy values
	int nModels;
	int mID1, mID2;
	fread(&nModels, sizeof(int), 1, fuzzyFile);
	assert(nModels==1);
	int nPnts1;
	fread(&mID1, sizeof(int), 1, fuzzyFile);
	assert(mID1==0);
	fread(&nPnts1, sizeof(int), 1, fuzzyFile);
	
	int vID1, vID2, nCorrs;
	double val;
	int numberofNullPnts=0;
	for (int i=0; i<nPnts1; i++)
	{
		fread(&vID1, sizeof(int), 1, fuzzyFile);
		fread(&nCorrs, sizeof(int), 1, fuzzyFile);
		fromSet.push_back(vID1);
		if (nCorrs==0)
			numberofNullPnts++;
		for (int cID=0; cID < nCorrs; cID++)
		{
			fread(&mID2, sizeof(int), 1, fuzzyFile);
			assert(mID2==1);
			fread(&vID2, sizeof(int), 1, fuzzyFile);
			fread(&val, sizeof(double), 1, fuzzyFile);
			fuzzyMap[vID1][vID2] = val;
			if (toSetAdded.find(vID2)==toSetAdded.end())
			{
				toSetAdded[vID2] = true;
				toSet.push_back(vID2);
			}
		}
	}
	fclose(fuzzyFile);
	if (numberofNullPnts>0)
		std::cout<<"[WARNING] number of non-defined pnts="<<numberofNullPnts<<std::endl;
	if (toSet.size()==0)
	{
		std::cout<<"[WARNING] Fuzzy correspondences are not defined"<<std::endl;
		m_VertexToFuzzySetID1 = new int[1];
		m_VertexToFuzzyVertexID2 = new int[1];
		m_aFuzzyValues = NULL;
		return;
	}
		
	DistanceOnTheFly * distance1 = GetSurface(0)->GetOnTheFlyDistanceMetric(-1, "default", -1);
	DistanceOnTheFly * distance2 = GetSurface(1)->GetOnTheFlyDistanceMetric(-1, "default", -1);

	int nVert1 = GetSurface(0)->GetMesh()->NVertices();
	int nVert2 = GetSurface(1)->GetMesh()->NVertices();
	
	m_VertexToFuzzySetID1 = new int[nVert1];
	m_VertexToFuzzyVertexID2 = new int[nVert2];
	
	for (int surfID=0; surfID < 2; surfID++)
	{
		std::vector<int> currSet = surfID==0 ? fromSet : toSet;
		int * currArrayMap = surfID==0 ? m_VertexToFuzzySetID1 : m_VertexToFuzzyVertexID2;
		DistanceOnTheFly * currDist = surfID==0 ? distance1 : distance2;
		int currNVert = surfID==0 ? nVert1 : nVert2;
		for (int i=0; i<(int)currSet.size(); i++)
			currDist->PrecomputeRow(currSet[i]);
		
		for (int i=0; i<currNVert; i++)
		{
			SurfaceSample currVtxSamp(i, GetSurface(surfID)->GetMesh());
			double minDist = FLT_MAX;
			int bestPnt = -1;
			for (int j=0; j<(int)currSet.size(); j++)
			{
				SurfaceSample mySamp = SurfaceSample(currSet[j], GetSurface(surfID)->GetMesh());
				double dist = currDist->Distance(mySamp, currVtxSamp);
				if (dist < minDist)
				{
					minDist = dist;
					bestPnt = j;
				}
			}
			assert(bestPnt!=-1);
			if (surfID==0)
				currArrayMap[i] = bestPnt;
			else
				currArrayMap[i] = currSet[bestPnt];
		}
	}	
	
	// fill fuzzy values
	m_aFuzzyValues = new std::map<int, double>[fromSet.size()];
	for (int i=0; i<(int)fromSet.size(); i++)
		m_aFuzzyValues[i] = fuzzyMap[fromSet[i]];
}

double PipelineExploreConfmaps::FuzzyConsistentConfidenceAtSample(SurfaceMap * surfMap, 
													const SurfaceSample & samp, 
													SurfaceMapConfidence * confidence)
{
	int v1, v2;
	if (samp.CheckMeshIsSame(GetSurface(0)->GetMesh()))
	{
		v1 = samp.NearestVertex();
		SurfaceSample forwSamp = surfMap->ForwardMap(samp);
		assert(!forwSamp.Invalid());
		v2 = forwSamp.NearestVertex();
	}
	else 
	{
		v2 = samp.NearestVertex();
		SurfaceSample invSamp = surfMap->InverseMap(samp);
		assert(!invSamp.Invalid());
		v1 = invSamp.NearestVertex();		
	}
	
	double fuzzyValue = GetFuzzyValue(v1, v2);
	
	return fuzzyValue * confidence->ConfidenceAtSample(samp);
}

double PipelineExploreConfmaps::FuzzyConsistentConfidence(SurfaceMap * surfMap,
														  const VKString & sampSetName,
														  const VKString & confidenceName)
{
	SurfaceMapConfidence * confidence = surfMap->GetMapConfidenceCalculator(confidenceName);

	double total=0;	
	int nSurf = 1;
	for (int i=0; i<nSurf; i++)
	{		
		SurfaceSampleSet * sampleSet = surfMap->GetSurface(i)->GetSampleSet(sampSetName);
		assert(sampleSet!=NULL);
		double lt = 0;
		for (int j=0; j<sampleSet->NumSamples(); j++)
		{
			const SurfaceSample & samp = sampleSet->GetSample(j);
			lt += FuzzyConsistentConfidenceAtSample(surfMap, samp, confidence);
		}
		
		total += lt / sampleSet->NumSamples();
	}		
	assert(nSurf>=1);
	return total / nSurf;
}


void PipelineExploreConfmaps::LoadAndProcessOneStepSpectralMap(const VKString & algorithmName)
{
	int MAX_EIGENVECTORS = m_Params.GetIntValue(algorithmName, "MAX_EIGENVECTORS", valid);
	assert(valid);
	
	WriteLog();	
	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("SurfaceCorrespondenceConfMap");	
	VKStringList loads = m_Params.GetStrValues("Pipeline", "Load", valid);
	VKStringList process = m_Params.GetStrValues("Pipeline", "Process", valid);		

	VKString useFuzzy = m_Params.GetStrValue("Pipeline", "UseFuzzy", valid);
	if (!valid)
		useFuzzy = "none";
	
	if (useFuzzy!="none")
	{
		VKString fuzzyFile = m_Params.GetStrValue("Pipeline", "FuzzyMapFile", valid);
		assert(valid);
		ReadFuzzyCorrespondences(fuzzyFile);
	}
	
	if (loads.contains("Maps") || process.contains("Maps"))
	{
		std::cout<<"Creating maps using one step spectral analysis"<<std::endl;
		VKString evenSetName = m_Params.GetStrValue(algorithmName, "SamplesFineCorr", valid);
		assert(valid);
		SurfaceSampleSet * evenSet1 = GetSurface(0)->GetSampleSet(evenSetName);
		assert(evenSet1!=NULL);
		int NEvenSamp = evenSet1->NumSamples();
		
		VKStringList mapProcess = m_Params.GetStrValues("Pipeline", "MapProcess", valid);
		
		// load confidences and similarities
		VKStringList confidencesNames = m_Params.GetStrValues(algorithmName, "LoadConfidences", valid);	assert(valid);
		VKStringList similaritiesNames = m_Params.GetStrValues(algorithmName, "LoadSimilarities", valid); assert(valid);		
		for (int i=0; i<confidencesNames.count(); i++)
			AddClonableConfidence(confidencesNames[i]);
		for (int i=0; i<similaritiesNames.count(); i++)
			AddClonableSimilarity(similaritiesNames[i]);
		
		// try to load conformal map and confidence values
		if (!process.contains("Maps") || !mapProcess.contains("ExploreConformals"))
		{	
			std::cout<<"Loading conformal maps"<<std::endl;					
			SurfaceMap * tempMap = SurfaceMap::CreateMap(LocalFilename(".conformal.collection"),
														 GetSurface(0), GetSurface(1));			
			if (tempMap!=NULL)
			{
				assert(tempMap->GetSurfaceMapType()=="MapScoredCollection");
				m_pCollectionOfConformalMaps = (MapScoredCollection*)tempMap;
				
				int NConfMaps = m_pCollectionOfConformalMaps->GetNumMaps();	
				m_pCachedConfidenceVals = new LinAlgMatrixReal(NConfMaps, NEvenSamp);
				if (!m_pCachedConfidenceVals->ReadMatrixBin(LocalFilename(".confidence.matrix").c_str()))
				{
					delete m_pCachedConfidenceVals;
					m_pCachedConfidenceVals = NULL;
				}
			}
		}
		
		// generate conformal randomized maps (remember c_i)
		if (m_pCollectionOfConformalMaps==NULL || m_pCachedConfidenceVals==NULL)
		{
			AnalysisStats::m_GlobalStats.m_Timing.startedProcess("AABasicCreateMapsSearch");
			VKString loadConfMaps = m_Params.GetStrValue(algorithmName, "LoadConfMaps", valid);
			std::cout<<"Processing Conformal Maps: Load="<<loadConfMaps.c_str()<<std::endl;
			if (valid && loadConfMaps!="none")
				OneStepSpectralLoadConformalMaps(algorithmName);
			else						
				OneStepSpectralGenerateConformalMaps(algorithmName);
			assert(m_pCachedConfidenceVals!=NULL);
			assert(m_pCollectionOfConformalMaps->GetNumMaps()>0);
			m_pCachedConfidenceVals->WriteMatrixBin(LocalFilename(".confidence.matrix").c_str());
			((SurfaceMap*)m_pCollectionOfConformalMaps)->SaveMap(LocalFilename(".conformal.collection"));
			AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("AABasicCreateMapsSearch");			
		}
		
		// fill integrated blending matrix
		if (!process.contains("Maps") || !mapProcess.contains("Aggregate"))
		{
			int NConfMaps = m_pCollectionOfConformalMaps->GetNumMaps();			
			std::cout<<"Loading Blending Matrix"<<std::endl;
			m_pIntegratedBlendingMatrix = new LinAlgMatrixSparseReal(NConfMaps, NConfMaps, 
																	 LinAlgMatrixSparseReal::SUPPORT_ALL_DYNAMIC_CAST);
			if (!m_pIntegratedBlendingMatrix->ReadMatrixBin(LocalFilename(".blending.matrix").c_str()))
			{
				delete m_pIntegratedBlendingMatrix;
				m_pIntegratedBlendingMatrix = NULL;
			}
			else
			{
				std::cout<<"Loading eigenstructure"<<std::endl;
				m_pEigenvectorsIBM = new LinAlgMatrixReal(NConfMaps, MAX_EIGENVECTORS);
				if (m_pEigenvectorsIBM->ReadMatrixBin(LocalFilename(".eigenvectors.matrix").c_str()))
				{
					m_pEigenvaluesIBM = new LinAlgVectorReal(MAX_EIGENVECTORS);
					if (!m_pEigenvaluesIBM->ReadVectorBin(LocalFilename(".eigenvalues.vector").c_str()))
					{
						delete m_pEigenvectorsIBM;
						delete m_pEigenvaluesIBM;
						m_pEigenvectorsIBM = NULL;
						m_pEigenvaluesIBM = NULL;
					}
				}
				else
				{
					delete m_pEigenvectorsIBM;
					m_pEigenvectorsIBM = NULL;
				}				
			}
		}
		
		if (m_pEigenvectorsIBM==NULL || m_pEigenvaluesIBM==NULL)
		{			
			AnalysisStats::m_GlobalStats.m_Timing.startedProcess("AABasicCreateMapsCluster");
			std::cout<<"Processing Blending Matrix"<<std::endl;
			if (m_pIntegratedBlendingMatrix!=NULL)
				OneStepSpectralFindEigenvalues(algorithmName);
			else
				OneStepSpectralFillIntegratedBlendedMatrix(algorithmName);
			assert(m_pIntegratedBlendingMatrix!=NULL);
			assert(m_pEigenvectorsIBM!=NULL);
			assert(m_pEigenvaluesIBM!=NULL);
			m_pIntegratedBlendingMatrix->WriteMatrixBin(LocalFilename(".blending.matrix").c_str());	
			m_pEigenvectorsIBM->WriteMatrixBin(LocalFilename(".eigenvectors.matrix").c_str());
			m_pEigenvaluesIBM->WriteVectorBin(LocalFilename(".eigenvalues.vector").c_str());						
			//m_pIntegratedBlendingMatrix->WriteMatrixASCII(LocalFilename(".blending.ascii.matrix.rows").c_str(),
			//											  LocalFilename(".blending.ascii.matrix.cols").c_str(),
			//											  LocalFilename(".blending.ascii.matrix.vals").c_str());
			//m_pEigenvectorsIBM->WriteMatrixASCII(LocalFilename(".eigenvectors.ascii.matrix").c_str());
			//m_pEigenvaluesIBM->WriteVectorASCII(LocalFilename(".eigenvalues.ascii.vector").c_str());						
			AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("AABasicCreateMapsCluster");
		}
		
		if (!process.contains("Maps") || !mapProcess.contains("Extrapolate"))
		{
			std::cout<<"Loading Final Maps"<<std::endl;
			SurfaceMap * tempCollect = SurfaceMap::CreateMap(LocalFilename(".finalmap.collection"), 
															 GetSurface(0), GetSurface(1));
			m_pTheFinalMap = SurfaceMap::CreateMap(LocalFilename(".theFinal.map"), 
												   GetSurface(0), GetSurface(1));
			if (tempCollect!=NULL)
				m_pCollectionOfMaps = (MapScoredCollection*)tempCollect;
		}
		
		if (m_pTheFinalMap==NULL || m_pCollectionOfMaps==NULL)
		{
			AnalysisStats::m_GlobalStats.m_Timing.startedProcess("AABasicCreateMapsWeights");
			std::cout<<"Generating final maps"<<std::endl;			
			OneStepSpectralCreateMultiConformalFromIBM(algorithmName);
			assert(m_pTheFinalMap!=NULL);
			assert(m_pCollectionOfMaps!=NULL);		
			((SurfaceMap*)m_pCollectionOfMaps)->SaveMap(LocalFilename(".finalmap.collection"));
			((SurfaceMap*)m_pTheFinalMap)->SaveMap(LocalFilename(".theFinal.map"));
			AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("AABasicCreateMapsWeights");
		}
	}
	WriteLog();	
	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("SurfaceCorrespondenceConfMap");
}



    ///////	************************	 ///////
   ///////							    ///////
  /////// PER POINT OPTIMIZATION STUFF ///////
 ///////					  		  ///////
///////    ************************  ///////
void PipelineExploreConfmaps::MapSelectedBasedOnIntegratedBlendedMatrix(ParamParser * drawParams)
{
	if (m_pEigenvectorsIBM==NULL)
	{
		std::cout<<"[WARNING] Eigenvectors of integrated blended matrix not defined"<<std::endl;
		return;
	}
	// some initialization
	int v = GetSurface(0)->GetSelectedVertex(drawParams);
	if (v==-1)
	{
		std::cout<<"\tVertex Not Selected"<<std::endl;
		return;
	}
	SurfaceSample sampV(v, GetSurface(0)->GetMesh());
	
	std::vector<double> weights;
	std::vector<SurfaceSample> samples;
	std::vector<Sortable> weightedSamples;
	
	// check top eigenvector
	double sign = ((*m_pEigenvectorsIBM)(0, 0) < 0) ? -1 : 1;
	for (int i=0; i<m_pCollectionOfConformalMaps->GetNumMaps(); i++)
	{
		//std::cout<<"i="<<i<<" sign="<<sign<<" val="<<(*m_pEigenvectorsIBM)(i, 0)<<std::endl;
		//assert( (*m_pEigenvectorsIBM)(i, 0)*sign >= 0);
		if ((*m_pEigenvectorsIBM)(i, 0)*sign < 0)
		{
			std::cout<<"[WARNING] inconsistent signs of the top eigenvector: ";
			std::cout<<(*m_pEigenvectorsIBM)(i,0)<<". i="<<i<<std::endl;
		}
	}
	
	int eigVID = 0;	
	for (int i=0; i<m_pCollectionOfConformalMaps->GetNumMaps(); i++)
	{
		double val = (*m_pEigenvectorsIBM)(i, eigVID);
		if (val < 0)
			val *= -1;
		weights.push_back(val);
		samples.push_back(m_pCollectionOfConformalMaps->GetMapByID(i)->ForwardMap(sampV));
		weightedSamples.push_back(Sortable(val, NULL, i));
	}
	
	std::sort(weightedSamples.begin(), weightedSamples.end());
	
	MapConfidenceMultiConf::NormalizeWeights(weights);
	m_MapsDueToTopEigenvectors.clear();
	m_EigenValues.clear();
	int showTopMaps = -1;
	//	double weightThreshold = 1. / m_pCollectionOfConformalMaps->GetNumMaps();
	double weightThreshold = weights[weightedSamples[(int)weightedSamples.size()-1].id] * .5;
	
	if (showTopMaps==-1)
		showTopMaps = m_pCollectionOfConformalMaps->GetNumMaps();
	
	for (int i=0; i<showTopMaps; i++)
	{
		int id = weightedSamples[(int)weightedSamples.size()-i-1].id;
		if (weights[id] <= weightThreshold)
			break;
		m_MapsDueToTopEigenvectors.push_back(samples[id]);
		m_EigenValues.push_back(weights[id]);		
	}
}

void PipelineExploreConfmaps::PrintIntegratedBlendingMatrix(ParamParser * drawParams)
{
	assert(false);	
}

void PipelineExploreConfmaps::PrintBlendingMatrixAtCurrentPoint(ParamParser * drawParams)
{
	assert(false);
}



