#include <iostream>
#include "SampleSetMapGenerator.h"
#include "TimeProfiler.h"
/////////////////// Generator abstract class ////////////////

SampleSetMapGenerator::SampleSetMapGenerator(SurfaceSampleSet * set1, SurfaceSampleSet * set2, 
											 int numPntsPerMap, int maxNumMaps,
											 GenerationMethod method,
											 double surface1Area, double surface2Area)
{
	m_bEnforceSameSampleType = true;
	m_pSet1 = set1;
	m_pSet2 = set2;	
	m_GenSearchMethod = method;
	m_iCurrentMapSequence = 0;
	m_fDistanceNormalization1 = sqrt(surface1Area);
	m_fDistanceNormalization2 = sqrt(surface2Area);
	
	m_iNumPointsPerMap = numPntsPerMap;
	m_iMaxNumMaps = maxNumMaps;
}

SampleSetMapGenerator::SampleSetMapGenerator()
{
	m_bEnforceSameSampleType = true;	
	m_pSet1 = NULL;
	m_pSet2 = NULL;	
	m_iCurrentMapSequence = 0;
	m_fDistanceNormalization1 = 0;
	m_fDistanceNormalization2 = 0;
	
	m_iNumPointsPerMap = 0;
	m_iMaxNumMaps = 0;
}

std::vector<int> & SampleSetMapGenerator::GenerateNextMap(bool * valid)
{
	if (valid!=NULL)
	{
		if (m_iCurrentMapSequence>=(int)m_MapSequencesAfterThresholding.size())
			*valid = false;
		else
			*valid = true;
	}
	
	if (m_iCurrentMapSequence>=0 && m_iCurrentMapSequence<(int)m_MapSequencesAfterThresholding.size())
		return m_MapSequencesAfterThresholding[m_iCurrentMapSequence++];
	
	return m_Empty;
}

int SampleSetMapGenerator::StartGeneration()
{
	m_iCurrentMapSequence = 0;
	int attempts=0;
	if (m_GenSearchMethod==GEN_RANDOM_SEARCH)
	{
		TimeProfiler profiler;
		std::vector<int> potentialMap;
		while ((int)m_MapSequencesAfterThresholding.size()<m_iMaxNumMaps
			   && attempts++ < m_iMaxNumMaps*100)
		{
			if (m_MapSequencesAfterThresholding.size() % 100000 == 0)
				profiler.WriteProgress("GeneratingSet", (int)m_MapSequencesAfterThresholding.size(), 
									   m_iMaxNumMaps);
			int A = rand()%m_pSet1->NumSamples();
			int B = rand()%m_pSet1->NumSamples();
			int C = rand()%m_pSet1->NumSamples();
			int Ap = rand()%m_pSet2->NumSamples();
			int Bp = rand()%m_pSet2->NumSamples();
			int Cp = rand()%m_pSet2->NumSamples();			
			potentialMap.clear();
			if (m_iNumPointsPerMap==3)
			{
				//if (FillTripletAndTest(A, B, C, Ap, Bp, Cp, potentialMap))
				FillTripletAndTest(A, B, C, Ap, Bp, Cp, potentialMap);	// add all - only used in mob votign
				m_MapSequencesAfterThresholding.push_back(potentialMap);	
			}
			else 
				assert(false);
		}
		std::cout<<std::endl;
	}
	else if (m_GenSearchMethod==GEN_EXHAUSTIVE_SEARCH)
	{
		int exploredMaps=0;
		std::vector<int> potentialMap;		
		if (m_iNumPointsPerMap==3)
		{
			for (int A=0; A<m_pSet1->NumSamples(); A++)
				for (int B=A+1; B<m_pSet1->NumSamples(); B++)
					for (int C=B+1; C<m_pSet1->NumSamples(); C++)
						for (int Ap=0; Ap<m_pSet2->NumSamples(); Ap++)
							for (int Bp=0; Bp<m_pSet2->NumSamples(); Bp++)
								for (int Cp=0; Cp<m_pSet2->NumSamples(); Cp++)
								{
									potentialMap.clear();
									if (Ap!=Bp && Ap!=Cp && Bp!=Cp)
										exploredMaps++;
									if (FillTripletAndTest(A, B, C, Ap, Bp, Cp, potentialMap))
										m_MapSequencesAfterThresholding.push_back(potentialMap);
									if ((int)m_MapSequencesAfterThresholding.size()>=m_iMaxNumMaps)
										return (int)m_MapSequencesAfterThresholding.size();
								}
		}
		else
			assert(false);
		
		std::cout<<"Potential Maps : "<<exploredMaps<<" Accepted: "<<m_MapSequencesAfterThresholding.size()<<std::endl;
	}
	
	// special case: cannot create a generator with given thresholds.
	if (m_MapSequencesAfterThresholding.size()==0)
	{
		bool nonzeroThresh=false;
		for (int i=0; i<(int)m_MinDistanceThresholds.size(); i++)
		{
			if (!nonzeroThresh)
				nonzeroThresh = m_MinDistanceThresholds[i]>0;
			m_MinDistanceThresholds[i] = 0;
		}
		for (int i=0; i<(int)m_DistanceToOtherRatioThreshold.size(); i++)
		{
			if (!nonzeroThresh)
				nonzeroThresh = m_DistanceToOtherRatioThreshold[i]>0;			
			m_DistanceToOtherRatioThreshold[i] = 0;
		}
		for (int i=0; i<(int)m_InvariantRatioThreshold.size(); i++)
		{
			if (!nonzeroThresh)
				nonzeroThresh = m_InvariantRatioThreshold[i]>0;			
			m_InvariantRatioThreshold[i] = 0;
		}
		m_bEnforceSameSampleType = false;
		if (nonzeroThresh)
		{
			std::cout<<"[WARNING] Could not create generator set with thresholds. Ignoring them. "<<std::endl;		
			if (StartGeneration()==0)
			{
				std::cout<<"[ERROR] Cannot create a generator set for triplet sampling"<<std::endl;
				assert(false);
			}
		}
		else
			return 0;
	}
	return (int)m_MapSequencesAfterThresholding.size();
}

void SampleSetMapGenerator::ScoreLastMap(double )
{
}

bool SampleSetMapGenerator::FillTripletAndTest(int A, int B, int C, int Ap, int Bp, int Cp,
											   std::vector<int> & potentialMap)
{
	bool valid = A!=B && A!=C && B!=C;	// exclude impossible mappings
	valid = valid && Ap!=Bp && Ap!=Cp && Bp!=Cp;	// exclude impossible mappings
	valid = valid && A < B && B < C;		// exclude redundant mappings
//	valid = valid && Ap < Bp && Bp < Cp;		// exclude redundant mappings	
	if (m_bEnforceSameSampleType && valid)
	{
//		std::cout<<m_pSet1->GetSample(A).GetSampleType().c_str()<<" vs ";
//		std::cout<<m_pSet2->GetSample(Ap).GetSampleType().c_str()<<std::endl;
		valid = valid && m_pSet1->GetSample(A).GetSampleType()==m_pSet2->GetSample(Ap).GetSampleType();
		valid = valid && m_pSet1->GetSample(B).GetSampleType()==m_pSet2->GetSample(Bp).GetSampleType();
		valid = valid && m_pSet1->GetSample(C).GetSampleType()==m_pSet2->GetSample(Cp).GetSampleType();
	}
	potentialMap.push_back(A);		potentialMap.push_back(Ap);				
	potentialMap.push_back(B);		potentialMap.push_back(Bp);				
	potentialMap.push_back(C);		potentialMap.push_back(Cp);
	
	if (!valid)
		return false;

	valid = valid && TestInvariantFunction(potentialMap);
	valid = valid && TestDistanceToOther(potentialMap);				
	valid = valid && TestMinDistanceThreshold(potentialMap);
	return valid;
}


// threshold minValue / maxValue > ratioThreshold
void SampleSetMapGenerator::AddInvariantFunctionThreshold(SurfaceFeature * feature1, 
														  SurfaceFeature * feature2,
														  double ratioThreshold)
{
	m_InvariantFunction1.push_back(feature1);
	m_InvariantFunction2.push_back(feature2);	
	m_InvariantRatioThreshold.push_back(ratioThreshold);
}

// threshold minValue / maxValue > ratioThreshold	
void SampleSetMapGenerator::AddDistancetoToOtherPoints(SurfaceDistance * distance1, double scale1, 
													   SurfaceDistance * distance2, double scale2,
													   double ratioThreshold)
{
	m_DistanceToOtherPoints1.push_back(distance1);
	m_DistanceToOtherPoints2.push_back(distance2);	
	if (scale1==-1)
		scale1 = 1. / m_fDistanceNormalization1;
	m_DistanceToOtherPointsScale1.push_back(scale1);
	if (scale2==-1)
		scale2 = 1. / m_fDistanceNormalization2;
	m_DistanceToOtherPointsScale2.push_back(scale2);	
	m_DistanceToOtherRatioThreshold.push_back(ratioThreshold);
	
}

// threshold D(sample1, sample2) < minThreshold
void SampleSetMapGenerator::AddMinDistanceThreshold(SurfaceDistance * distance1, double scale1, 
													SurfaceDistance * distance2, double scale2,
													double minThreshold)
{	
	m_MinDistances1.push_back(distance1);
	m_MinDistances2.push_back(distance2);	
	if (scale1==-1)
		scale1 = sqrt(3.14159265) / m_fDistanceNormalization1;
	m_MinDistanceScale1.push_back(scale1);
	if (scale2==-1)
		scale2 = sqrt(3.14159265) / m_fDistanceNormalization2;
	m_MinDistanceScale2.push_back(scale2);	
	m_MinDistanceThresholds.push_back(sqrt(minThreshold));
}


bool SampleSetMapGenerator::TestInvariantFunction(std::vector<int> & candidateMap)
{
	for (int f=0; f<(int)m_InvariantFunction1.size(); f++)
	{
		SurfaceFeature * feature1 = m_InvariantFunction1[f];
		SurfaceFeature * feature2 = m_InvariantFunction2[f];	
		double ratioThresh = m_InvariantRatioThreshold[f];
		for (int i=0; i<(int)candidateMap.size(); i+=2)
		{
			double val1 = feature1->Value(m_pSet1->GetSample(candidateMap[i+0]));
			double val2 = feature2->Value(m_pSet2->GetSample(candidateMap[i+1]));
			double ratio = val1 / val2;
			if (val2 < val1)
				ratio = val2 / val1;
			if (ratio < ratioThresh)
				return false;
		}
	}
	return true;
}

bool SampleSetMapGenerator::TestDistanceToOther(std::vector<int> & candidateMap)
{	
	for (int d=0; d<(int)m_DistanceToOtherPoints1.size(); d++)
	{
		SurfaceDistance * distance1 = m_DistanceToOtherPoints1[d];
		SurfaceDistance * distance2 = m_DistanceToOtherPoints2[d];
		double scale1 = m_DistanceToOtherPointsScale1[d];
		double scale2 = m_DistanceToOtherPointsScale2[d];
		double distanceRatioThresh = m_DistanceToOtherRatioThreshold[d];
		
		for (int i=0; i<(int)candidateMap.size(); i+=2)
		{
			SurfaceSample s1 = m_pSet1->GetSample(candidateMap[i+0]);
			SurfaceSample s2 = m_pSet2->GetSample(candidateMap[i+1]);			
			for (int j=0; j<(int)candidateMap.size(); j+=2)
			{
				if (i==j)
					continue;
				
				SurfaceSample corrS1 = m_pSet1->GetSample(candidateMap[j+0]);
				SurfaceSample corrS2 = m_pSet2->GetSample(candidateMap[j+1]);			
				
				double d1 = distance1->Distance(s1, corrS1) * scale1;
				double d2 = distance2->Distance(s2, corrS2) * scale2;				
				double ratio = d1/d2;
				if (d1>d2)
					ratio = d2/d1;
				if (ratio < distanceRatioThresh)
					return false;
			}
		}
	}
	return true;
}

bool SampleSetMapGenerator::TestMinDistanceThreshold(std::vector<int> & candidateMap)
{		
	for (int d=0; d<(int)m_MinDistances1.size(); d++)
	{
		SurfaceDistance * distance1 = m_MinDistances1[d];
		SurfaceDistance * distance2 = m_MinDistances2[d];		
		double scale1 = m_MinDistanceScale1[d];
		double scale2 = m_MinDistanceScale2[d];		
		double minDistThreshold = m_MinDistanceThresholds[d];

		for (int i=0; i<(int)candidateMap.size()/2; i++)
		{
			for (int j=i+1; j<(int)candidateMap.size()/2; j++)
			{
				if (distance1->Distance(m_pSet1->GetSample(candidateMap[2*i+0]), 
										m_pSet1->GetSample(candidateMap[2*j+0]))*scale1 < minDistThreshold)
					return false;
				
				if (distance2->Distance(m_pSet2->GetSample(candidateMap[2*i+1]), 
										m_pSet2->GetSample(candidateMap[2*j+1]))*scale2 < minDistThreshold)
					return false;
			}
				
		}
	}
	return true;
}


void SampleSetMapGenerator::GetProgress(int & currProg, int & totalNum)
{
	totalNum = (int)m_MapSequencesAfterThresholding.size();
	currProg = m_iCurrentMapSequence;
}


/////////////////// SYMMETRIC ////////////////
MapGeneratorSymmetric::MapGeneratorSymmetric(SurfaceSampleSet * set, int numPntsPerMap, int maxNumMaps, 
											 GenerationMethod method, double surfaceArea)
: SampleSetMapGenerator(set, set, numPntsPerMap, maxNumMaps, method, surfaceArea, surfaceArea)
{
}

int MapGeneratorSymmetric::StartGeneration()
{
	m_iCurrentMapSequence = 0;
	int attempts=0;
	if (m_GenSearchMethod==GEN_RANDOM_SEARCH)
	{
		std::vector<int> potentialMap;
		while ((int)m_MapSequencesAfterThresholding.size()<m_iMaxNumMaps
			&& attempts++ < m_iMaxNumMaps*100)
		{
			int A = rand()%m_pSet1->NumSamples();
			int B = rand()%m_pSet1->NumSamples();
			int C = rand()%m_pSet1->NumSamples();
			
			potentialMap.clear();
			if (m_iNumPointsPerMap==3)
			{
				if (FillTripletAndTest(A, B, C, potentialMap))
					m_MapSequencesAfterThresholding.push_back(potentialMap);	
			}
			else if (m_iNumPointsPerMap==4)
			{
				int D = rand()%m_pSet1->NumSamples();				
				if (FillQuadrupletAndTest(A, B, C, D, potentialMap))
					m_MapSequencesAfterThresholding.push_back(potentialMap);
			}
		}
	}
	else if (m_GenSearchMethod==GEN_EXHAUSTIVE_SEARCH)
	{
		std::vector<int> potentialMap;		
		for (int A=0; A<m_pSet1->NumSamples(); A++)
			for (int B=0; B<m_pSet1->NumSamples(); B++)
				for (int C=0; C<m_pSet1->NumSamples(); C++)
					if (m_iNumPointsPerMap==3)
					{
						potentialMap.clear();
						if (FillTripletAndTest(A, B, C, potentialMap))
							m_MapSequencesAfterThresholding.push_back(potentialMap);
						if ((int)m_MapSequencesAfterThresholding.size()>=m_iMaxNumMaps)
							return (int)m_MapSequencesAfterThresholding.size();
					}
					else if (m_iNumPointsPerMap==4)
					{
						for (int D=0; D<m_pSet1->NumSamples(); D++)
						{
							potentialMap.clear();
							if (FillQuadrupletAndTest(A, B, C, D, potentialMap))
								m_MapSequencesAfterThresholding.push_back(potentialMap);							
							if ((int)m_MapSequencesAfterThresholding.size()>=m_iMaxNumMaps)
								return (int)m_MapSequencesAfterThresholding.size();
						}
					}
	}
	
	if (m_MapSequencesAfterThresholding.size()==0)
	{
		for (int i=0; i<(int)m_MinDistanceThresholds.size(); i++)
			m_MinDistanceThresholds[i] = 0;
		for (int i=0; i<(int)m_DistanceToOtherRatioThreshold.size(); i++)
			m_DistanceToOtherRatioThreshold[i] = 0;
		for (int i=0; i<(int)m_InvariantRatioThreshold.size(); i++)
			m_InvariantRatioThreshold[i] = 0;
		
		if (StartGeneration()==0)
		{
			std::cout<<"[ERROR] Cannot create maps"<<std::endl;
			assert(false);
		}
	}
	return (int)m_MapSequencesAfterThresholding.size();
}

bool MapGeneratorSymmetric::FillTripletAndTest(int A, int B, int C, 
											   std::vector<int> & potentialMap)
{
	bool valid = A!=B && A!=C && B!=C;	// exclude impossible mappings
	valid = valid && B<C;		// exclude redundant mappings
	potentialMap.push_back(A);		potentialMap.push_back(A);				
	potentialMap.push_back(B);		potentialMap.push_back(C);				
	potentialMap.push_back(C);		potentialMap.push_back(B);
	if (!valid)
		return false;
	valid = valid && TestInvariantFunction(potentialMap);
	valid = valid && TestDistanceToOther(potentialMap);				
	valid = valid && TestMinDistanceThreshold(potentialMap);	
	return valid;
}

bool MapGeneratorSymmetric::FillQuadrupletAndTest(int A, int B, int C, int D,
												  std::vector<int> & potentialMap)
{
	bool valid = A!=B && A!=C && A!=D && B!=C && B!=D && C!=D;	// exclude impossible mappings
	valid = valid && A<B && C<D;		// exclude redundant
	potentialMap.push_back(A);		potentialMap.push_back(B);
	potentialMap.push_back(B);		potentialMap.push_back(A);
	potentialMap.push_back(C);		potentialMap.push_back(D);
	potentialMap.push_back(D);		potentialMap.push_back(C);
	valid = valid && TestInvariantFunction(potentialMap);
	valid = valid && TestDistanceToOther(potentialMap);				
	valid = valid && TestMinDistanceThreshold(potentialMap);	
	return valid;
}


void MapGeneratorSymmetric::AddInvariantFunctionThreshold(SurfaceFeature * feature, double ratioThreshold)
{
	((SampleSetMapGenerator*)this)->AddInvariantFunctionThreshold(feature, feature, ratioThreshold);
}

void MapGeneratorSymmetric::AddDistancetoToOtherPoints(SurfaceDistance * distance, double ratioThreshold)
{
	((SampleSetMapGenerator*)this)->AddDistancetoToOtherPoints(distance, -1, distance, -1, ratioThreshold);
}

void MapGeneratorSymmetric::AddMinDistanceThreshold(SurfaceDistance * distance, double minDistance)
{
	((SampleSetMapGenerator*)this)->AddMinDistanceThreshold(distance, -1, distance, -1, minDistance);
}


//////////// MapGeneratorAdaptive ///////////////
MapGeneratorAdaptive::MapGeneratorAdaptive(SampleSetMapGenerator * otherMap, int maxNumMaps)
{
	m_pDefaultGenerator = otherMap;
	m_bSimpleExhaustive = true;
	m_iMaxNumMaps = maxNumMaps;
}

int MapGeneratorAdaptive::StartGeneration()
{
	int numMaps = m_pDefaultGenerator->m_MapSequencesAfterThresholding.size();

	if (numMaps > m_iMaxNumMaps)
	{
		m_bSimpleExhaustive = false;
		for (int i=0; i<numMaps; i++)
			m_ExploredMap.push_back(false);
		m_iNumMaps = 0;
	}
	else
	{
		m_bSimpleExhaustive = true;
	}
	return numMaps;
}

std::vector<int> & MapGeneratorAdaptive::GenerateNextMap(bool * valid)
{
//	std::cout<<"Generating : "<<m_bSimpleExhaustive<<" num="<<m_iNumMaps<<" vs "<<m_iMaxNumMaps<<std::endl;
	if (m_bSimpleExhaustive)
		return m_pDefaultGenerator->GenerateNextMap(valid);
	else
	{
		if (m_iNumMaps >= m_iMaxNumMaps)
		{
			*valid = false;
			return m_Empty;			
		}
		
		*valid = true;
		int randVal = rand();
		if (randVal<0)
			randVal *= -1;
		randVal = randVal % (int)m_ExploredMap.size();
		while (m_ExploredMap[randVal])
			randVal = (randVal+1) % m_ExploredMap.size();
		m_ExploredMap[randVal] = true;
		m_iNumMaps++;
		return m_pDefaultGenerator->m_MapSequencesAfterThresholding[randVal];
	}
}

void MapGeneratorAdaptive::ScoreLastMap(double score)
{
	if (!m_bSimpleExhaustive)
	{
		//todo
	}
	assert(false);
}

void MapGeneratorAdaptive::GetProgress(int & currProg, int & totalNum)
{
	if (m_bSimpleExhaustive)
	{
		m_pDefaultGenerator->GetProgress(currProg, totalNum);
	}
	else
	{
		currProg = m_iNumMaps;
		totalNum = m_iMaxNumMaps;
	}
}










