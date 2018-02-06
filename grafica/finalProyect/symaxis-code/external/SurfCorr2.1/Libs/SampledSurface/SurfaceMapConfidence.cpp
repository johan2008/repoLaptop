#include "SurfaceMapConfidence.h"
#include "MapScoredCollection.h"
#include "SurfaceMapSimilarity.h"
#include "MapConfidenceMultiConf.h"
#include "MapMultiConformal.h"
#include "SurfaceMap.h"
#include "Sortable.h"

std::map<VKString, SurfaceMapConfidence*> SurfaceMapConfidence::m_mNameToSeed;

SurfaceMapConfidence::SurfaceMapConfidence(SurfaceMap * correspondingMap, 
										   bool smoothConfidence, bool smoothError,
										   double smoothingSigma, 
										   const VKString & testSamples1,
										   const VKString & testSamples2,
										   bool cacheConfidence, bool cacheError)
{
	m_bCacheConfidence = cacheConfidence;
	m_bCacheError = cacheError;
	
	m_SmoothSigma = smoothingSigma;
	m_TotalConfidence = -1;
	m_TotalError = -1;

	m_pSurfaceMap = correspondingMap;
	m_SampleSet1 = testSamples1;
	m_SampleSet2 = testSamples2;	
	m_bSmoothConfidence = smoothConfidence;
	m_bSmoothError = smoothError;
	if (m_bSmoothConfidence)
		assert(m_bCacheConfidence);
	if (m_bSmoothError)
		assert(m_bCacheError);
	
	if (m_bCacheConfidence && m_pSurfaceMap!=NULL)
	{
		for (int i=0; i<2; i++)
		{
			SurfacePerVertexValue * values = new SurfacePerVertexValue(m_pSurfaceMap->GetSurface(i), true);
			m_CachedConfidenceOnSurf[m_pSurfaceMap->GetSurface(i)] = values;
		}
	}
	
	if (m_bCacheError && m_pSurfaceMap!=NULL)
	{
		for (int i=0; i<2; i++)
		{
			SurfacePerVertexValue * values = new SurfacePerVertexValue(m_pSurfaceMap->GetSurface(i), true);
			m_CachedErrorOnSurf[m_pSurfaceMap->GetSurface(i)] = values;
		}
	}	
}

SurfaceMapConfidence::~SurfaceMapConfidence()
{
	if (m_bCacheConfidence)
	{	
		for (int i=0; i<2; i++)
			delete m_CachedConfidenceOnSurf[m_pSurfaceMap->GetSurface(i)];
	}
	
	if (m_bCacheError)
	{
		for (int i=0; i<2; i++)
			delete m_CachedErrorOnSurf[m_pSurfaceMap->GetSurface(0)];
	}
}

MapCoarse * SurfaceMapConfidence::GetBestSparseCorrs(SurfaceSampleSet * sampleSetFrom, 
													 int numCorrs, SurfaceSampleSet ** bestSparse)
{
	assert(sampleSetFrom->NumSamples() >= numCorrs);
	std::vector<Sortable> scoredSamples;
	for (int i=0; i<sampleSetFrom->NumSamples(); i++)
	{
		double value = ConfidenceAtSample(sampleSetFrom->GetSample(i));
		scoredSamples.push_back(Sortable(value, NULL, i));
	}
	
	std::sort(scoredSamples.begin(), scoredSamples.end());
	
	int tail = (int)scoredSamples.size()-1;
	*bestSparse = new SurfaceSampleSet();
	while ((**bestSparse).NumSamples() < numCorrs)
	{
		(**bestSparse).AddSample(sampleSetFrom->GetSample(scoredSamples[tail].id));
		tail--;
		assert(tail>=0);
	}
	
	return new MapCoarse(*bestSparse, NULL, m_pSurfaceMap, MapCoarse::F2C_FORWARD_ONLY);
}

bool SurfaceMapConfidence::DisplayProgressBar()
{
	return false;
}

void SurfaceMapConfidence::ClearConfidenceValue()
{
	m_TotalConfidence = -1;
	m_TotalError = -1;
	for (int i=0; i<2; i++)
		if (m_CachedConfidenceOnSurf.find(m_pSurfaceMap->GetSurface(i))!=m_CachedConfidenceOnSurf.end())
			m_CachedConfidenceOnSurf[m_pSurfaceMap->GetSurface(i)]->ClearValues();
		
	for (int i=0; i<2; i++)
		if (m_CachedErrorOnSurf.find(m_pSurfaceMap->GetSurface(i))!=m_CachedErrorOnSurf.end())
			m_CachedErrorOnSurf[m_pSurfaceMap->GetSurface(i)]->ClearValues();
}

double SurfaceMapConfidence::Confidence()
{
	if (m_TotalConfidence==-1)
		m_TotalConfidence = CalculateTotalConfidence();
	return m_TotalConfidence;
}

double SurfaceMapConfidence::Error()
{
	if (m_TotalError==-1)
		m_TotalError = CalculateTotalError();
	return m_TotalError;
}

void SurfaceMapConfidence::GetPerVertexConfidence(SampledSurface * surf, 
													std::vector<double> & vals)
{	
//	std::cout<<"Getting per vertex confidence "<<std::flush;
	m_CachedConfidenceOnSurf[surf]->GetPerVertexValues(vals, false);
//	std::cout<<" - DONE"<<std::endl;	
}

bool SurfaceMapConfidence::GetPerVertexErrors(SampledSurface * surf, 
												std::vector<double> & vals, 
												bool normalize)
{
	return m_CachedErrorOnSurf[surf]->GetPerVertexValues(vals, -FLT_MAX, FLT_MAX, normalize);
}

void SurfaceMapConfidence::GetPerSampleConfidence(SampledSurface * surf, 
												  std::vector<double> & vals)
{
	SurfaceSampleSet * sampleSet = NULL;
	if (m_pSurfaceMap->GetSurface(0)==surf)
		sampleSet = surf->GetSampleSet(m_SampleSet1);
	else
		sampleSet = surf->GetSampleSet(m_SampleSet2);
	
	for (int i=0; i<sampleSet->NumSamples(); i++)
		vals.push_back(ConfidenceAtSample(sampleSet->GetSample(i)));
}

void SurfaceMapConfidence::GetPerSampleErrors(SampledSurface * surf, 
											  std::vector<double> & vals, 
											  bool normalize)
{
	SurfaceSampleSet * sampleSet = NULL;
	if (m_pSurfaceMap->GetSurface(0)==surf)
		sampleSet = surf->GetSampleSet(m_SampleSet1);
	else
		sampleSet = surf->GetSampleSet(m_SampleSet2);
	
	double norm=-FLT_MAX;
	double maxErr= ErrorThreshold(surf);
	for (int i=0; i<sampleSet->NumSamples(); i++)
	{
		double err = ErrorAtSample(sampleSet->GetSample(i));
		if (normalize && maxErr>0)
			err = err > maxErr ? 1. : (err / maxErr);
		else if (normalize)
			norm = vkMax(norm, err);
		
		vals.push_back(err);
	}
	
	if (normalize && maxErr==-1)
	{
		for (int i=0; i<(int)vals.size(); i++)
			vals[i]/=norm;
	}
}

double SurfaceMapConfidence::ErrorThreshold(SampledSurface * )
{
	return -1;
}

double SurfaceMapConfidence::ConfidenceAtSampleNoSmoothing(const SurfaceSample & sample)
{
	bool disableSmoothing = m_bSmoothConfidence;
	m_bSmoothConfidence = false;
	double retVal = ConfidenceAtSample(sample);
	m_bSmoothConfidence = disableSmoothing;
	return retVal;
}

double SurfaceMapConfidence::ConfidenceAtSample(const SurfaceSample & sample)
{
	if (m_bSmoothConfidence )
	{
		PrecomputeSmoothConfidence();
	}

	bool atVertex;
	sample.NearestVertex(&atVertex);
	SampledSurface * surface = m_pSurfaceMap->GetSurface(sample);
	R3Mesh * mesh = surface->GetMesh();

	if (!m_bCacheConfidence)
	{
		if (atVertex)
			return CalculateConfidenceAtSample(sample);
		else
			return sample.Interpolate(CalculateConfidenceAtSample(SurfaceSample(sample.VertID(0), mesh)),
									  CalculateConfidenceAtSample(SurfaceSample(sample.VertID(1), mesh)),
									  CalculateConfidenceAtSample(SurfaceSample(sample.VertID(2), mesh)));
	}
	
	
	SurfaceSample getValues = m_CachedConfidenceOnSurf[surface]->NextAtVertexSampleRequiredForInterpolation(sample);
	while (!getValues.Invalid())
	{
		m_CachedConfidenceOnSurf[surface]->SetVertexValue(getValues,
														  CalculateConfidenceAtSample(getValues));
		getValues = m_CachedConfidenceOnSurf[surface]->NextAtVertexSampleRequiredForInterpolation(sample);
	}
	
	return m_CachedConfidenceOnSurf[surface]->GetInterpolatedValue(sample);
}

double SurfaceMapConfidence::ErrorAtSample(const SurfaceSample & sample)
{	
	if (m_bSmoothError)
		PrecomputeSmoothError();
	bool atVertex;
	sample.NearestVertex(&atVertex);
	SampledSurface * surface = m_pSurfaceMap->GetSurface(sample);
	R3Mesh * mesh = surface->GetMesh();
	
	if (!m_bCacheError)
	{
		if (atVertex)
			return CalculateErrorAtSample(sample);
		else
			return sample.Interpolate(CalculateErrorAtSample(SurfaceSample(sample.VertID(0), mesh)),
									  CalculateErrorAtSample(SurfaceSample(sample.VertID(1), mesh)),
									  CalculateErrorAtSample(SurfaceSample(sample.VertID(2), mesh)));
	}
	
	SurfaceSample getValues = m_CachedErrorOnSurf[surface]->NextAtVertexSampleRequiredForInterpolation(sample);
	while (!getValues.Invalid())
	{
		m_CachedErrorOnSurf[surface]->SetVertexValue(getValues,
											   CalculateErrorAtSample(getValues));
		getValues = m_CachedErrorOnSurf[surface]->NextAtVertexSampleRequiredForInterpolation(sample);
	}
	
	return m_CachedErrorOnSurf[surface]->GetInterpolatedValue(sample);
}

void SurfaceMapConfidence::PrecomputeSmoothConfidence()
{
	m_bSmoothConfidence = false;

	if (m_SampleSet1!="none")	
		PrecomputeSmoothConfidence(m_pSurfaceMap->GetSurface(0));
	if (m_SampleSet2!="none")
		PrecomputeSmoothConfidence(m_pSurfaceMap->GetSurface(1));	
}

void SurfaceMapConfidence::PrecomputeSmoothConfidence(SampledSurface * surface)
{
	VKString sampleSetName;
	SurfaceSampleSet * sampleSet = NULL;
	if (surface==m_pSurfaceMap->GetSurface(0))
		sampleSetName = m_SampleSet1;
	else
	{
		assert(surface==m_pSurfaceMap->GetSurface(1));
		sampleSetName = m_SampleSet2;
	}
	
	if (sampleSetName=="none")
	{
		std::cout<<"[ERROR] Cannot smooth confidence, sample sets: "<<sampleSetName.c_str()<<" in [";
		std::cout<<m_SampleSet1.c_str()<<" "<<m_SampleSet2.c_str()<<"]"<<std::endl;
		assert(sampleSetName!="none");
	}
		
	sampleSet = surface->GetSampleSet(sampleSetName);
	assert(sampleSet!=NULL);
	
	//std::cout<<"[WARNING] Precomputing smooth confidence ("<<sampleSet->NumSamples()<<" / ";
	//std::cout<<surface->GetMesh()->NVertices()<<") : "<<std::flush;
	
	std::vector<double> samplesInSet;
	for (int i=0; i<sampleSet->NumSamples(); i++)
	{
		//if (i%100==0)
		//	std::cout<<i<<" "<<std::flush;
		
		const SurfaceSample & currSamp = sampleSet->GetSample(i);
		if (m_CachedConfidenceOnSurf[surface]->NextVertexRequiredForInterpolation(currSamp)!=-1)
			m_CachedConfidenceOnSurf[surface]->SetVertexValue(currSamp, 
															  CalculateConfidenceAtSample(currSamp));
	}
	//std::cout<<std::endl;
	
	m_CachedConfidenceOnSurf[surface]->SmoothGaussianFromKnown(m_SmoothSigma);
	
}

void SurfaceMapConfidence::PrecomputeSmoothError()
{	
	m_bSmoothError = false;
	PrecomputeSmoothError(m_pSurfaceMap->GetSurface(0));
	PrecomputeSmoothError(m_pSurfaceMap->GetSurface(1));	
}

void SurfaceMapConfidence::PrecomputeSmoothError(SampledSurface * surface)
{
	VKString sampleSetName;
	SurfaceSampleSet * sampleSet = NULL;
	if (surface==m_pSurfaceMap->GetSurface(0))
		sampleSetName = m_SampleSet1;
	else
		sampleSetName = m_SampleSet2;
	
	assert(sampleSetName!="none" && sampleSetName!="AllVertices");
	
	sampleSet = surface->GetSampleSet(sampleSetName);
	assert(sampleSet!=NULL);
	
//	std::cout<<"[WARNING] Precomputing smooth error ("<<sampleSet->NumSamples()<<" / ";
//	std::cout<<surface->GetMesh()->NVertices()<<") : "<<std::flush;
	
	std::vector<double> samplesInSet;
	for (int i=0; i<sampleSet->NumSamples(); i++)
	{
//		if (i%100==0)
//			std::cout<<i<<" "<<std::flush;
		
		const SurfaceSample & currSamp = sampleSet->GetSample(i);
		if (m_CachedErrorOnSurf[surface]->NextVertexRequiredForInterpolation(currSamp)!=-1)
			m_CachedErrorOnSurf[surface]->SetVertexValue(currSamp, 
														 CalculateErrorAtSample(currSamp));
	}
	
	m_CachedErrorOnSurf[surface]->SmoothGaussianFromKnown(sampleSetName, m_SmoothSigma);
}


double SurfaceMapConfidence::CalculateTotalConfidence()
{
	bool disableSmoothing = m_bSmoothConfidence;
	m_bSmoothConfidence = false;
	
	double totalConfidence=0;
	int finalNormalization=0;
	
	int MAX_EVALS=0;
	TimeProfiler profiler;
	if (DisplayProgressBar())
	{
		for (int i=0; i<2; i++)
		{
			VKString setname = (i==0) ? m_SampleSet1 : m_SampleSet2;
			if (setname!="none")
			{
				SurfaceSampleSet * sampleSet = m_pSurfaceMap->GetSurface(i)->GetSampleSet(setname);
				assert(sampleSet!=NULL);	
				MAX_EVALS+=sampleSet->NumSamples();
			}
		}
	}
	int currEvals=0;
	
	for (int i=0; i<2; i++)
	{
		VKString sampleSetName = (i==0) ? m_SampleSet1 : m_SampleSet2;
		
		if (sampleSetName!="none")
		{
			double locTotalConfidence=0;
			SurfaceSampleSet * sampleSet = m_pSurfaceMap->GetSurface(i)->GetSampleSet(sampleSetName);
			assert(sampleSet!=NULL);
			for (int j=0; j<sampleSet->NumSamples(); j++)
			{
				if (DisplayProgressBar())
					profiler.WriteProgress("ConfidenceCalculator", currEvals++, MAX_EVALS);
				
				locTotalConfidence += ConfidenceAtSample(sampleSet->GetSample(j));
			}
			totalConfidence += locTotalConfidence / sampleSet->NumSamples();
			finalNormalization++;
		}
	}	
	if (DisplayProgressBar())
		std::cout<<std::endl;
	
	
	m_bSmoothConfidence = disableSmoothing;
	assert(finalNormalization>=1);
	return totalConfidence / finalNormalization;
}

double SurfaceMapConfidence::CalculateTotalError()
{	
	bool disableSmoothing = m_bSmoothError;
	m_bSmoothError = false;
	
	double totalError=0;
	int finalNormalization=0;
	
	int MAX_EVALS=0;
	TimeProfiler profiler;
	if (DisplayProgressBar())
	{
		for (int i=0; i<2; i++)
		{
			VKString setname = (i==0) ? m_SampleSet1 : m_SampleSet2;
			if (setname!="none")
			{
				SurfaceSampleSet * sampleSet = m_pSurfaceMap->GetSurface(i)->GetSampleSet(setname);
				assert(sampleSet!=NULL);	
				MAX_EVALS+=sampleSet->NumSamples();
			}
		}
	}
	
	int currEvals=0;
	for (int i=0; i<2; i++)
	{
		VKString sampleSetName = (i==0) ? m_SampleSet1 : m_SampleSet2;
		if (sampleSetName!="none")
		{
			double locTotalError=0;
			SurfaceSampleSet * sampleSet = m_pSurfaceMap->GetSurface(i)->GetSampleSet(sampleSetName);
			assert(sampleSet!=NULL);
			int currNorm=0;
			for (int j=0; j<sampleSet->NumSamples(); j++)
			{
				if (DisplayProgressBar())
					profiler.WriteProgress("ConfidenceCalculator", currEvals++, MAX_EVALS);

				double err = ErrorAtSample(sampleSet->GetSample(j));
				if (err>=0)
				{
					locTotalError += err;
					currNorm++;
				}
			}
			totalError += locTotalError / currNorm;
			finalNormalization++;
		}
	}
	if (DisplayProgressBar())
		std::cout<<std::endl;
	
	assert(finalNormalization>0);
	m_bSmoothError = disableSmoothing;	
	return totalError / finalNormalization;	
}


double SurfaceMapConfidence::CalculateErrorAtSample(const SurfaceSample & sample)
{
	return 1. - ConfidenceAtSample(sample);
}

SurfaceMapConfidence * SurfaceMapConfidence::GetClonableSeed(const VKString & seedName)
{
	if (m_mNameToSeed.find(seedName)!=m_mNameToSeed.end())
		return m_mNameToSeed[seedName];

	return NULL;
}

void SurfaceMapConfidence::SetClonableSeed(const VKString & seedName, SurfaceMapConfidence * seed)
{
	m_mNameToSeed[seedName] = seed;
}

void SurfaceMapConfidence::SaveConfidenceValues(std::ofstream & textStream)
{	
	textStream<<"ConfidenceValues "<<m_TotalConfidence<<" "<<m_TotalError<<"\n";
	textStream<<"PerVertexConfidence "<<m_bCacheConfidence<<"\n";
	
	for (int surfID=0; surfID<2; surfID++)
	{
		SampledSurface * surface = m_pSurfaceMap->GetSurface(surfID);
		bool exists =  (m_CachedConfidenceOnSurf.find(surface)!=m_CachedConfidenceOnSurf.end());
		textStream<<"SurfaceExists "<<exists<<"\n";
		if (exists)
			m_CachedConfidenceOnSurf[surface]->WritePerVertexValues(textStream);
	}
	
	textStream<<"PerVertexError "<<m_bCacheError<<"\n";
	
	for (int surfID=0; surfID<2; surfID++)
	{
		SampledSurface * surface = m_pSurfaceMap->GetSurface(surfID);
		bool exists =  (m_CachedErrorOnSurf.find(surface)!=m_CachedErrorOnSurf.end());
		textStream<<"SurfaceExists "<<exists<<"\n";
		if (exists)
			m_CachedErrorOnSurf[surface]->WritePerVertexValues(textStream);
	}
}

void SurfaceMapConfidence::LoadConfidenceValues(std::ifstream & textStream)
{
	std::string tempStr;
	textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "ConfidenceValues")==0);
	textStream>>m_TotalConfidence>>m_TotalError;
	
	textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "PerVertexConfidence")==0);
	bool cacheConfidence;
	textStream>>cacheConfidence;
	assert(cacheConfidence==m_bCacheConfidence);
		
	for (int surfID=0; surfID<2; surfID++)
	{
		textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "SurfaceExists")==0);
		bool exists;			textStream>>exists;
		SampledSurface * surface = m_pSurfaceMap->GetSurface(surfID);

		if (exists)
		{
			if (m_CachedConfidenceOnSurf[surface]==NULL)
				m_CachedConfidenceOnSurf[surface] = new SurfacePerVertexValue(surface, true);
			bool ok = m_CachedConfidenceOnSurf[surface]->LoadPerVertexValues(textStream);
			assert(ok);
		}
	}
	
	textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "PerVertexError")==0);
	bool cacheError;
	textStream>>cacheError;
	assert(cacheError==m_bCacheError);
	
	for (int surfID=0; surfID<2; surfID++)
	{
		textStream>>tempStr;	assert(strcmp(tempStr.c_str(), "SurfaceExists")==0);
		bool exists;			textStream>>exists;
		SampledSurface * surface = m_pSurfaceMap->GetSurface(surfID);
		
		if (exists)
		{
			if (m_CachedErrorOnSurf[surface]==NULL)
				m_CachedErrorOnSurf[surface] = new SurfacePerVertexValue(surface, true);			
			bool ok = m_CachedErrorOnSurf[surface]->LoadPerVertexValues(textStream);
			assert(ok);
		}
	}	
}



////////////////////////// MapConfidenceConstant //////////////////////////
//MapConfidenceConstant::MapConfidenceConstant(SurfaceMap * correspondingMap)
//: SurfaceMapConfidence(correspondingMap, false, false, -1., "none", "none", false, false)
//{
//}
//
//double MapConfidenceConstant::Confidence()
//{
//	return 1.;
//}
//double MapConfidenceConstant::Error()
//{
//	return 0;
//}
//double MapConfidenceConstant::ConfidenceAtSampleNoSmoothing(const SurfaceSample & sample)
//{
//	return 1.;
//}
//
//double MapConfidenceConstant::ConfidenceAtSample(const SurfaceSample & sample)
//{
//	return 1.;
//}
//
//double MapConfidenceConstant::ErrorAtSample(const SurfaceSample & sample)
//{
//	return 0.;
//}
//
//VKString MapConfidenceConstant::GetMapConfidenceType()
//{
//	return "MapConfidenceConstant";
//}
//
//SurfaceMapConfidence * MapConfidenceConstant::Clone(SurfaceMap * currMap)
//{
//	return new MapConfidenceConstant(currMap);
//}



