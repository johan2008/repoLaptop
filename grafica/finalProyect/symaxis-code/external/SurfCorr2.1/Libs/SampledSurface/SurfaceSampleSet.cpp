#include "SurfaceSampleSet.h"
#include "SurfaceFeature.h"
#include "FeatureMGD.h"
#include "FeatureAGD.h"
#include "SurfaceDistance.h"
#include <iostream>
#include <sstream>

SurfaceSampleSet::SurfaceSampleSet()
{
}

SurfaceSampleSet::SurfaceSampleSet(R3Mesh * mesh, const std::vector<int> & vertexIDs)
{
	if (vertexIDs.size()==0)
		for (int i=0; i<mesh->NVertices(); i++)
			AddSample(SurfaceSample(i, mesh));
	else
		for (int i=0; i<(int)vertexIDs.size(); i++)
			AddSample(SurfaceSample(vertexIDs[i], mesh));
}

SurfaceSampleSet::SurfaceSampleSet(SurfaceFeature * takeExtremesOfFeature, double ring,
								   int minNumPnts, int maxNumPnts,
								   ExtremaTypes extremaTypes)
{
	SurfaceSampleSet maxSet;
	SurfaceSampleSet minSet;
	takeExtremesOfFeature->Extremas(&minSet, &maxSet, ring);
	if (minSet.NumSamples()+maxSet.NumSamples() < minNumPnts)
	{
		maxSet.m_SurfaceSamples.clear();
		minSet.m_SurfaceSamples.clear();	
		ring /= 2.;
		takeExtremesOfFeature->Extremas(&minSet, &maxSet, ring);
	}
	switch(extremaTypes)
	{
		case FEAT_MAXIMA_ONLY:
			Union(maxSet);
			break;
		case FEAT_MINIMA_ONLY:
			Union(minSet);			
			break;
		case FEAT_ALL_EXTREMA:
			Union(maxSet);
			Union(minSet);
			break;
		case FEAT_GLOBAL_MAXIMA:
			assert(maxSet.NumSamples()>0);
			AddSample(maxSet.GetSample(0));
			break;
		case FEAT_GLOBAL_MINIMA:
			assert(minSet.NumSamples()>0);			
			AddSample(minSet.GetSample(0));			
			break;
	}	
	if (NumSamples()<minNumPnts)
	{
		Union(minSet);
		Union(maxSet);
	}
	if (NumSamples()>maxNumPnts)
	{
		while (NumSamples()>maxNumPnts)
			m_SurfaceSamples.erase((m_SurfaceSamples.begin()+((int)m_SurfaceSamples.size()-1)));
	}
	//if (NumSamples() < minNumPnts || NumSamples() > maxNumPnts)
	//{
	//	std::cout<<"[WARNING] Incorrect number of samples: ";
	//	std::cout<<minNumPnts<<" <= "<<NumSamples()<<" <= "<<maxNumPnts<<std::endl;
	//	std::cout<<"MinSet = "<<minSet.NumSamples()<<" MaxSet = "<<maxSet.NumSamples()<<std::endl;
	//	assert(NumSamples()>=minNumPnts && NumSamples()<=maxNumPnts);
	//}
	assert(NumSamples() >= 3);	
}

void SurfaceSampleSet::ClearSamples()
{
	m_SurfaceSamples.clear();
}

void SurfaceSampleSet::WriteSet(const VKString & setFile)
{
	std::ofstream dataStream(setFile.c_str(), std::ios::out);
	if (!dataStream.is_open())
	{
		std::cout<<"[ERROR] Could not write (sample set) file: "<<setFile.c_str()<<std::endl;
		assert(dataStream.is_open());
	}
	WriteSet(dataStream);
	dataStream.close();		
}

bool SurfaceSampleSet::LoadSet(const VKString & setFile, R3Mesh * mesh)
{
	std::ifstream dataStream(setFile.c_str());	
	if (!dataStream.is_open())
		return false;
	LoadSet(dataStream, mesh);
	dataStream.close();	
	return true;
}

void SurfaceSampleSet::WriteSet(std::ofstream & dataStream)
{
	dataStream<<NumSamples()<<"\n";
	for (int i=0; i<NumSamples(); i++)
	{
		SurfaceSample samp = GetSample(i);
		dataStream<<samp.GetSampleType().c_str()<<" "<<samp.TriID()<<" ";
		dataStream<<samp.B(0)<<" "<<samp.B(1)<<" "<<samp.B(2)<<"\n";
	}
}

void SurfaceSampleSet::LoadSet(std::ifstream & dataStream, R3Mesh * mesh)
{
	int numSamples;
	dataStream>>numSamples;
	for (int i=0; i<numSamples; i++)
	{
		VKString type="";
		type.readToken(dataStream);
		int triID;
		double b1, b2, b3;
		dataStream>>triID>>b1>>b2>>b3;
		AddSample(SurfaceSample(triID, b1, b2, b3, mesh, type));		
	}
}

SurfaceSampleSet * SurfaceSampleSet::GetCopyAnotherMesh(R3Mesh * meshCopy)
{
	SurfaceSampleSet * newSet = new SurfaceSampleSet();
	for (int i=0; i<NumSamples(); i++)
	{
		SurfaceSample samp = GetSample(i);
		samp.m_pMesh = meshCopy;
		newSet->AddSample(samp);
	}
	return newSet;	
}

SurfaceSampleSet * SurfaceSampleSet::GetSamples(VKString type)
{
	SurfaceSampleSet * newSet = new SurfaceSampleSet();
	for (int i=0; i<NumSamples(); i++)
		if (GetSample(i).GetSampleType()==type)
			newSet->AddSample(GetSample(i));
	return newSet;
}

const SurfaceSample & SurfaceSampleSet::GetSample(int sampleID) const
{
	if (sampleID >= NumSamples() || sampleID<0)
	{
		std::cout<<"[ERROR] Cannot get sample with ID = "<<sampleID;
		std::cout<<" numsamples="<<NumSamples()<<std::endl;
		assert(false);
	}
	return m_SurfaceSamples[sampleID];
}

int SurfaceSampleSet::AddSample(const SurfaceSample & sample)
{
	m_SurfaceSamples.push_back(sample);
	return ((int)m_SurfaceSamples.size()-1);
}

void SurfaceSampleSet::ReplaceSample(int id, const SurfaceSample & newSample)
{
	assert(id < (int)m_SurfaceSamples.size());
	m_SurfaceSamples[id] = newSample;
}

void SurfaceSampleSet::SetSamplesType(const VKString & type)
{
	for (int i=0; i<(int)m_SurfaceSamples.size(); i++)
		m_SurfaceSamples[i].SetSampleType(type); 
}

int SurfaceSampleSet::NumSamples() const
{
	return (int)m_SurfaceSamples.size();
}

SurfaceSampleSet & SurfaceSampleSet::Union(const SurfaceSampleSet & anotherSet)
{
	int origNumSamples = NumSamples();
	for (int i=0; i<anotherSet.NumSamples(); i++)
	{
		bool exists = false;
		for (int j=0; j<origNumSamples && !exists; j++)
			if (GetSample(j)==anotherSet.GetSample(i))
				exists = true;
			
		if (!exists)
			AddSample(anotherSet.GetSample(i));
	}
	return *this;
}

SurfaceSampleSet & SurfaceSampleSet::Union(const SurfaceSampleSet & anotherSet, double distThreshold, 
										   SurfaceDistance * distance, double distNorm)
{
	if (distance==NULL || distThreshold==0)
		return Union(anotherSet);
	
	int origNumSamples = NumSamples();
	for (int i=0; i<anotherSet.NumSamples(); i++)
	{
		bool exists = false;
		for (int j=0; j<origNumSamples && !exists; j++)
		{
			double dist = distance->Distance(GetSample(j), anotherSet.GetSample(i)) / distNorm;
			if (dist < distThreshold)
				exists = true;
		}
		
		if (!exists)
			AddSample(anotherSet.GetSample(i));
	}
	return *this;
}

bool SurfaceSampleSet::ContainsExact(const SurfaceSample & sample, int * localID)
{
	if (localID!=NULL)
		*localID=-1;
	for (int i=0; i<NumSamples(); i++)
		if (GetSample(i)==sample)
		{
			if (localID!=NULL)
				*localID = i;
			return true;
		}
			
	return false;
}
