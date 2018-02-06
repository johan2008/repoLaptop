#include "gaps.h"
#include "SurfaceSample.h"
#include <vector>

#ifndef __SURFACE_SAMPLE_SET_H
#define __SURFACE_SAMPLE_SET_H

class SurfaceFeature;
class FeatureMGD;
class SurfaceDistance;

class SurfaceSampleSet
{
	public:
		enum ExtremaTypes
		{
			FEAT_MAXIMA_ONLY,
			FEAT_MINIMA_ONLY,
			FEAT_ALL_EXTREMA,
			FEAT_GLOBAL_MAXIMA,
			FEAT_GLOBAL_MINIMA,
		};
		SurfaceSampleSet();
		SurfaceSampleSet(SurfaceFeature * takeExtremesOfFeature, double ring,
						 int minNumPnts, int maxNumPnts,
						 ExtremaTypes extremaTypes);
		SurfaceSampleSet(R3Mesh * mesh, const std::vector<int> & vertexIDs=std::vector<int>());	
		virtual ~SurfaceSampleSet(){}
		virtual void ClearSamples();
	
		virtual void WriteSet(const VKString & setFile);
		virtual bool LoadSet(const VKString & setFile, R3Mesh * mesh);
		virtual void WriteSet(std::ofstream & textStream);
		virtual void LoadSet(std::ifstream & textStream, R3Mesh * mesh);
	
		SurfaceSampleSet * GetCopyAnotherMesh(R3Mesh * meshCopy);
		SurfaceSampleSet * GetSamples(VKString type);
		void SetSamplesType(const VKString & type);
	
		bool ContainsExact(const SurfaceSample & sample, int * localID = NULL);
		//TODO: approximate w/ threshold and speedup

		const SurfaceSample & GetSample(int sampleID) const;
		int AddSample(const SurfaceSample & sample);
		void ReplaceSample(int id, const SurfaceSample & newSample);
	
		int NumSamples() const;
		SurfaceSampleSet & Union(const SurfaceSampleSet & anotherSet);
		SurfaceSampleSet & Union(const SurfaceSampleSet & anotherSet, double distThreshold, 
								 SurfaceDistance * distance, double distNorm);
	
	protected:
		std::vector<SurfaceSample> m_SurfaceSamples;
};

#endif
