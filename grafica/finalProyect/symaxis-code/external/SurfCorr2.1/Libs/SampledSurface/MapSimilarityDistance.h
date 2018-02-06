#ifndef __MAP_SIMILARITY_DISTANCE_H
#define __MAP_SIMILARITY_DISTANCE_H

#include "SurfaceMapSimilarity.h"

/**
 * Inverse of sum of geodesic distances between mapped points
 */
class MapSimilarityDistance : public SurfaceMapSimilarity
{
	public:
		enum MapSimInvDistApproxType
		{
			MAPSIM_FORWARD_ONLY,
			MAPSIM_BIDIR,
			MAPSIM_FAST_DIST_SUBSETS,
			MAPSIM_FAST_DIST_SUBSETS_BIDIR
		};
		
		enum SimilarityFlags
		{
			SIMFLAG_NONE = 0,
			SIMFLAG_NO_INCONSISTENT_GENS=1,
			SIMFLAG_SCALE_BY_CONFIDENCE=2
		};
		
		// to compare any maps
		MapSimilarityDistance(SurfaceMap * M1, SurfaceMap * M2,
							  MapSimInvDistApproxType type, int simFlags, double maxDistance, 
							  double similarityThreshold, double sigma,
							  const VKString & sampleSet1, const VKString & sampleSet2,
							  const VKString & distance1, const VKString & distance2);
		virtual ~MapSimilarityDistance();
		virtual VKString GetSurfaceMapSimilarityType();
	
		virtual SurfaceMapSimilarity * Clone(SurfaceMap * M1, SurfaceMap * M2);
	
		virtual double Similarity();

		virtual bool HaveInconsistentGenerators();	
		static int NumSharedCorrespondences(MapConformal * M1, MapConformal * M2);
		static double FastSimilarity(double sigma2, MapConformal * M1, MapConformal * M2,
									 const SurfaceSample & samp, SurfaceDistance * dist1,
									 double distNorm);
	
	protected:
	
		bool IsFlagSet(SimilarityFlags flag);
		
		MapSimInvDistApproxType m_ApproximationType;
		double m_fMaxDistance;
		
		int m_SimFlags;
	
		virtual double CalculateSimilarityAtSample(const SurfaceSample & samp);
		virtual double CalculateDissimilarityAtSample(const SurfaceSample & samp);	
	
		
		double GetDistFastError(SurfaceMap * M1, SurfaceMap * M2, bool forward, 
								std::vector<double> * perSampleErr=NULL);
		double GetOneDirError(SurfaceMap * M1, SurfaceMap * M2, bool forward, 
							  std::vector<double> * perSampleErr=NULL);
		
		void PrepareSamplesIteratorOneDir(SurfaceDistance ** dist, double * norm, double * maxDist, 
										  SurfaceSampleSet ** sampSet, SurfaceMap * M1, SurfaceMap * M2, 
										  bool forward);
		void PrepareSamplesIteratorDistFast(SurfaceDistance ** dist, double * norm, double * maxDist, 
											SurfaceSampleSet ** sampSet, SurfaceMap * M1, SurfaceMap * M2, 
											bool forward);
		SurfaceDistance * GetDistanceForSurface(SampledSurface * surf, const VKString & metricName);
		
		double GetSampleOneDirError(SurfaceMap * M1, SurfaceMap * M2, const SurfaceSample & samp, bool forward,
									SurfaceDistance * dist, double norm, double maxDist);
		double GetSampleDistFastError(SurfaceMap * M1, SurfaceMap * M2, const SurfaceSample & samp, bool forward,
									  SurfaceDistance * dist, double norm, double maxDist);
		
		double m_fSigma;
};

#endif
