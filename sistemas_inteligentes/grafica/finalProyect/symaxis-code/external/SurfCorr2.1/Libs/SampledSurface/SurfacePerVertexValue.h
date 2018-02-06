#ifndef __SURFACE_PER_VERTEX_VALUE_H
#define __SURFACE_PER_VERTEX_VALUE_H

#include "SurfaceSample.h"
#include "SurfaceSampleSet.h"
#include <map>

class SampledSurface;
class DistanceOnTheFly;

/**
 * This class caches per vertex values. 
 * It can interpolate, smooth, etc...
 * Used by: SurfaceFeature, SurfaceMapConfidence, SurfaceMapSimilarity
 */
class SurfacePerVertexValue
{
	public:
		SurfacePerVertexValue(SampledSurface * surface, bool allocateAsNeeded);
		~SurfacePerVertexValue();

		void ClearValues();
		void SetVertexValue(int vertexID, double value);
		void SetVertexValue(const SurfaceSample & atVertexSample, double value);
		double GetInterpolatedValue(const SurfaceSample & samp);
		SurfaceSample NextAtVertexSampleRequiredForInterpolation(const SurfaceSample & samp);
		int NextVertexRequiredForInterpolation(const SurfaceSample & samp);
			
		void SmoothSimpleLaplacian(int numIterations);
		void SmoothGaussianFromKnown(double sigma);	
		void SmoothGaussianFromKnown(const VKString & knownSet, double sigma);
		void SmoothGaussianFromKnown(SurfaceSampleSet * knownSet, double sigma);
		
		bool SmoothGaussianFromKnownBySplatting(SurfaceSampleSet * knownSet, double sigma, double splatRadius);
		void SmoothGaussianFromKnownExact(SurfaceSampleSet * knownSet, double sigma);
	
		void SmoothPickNearestNeighborsForUnassigned();
	
		bool GetPerVertexValues(std::vector<double> & saveValues, double setMin = -FLT_MAX, 
								double setMax = FLT_MAX, bool maxToOne=true);
	
		void LocalExtrema(SurfaceSampleSet * minSet, SurfaceSampleSet * maxSet, double minGeodesic);
	
		void CutToMaximum(double maxValue);
		void CutToMinimum(double minValue);	
	
		double MinValue();
		double MaxValue();
	
		void WritePerVertexValuesInArff(const VKString & filename,
										const VKString & relation, 
										const VKString & attributeName);
		void WritePerVertexValues(std::ofstream & textStream);
		bool LoadPerVertexValues(std::ifstream & textStream);
	
		bool AllValuesSet();
	
		DistanceOnTheFly * GetSurfaceDistance(double maxDistance=-1, 
											  const VKString & checkDistance="default");
	
	
		static double GetPatchArea(R3Mesh * mesh, R3MeshVertex * vertex);
		static double GetPatchArea(R3Mesh * mesh, int vertexID);

	protected:
	
		void AllocatePerVertexStorage();
		void AllocatePerVertexStorageFillFromKnown();
	
		struct VertexData 
		{
			R3MeshVertex *vertex;
			double distance;
			RNBoolean selected;
			VertexData **heappointer;
		};
	
		double ReestimateSigma(double sigma, int numSamples);
		double EstimateRadiusForSplatting(int numSamples);
	
		SampledSurface * m_pSurface;
	
		double m_fMinValue;
		double m_fMaxValue;
	
		bool m_bAllVertexValuesAllocated;
	
		double * m_PerVertexValues;
		bool * m_PerVertexValuesExist;
	
		int m_iNumValuesNotSet;
	
		std::map<int, double> m_VertexToValueMap;
	
		static void IsExtrema(R3Mesh * mesh, VertexData * vertex_data, bool & isMin, bool &isMax, 
							  R3MeshVertex * vOrg, double ring);
		
		static double * GetKnownValueArray(int id);
		static double * GetNormalizationArray(int id);	
		static std::map<int, double *> s_KnownValuesArrays;
		static std::map<int, double *> s_NormalizationArrays;
	
		struct PrecomputedWeights
		{	
			// note pointers to a row of data is stored
			void GetSampleWeighting(int vID, double ** savedWeights, 
									int ** savedIDs, int * savedNumPnts);
			double m_fSigma;
			VKString m_sSampleSet;
			double ** m_WeightArray;
			int ** m_SampleIDs;
			int * m_NumPointsInNhd;
			SampledSurface * m_pSurface;
		};
		// surface 
		static void ClearPrecomputedWeights();
		static PrecomputedWeights * GetInterpolationWeighting(const VKString & sampleSetName,
															  SampledSurface * surf, double sigma);
		static std::map<SampledSurface *, std::vector<PrecomputedWeights *> > s_InterpolationWeights;
};

#endif

