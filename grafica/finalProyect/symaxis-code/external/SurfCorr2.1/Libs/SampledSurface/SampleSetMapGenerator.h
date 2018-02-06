#include "gaps.h"
#include "VKStringList.h"
#include <vector>

#include "SurfaceSampleSet.h"
#include "SurfaceDistance.h"
#include "SurfaceFeature.h"

#ifndef __SAMPLE_SET_MAP_GENERATOR_H
#define __SAMPLE_SET_MAP_GENERATOR_H

/**
 * Explores potential mapps from set1 to set2.
 */ 
class SampleSetMapGenerator
{
	public:
		friend class MapGeneratorAdaptive;
	
		enum GenerationMethod
		{
			GEN_EXHAUSTIVE_SEARCH,
			GEN_RANDOM_SEARCH
		};
	
		SampleSetMapGenerator(SurfaceSampleSet * set1, SurfaceSampleSet * set2, 
							  int numPntsPerMap, int maxNumMaps, GenerationMethod method, 
							  double surface1Area, double surface2Area);
		virtual ~SampleSetMapGenerator(){}
	
		virtual int StartGeneration();
		virtual std::vector<int> & GenerateNextMap(bool * valid);
		virtual void ScoreLastMap(double score);

		// threshold minValue / maxValue > ratioThreshold
		virtual void AddInvariantFunctionThreshold(SurfaceFeature * feature1, SurfaceFeature * feature2,
												   double ratioThreshold);
		
		// threshold minValue / maxValue > ratioThreshold	
		virtual void AddDistancetoToOtherPoints(SurfaceDistance * distance1, double scale1, 
												SurfaceDistance * distance2, double scale2,
												double ratioThreshold);

		// threshold minValue / maxValue > ratioThreshold	
		virtual void AddMinDistanceThreshold(SurfaceDistance * distance1, double scale1, 
											 SurfaceDistance * distance2, double scale2,
											 double ratioThreshold);
	
		virtual void GetProgress(int & currProg, int & totalNum);
		
	protected:
		SampleSetMapGenerator();
		virtual bool FillTripletAndTest(int A, int B, int C, int Ap, int Bp, int Cp, 
										std::vector<int> & mapToFill);
	
		virtual bool TestDistanceToOther(std::vector<int> & candidateMap);
		virtual bool TestInvariantFunction(std::vector<int> & candidateMap);
		virtual bool TestMinDistanceThreshold(std::vector<int> & candidateMap);
	
		int m_iNumPointsPerMap;
		int m_iMaxNumMaps;
	
		int m_iCurrentMapSequence;
		std::vector<std::vector<int> > m_MapSequencesAfterThresholding;
		std::vector<int> m_Empty;
	
		SurfaceSampleSet * m_pSet1;
		SurfaceSampleSet * m_pSet2;
	
		std::vector<SurfaceFeature*> m_InvariantFunction1;
		std::vector<SurfaceFeature*> m_InvariantFunction2;	
		std::vector<double> m_InvariantRatioThreshold;

		std::vector<SurfaceDistance*> m_DistanceToOtherPoints1;
		std::vector<SurfaceDistance*> m_DistanceToOtherPoints2;	
		std::vector<double> m_DistanceToOtherPointsScale1;
		std::vector<double> m_DistanceToOtherPointsScale2;	
		std::vector<double> m_DistanceToOtherRatioThreshold;
	
		std::vector<SurfaceDistance * > m_MinDistances1;
		std::vector<SurfaceDistance * > m_MinDistances2;	
		std::vector<double> m_MinDistanceScale1;
		std::vector<double> m_MinDistanceScale2;	
		std::vector<double> m_MinDistanceThresholds;
	
		GenerationMethod m_GenSearchMethod;
		double m_fDistanceNormalization1;
		double m_fDistanceNormalization2;
		bool m_bEnforceSameSampleType;
};

class MapGeneratorSymmetric : public SampleSetMapGenerator
{
	public:
		MapGeneratorSymmetric(SurfaceSampleSet * set1, int numPntsPerMap, int maxNumMaps,
							  GenerationMethod method, double surfaceArea);
		virtual ~MapGeneratorSymmetric(){}
		// threshold minValue / maxValue > ratioThreshold
		virtual void AddInvariantFunctionThreshold(SurfaceFeature * feature, double ratioThreshold);
		// threshold minValue / maxValue > ratioThreshold	
		virtual void AddDistancetoToOtherPoints(SurfaceDistance * distance, double ratioThreshold);
		// threshold distance > minDistance
		virtual void AddMinDistanceThreshold(SurfaceDistance * distance, double minDistance);
	
		virtual int StartGeneration();
	
	protected:
		VKStringList m_SymmetrySubsets;
	
		virtual bool FillTripletAndTest(int stationary, int corrA1, int corrA2, 
										std::vector<int> & mapToFill);
		virtual bool FillQuadrupletAndTest(int corrA1, int corrA2, int corrB1, int corrB2,
										   std::vector<int> & mapToFill);	
};

class MapGeneratorAdaptive : public SampleSetMapGenerator
{
	public:
		MapGeneratorAdaptive(SampleSetMapGenerator * otherMap, int maxNumMaps);
		virtual ~MapGeneratorAdaptive(){}
		virtual int StartGeneration();
		virtual std::vector<int> & GenerateNextMap(bool * valid);
		virtual void ScoreLastMap(double score);
		virtual void GetProgress(int & currProg, int & totalNum);
	
	protected:
		SampleSetMapGenerator * m_pDefaultGenerator;
		bool m_bSimpleExhaustive;
		std::vector<bool> m_ExploredMap;
		int m_iNumMaps;
};


#endif

