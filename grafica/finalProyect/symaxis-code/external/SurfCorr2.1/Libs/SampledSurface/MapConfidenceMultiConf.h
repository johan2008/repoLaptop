#ifndef __MAP_CONFIDENCE_MULTI_CONF_H
#define __MAP_CONFIDENCE_MULTI_CONF_H

#include "SurfaceMapConfidence.h"

class MapMultiConformal;

class MapConfidenceMultiConf : public SurfaceMapConfidence
{
	public:
		enum WeightCalculation
		{
			MULTICONF_WEIGHT_CONFIDENCE,
			MULTICONF_WEIGHT_CONFIDENCE_AND_CONSISTENCY,
			MULTICONF_WEIGHT_PER_POINT_EIGENANALYSIS,
			MULTICONF_WEIGHT_PER_POINT_EIGENANALYSIS_CONFIDENCE,
			MULTICONF_WEIGHT_GLOBAL_EIGENANALYSIS			
		};
	
		MapConfidenceMultiConf(SurfaceMap * currMap, 
							   const VKString & setNameFrom,
							   const VKString & setNameTo,
							   const VKStringList & confidences, const VKString & similarity,
							   WeightCalculation weightCalculation);
		virtual ~MapConfidenceMultiConf();
		
		virtual void ClearMultiConfWeights();
		virtual VKString GetMapConfidenceType();
		virtual void FillFinalWeights(std::vector< std::map<int, double> > & weights);
	
		virtual SurfaceMapConfidence * Clone(SurfaceMap * currMap);
	
		virtual void PrepareFinalWeights();
		virtual void SmoothFinalWeights();
		
		virtual double Weight(int mapID, const SurfaceSample & sample);
	
		virtual double GetConfidence(int mapID, int confID, const SurfaceSample & sample);
		virtual double GetTotalConfidence(int mapID, const SurfaceSample & sample);
		virtual double GetConsistency(int mapID, const SurfaceSample & sample);	
		virtual double GetSimilarity(int i, int j, const SurfaceSample & sample);		
	
		virtual void SaveConfidence(std::ofstream & textStream);
		virtual void LoadConfidence(std::ifstream & textStream);
	
		virtual bool WeightsReady();
	
		static WeightCalculation StrToWeightCalculation(const VKString & str);

		static void NormalizeWeights(std::vector<double> & weights);
		static void NormalizeWeights2(std::vector<double> & weights);	
	
	// ad-hoc threshold for extrapolation: if there is no good conformal map for a point,
	// map a point same as the nearest point with acceptable map value.
		void SetIfNoGoodConformalMapThreshold(double threshold);

	//protected:
		double m_fIfNoGoodConformalMapTheshold;
	
		bool m_bWeightsReady;
		bool m_bInterpolatedWeights;	
	
		WeightCalculation m_WeightCalculation;
		bool m_bPresmoothSimConf;
	
		virtual void CalculateWeightsAtSample(std::vector<double> & weights, 
											  const SurfaceSample & sample);
		virtual void CalcWeightsConfidence(std::vector<double> & weights, 
										  const SurfaceSample & sample);
		virtual void CalcWeightsConstant(std::vector<double> & weights);
		virtual void CalcWeightsConfidenceTakeBestConf(std::vector<double> & weights, 
													  const SurfaceSample & sample);
		virtual void CalcWeightsConfidenceConsistency(std::vector<double> & weights, 
													 const SurfaceSample & sample);
		virtual void CalcWeightsLocalEigenAnalysis(std::vector<double> & weights, 
												  const SurfaceSample & sample);
		virtual void FillWeightsUsingGlobalEigenAnalysis(SampledSurface * surface);
	

		virtual double CalculateConfidenceAtSample(const SurfaceSample & sample);
		
		MapMultiConformal * m_pMultiConformalMap;
		bool m_bDeleteMultiConformal;
	
		std::map<SampledSurface *, std::vector<SurfacePerVertexValue * > > m_SuggestedWeights;
		std::map<SampledSurface *, std::map<int, std::vector<double> > > m_mVertexToOptimalWeights;
	
		VKStringList m_Confidences;
		VKString m_Similarity;
		bool m_bCacheWeights;
		
		LinAlgMatrixReal * m_pSimilarityMatrix;
		LinAlgMatrixReal * m_pEigenVectors;
		LinAlgVectorReal * m_pEigenValues;	
};

#endif

