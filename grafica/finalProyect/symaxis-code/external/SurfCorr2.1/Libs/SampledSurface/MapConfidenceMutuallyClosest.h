#ifndef __MAP_CONFIDENCE_MUTUALLY_CLOSEST_H
#define __MAP_CONFIDENCE_MUTUALLY_CLOSEST_H

#include "SurfaceMapConfidence.h"

/**
 * Approximately tells how good does the map preserve area. 
 * Based on pairwise assigment
 */
class MapConfidenceMutuallyClosest : public SurfaceMapConfidence
{
	public:
		enum ValueType
		{
			CMV_FRAC_MCN,
			CMV_INV_DIST
		};	
	
		MapConfidenceMutuallyClosest(SurfaceMap * corrMap, const VKString & set1, 
									 const VKString & set2, ValueType valueType, 
									 const VKString & dist1="default", 
									 const VKString & dist2="default");	
		virtual ~MapConfidenceMutuallyClosest();
	
		virtual VKString GetMapConfidenceType();	
		virtual double Confidence();
	
		static double GetFractionOfMutuallyClosest(SurfaceMidEdgeConf * mc1, 
												   SurfaceMidEdgeConf * mc2,
												   SurfaceSampleSet * sampleSet1, 
												   SurfaceSampleSet * sampleSet2,
												   LinAlgComplex * origCoords1,
												   LinAlgComplex * origCoords2,
												   LinAlgComplex * transCoords1, 
												   LinAlgComplex * transCoords2,												   
												   MapConformal * confMap);
	
		static void ConstructCoordsHolderForFastMCN(SurfaceMidEdgeConf * mc1, 
													SurfaceMidEdgeConf * mc2,
													SurfaceSampleSet * sampleSet1, 
													SurfaceSampleSet * sampleSet2,
													LinAlgComplex ** origCoords1, 
													LinAlgComplex ** origCoords2,
													LinAlgComplex ** transCoords1, 
													LinAlgComplex ** transCoords2);
	
	protected:	
		virtual SurfaceMapConfidence * Clone(SurfaceMap * currMap);
	
		virtual double CalculateConfidenceAtSample(const SurfaceSample & sample);
		virtual double CalculateErrorAtSample(const SurfaceSample & sample);
	
		double m_fEpsilon;
		ValueType m_ValueType;
		VKString m_Distance1Name;
		VKString m_Distance2Name;
		
		MapCoarse * m_pCorrespondingMCNMap;
			
	protected:
};

#endif


