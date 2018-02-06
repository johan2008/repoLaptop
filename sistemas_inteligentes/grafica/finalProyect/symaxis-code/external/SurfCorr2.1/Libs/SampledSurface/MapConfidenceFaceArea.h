#ifndef __MAP_CONFIDENCE_FACE_AREA_H
#define __MAP_CONFIDENCE_FACE_AREA_H

#include "SurfaceMapConfidence.h"

/**
 * Compares areas of a triangle before and after maping (suppose it maps to a plane)
 */
class MapConfidenceFaceArea : public SurfaceMapConfidence
{
	public:
		enum ValueType
		{
			TRIAREA_VAL_DIFF,
			TRIAREA_VAL_MIN_OVER_MAX,
			TRISTRETCH_SV_RATIO,
			TRISTRETCH_SV_L2STRETCH,
			TRISTRETCH_SV_AREA
		};
		
		MapConfidenceFaceArea(SurfaceMap * correspondingMap, ValueType valType, 
							  const VKString & samplesFrom, const VKString & samplesTo);
		virtual ~MapConfidenceFaceArea();
		virtual VKString GetMapConfidenceType();
	
//	protected:
		virtual SurfaceMapConfidence * Clone(SurfaceMap * currMap);
	
		virtual double CalculateConfidenceAtSample(const SurfaceSample & sample);
	
		static void FillSamples(bool forward, int triangleID, SurfaceMap * map1,
								double * A1save, double * A2save,
								SurfaceSample & s11, SurfaceSample & s12, SurfaceSample & s13,
								SurfaceSample & s21, SurfaceSample & s22, SurfaceSample & s23,
								double * Anorm1, double * Anorm2);
		static double GetPerTriangleErrorStretch(ValueType tp, bool forward, int triangleID, 
												 SurfaceMap * map1, double * A1=NULL, double * A2=NULL);
		static double GetPerTriangleErrorArea(ValueType tp, bool forward, int triangleID, 
											  SurfaceMap * map1, double * A1=NULL, double * A2=NULL);
		static double GetPerVertexConfidence(ValueType tp, bool forward, int vertexID, SurfaceMap * map1);
		static double GetPerVertexError(ValueType tp, bool forward, int vertexID, SurfaceMap * map1);
	
		ValueType m_ValueType;

};

#endif


