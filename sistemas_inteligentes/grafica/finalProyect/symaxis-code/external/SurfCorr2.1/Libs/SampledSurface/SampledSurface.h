#ifndef __SAMPLED_SURFACE_H
#define __SAMPLED_SURFACE_H

#include "SurfaceSampleSet.h"
#include "SurfaceSample.h"
#include "ParamParser.h"
#include "AnalysisWindow.h"


class SurfaceFeature;
class SurfaceDistance;
class SurfaceMap;
class SurfaceMapSimilarity;
class SurfaceMapValue;
class DistanceOnTheFly;

/**
 * This is surface class: surface is defined by a 3D mesh, and can optionally
 * contain SuraceDistance(s), SurfaceSampleSet(s), SurfaceFeature(s), each uniquely identified
 * by a name. Note that name "default" can also be mapped to a set, a feature, or a 
 * distance metric for rendering and for most algoirthms.
 */ 
class SampledSurface
{
	public:
		enum GeodesicCentroidAlgorithmType
		{
			GEO_CENTROID_EVERY_VERTEX,
			GEO_CENTROID_NEIGHBORHOOD,
			GEO_CENTROID_BEST_AMONG_GIVEN,
			GEO_CENTROID_EUCLIDEAN_APPROX
		};
	
		SampledSurface(R3Mesh * mesh);
		virtual ~SampledSurface(){}

		/** 
		* Copy Sample Set, Features and Distances (e.g. if surfaces are for symmtery detection, they
		*	must be same, thus features, distances and sample sets are same, however, conformal
		*	flatening (described in surface), are different).
		*/
		virtual void CopyFromSameSurface(SampledSurface * sameSurface);

		virtual void AddMapActingAsFrom(const VKString & mapName, SurfaceMap * surfMap );
		virtual SurfaceMap * GetMapActingAsFrom(const VKString & mapName);
		virtual void AddMapActingAsTo(const VKString & mapName, SurfaceMap * surfMap );
		virtual SurfaceMap * GetMapActingAsTo(const VKString & mapName);
	
		virtual void AddSampleSet(const VKString & setName, SurfaceSampleSet * sampSet );
		virtual SurfaceSampleSet * GetSampleSet(const VKString & setName);

		virtual void AddFeature(const VKString & featureName, SurfaceFeature * feature);
		virtual SurfaceFeature * GetFeature(const VKString & featureName);

		virtual void AddDistanceMetric(const VKString & metricName, SurfaceDistance * distances);
		virtual SurfaceDistance * GetDistanceMetric(const VKString & metricName);
		virtual DistanceOnTheFly * GetOnTheFlyDistanceMetric(double maxDistanceUnscaled, 
															 const VKString & tryMeFirstMetric, int maxCache);
	
		virtual VKStringList GetAllFeatureNames();
		virtual VKStringList GetAllSampleSetNames();
		virtual VKStringList GetAllDistanceMetricNames();

		virtual R3Mesh * GetMesh();
	
		virtual SurfaceSample WeightedGeodesicCentroid(std::vector<SurfaceSample> & samples,
													   std::vector<double> & weights,
													   GeodesicCentroidAlgorithmType type = GEO_CENTROID_NEIGHBORHOOD, 
													   double outlierThreshold = 0);	
		virtual double GetWeightedGeoCentroidError(const SurfaceSample & s, 
												   std::vector<SurfaceSample> & samples, 
												   std::vector<double> & weights, 
												   SurfaceDistance * distanceMetric);
		virtual SurfaceSample WeightedGeoCentrEveryVertex(std::vector<SurfaceSample> & samples,
														  std::vector<double> & weights);
		virtual SurfaceSample WeightedGeoCentrNhdOnly(std::vector<SurfaceSample> & samples,
													  std::vector<double> & weights, double thresh);
		virtual SurfaceSample WeightedGeoCentrBestAmongGiven(std::vector<SurfaceSample> & samples,
															 std::vector<double> & weights);
		virtual SurfaceSample WeightedGeoCentrEuclideanApprox(std::vector<SurfaceSample> & samples,
															  std::vector<double> & weights, 
															  double outlierThreshold);
		virtual SurfaceSample GetNearestSample(const R3Point & queryPnt, 
											   double minDist = 0, double maxDist = FLT_MAX);
	
		virtual void SetParentWindow(AnalysisWindow * window);
	
		virtual void SelectFeature(const VKString & featureName);
		virtual void SelectDistanceMetric(const VKString & distanceName);
		virtual void SelectSampleSet(const VKString & sampleName);	
	
		virtual void LoadCamera();
	
		virtual void ClearCachedSurfaceColors();
	
		virtual void SetDrawSurfaceColors(ParamParser * params, 
										  const VKString & mainSet="none", SurfaceMap * currMap = NULL,
										  SurfaceMap * otherMap=NULL, const VKString & similarityName="none",
										  const VKString & confidenceName="none",
										  const VKString & renderingParams="RendererDefault",
										  const VKString & surfaceName="none");
		virtual void Draw(ParamParser * drawParams,
						  const VKString & renderingParams="RendererDefault", 
						  const VKString & surfaceName="none");

		virtual void DrawSurface(ParamParser * drawParams,
								 const VKString & renderingParams="RendererDefault",
								 const VKString & surfaceName="none");
		
		virtual SurfaceSampleSet * GetRenderableSampleSet(ParamParser * drawParams,
														  const VKString & renderingParams="RendererDefault",
														  const VKString & surfaceName="none");
		virtual void DrawSamples(ParamParser * drawParams,
								 const VKString & renderingParams="RendererDefault",
								 const VKString & surfaceName="none");
	
		virtual R3Point GetDrawablePosition(const SurfaceSample & sample, 
											ParamParser * drawParams,
											const VKString & renderingParams="RendererDefault",
											const VKString & surfaceName="none");
	
		virtual double AdjustedRadius(double fractionOfArea);
		virtual double Area();
		virtual double GetStandardRadius(ParamParser * params, 
										 const VKString & renderingParams="RendererDefault", 
										 const VKString & surfaceName="none");
	
		virtual VKString GetSurfaceType();
	
		virtual VKString GetName();
		virtual void SetName(const VKString & name);
	
		static R3Mesh * CreateCopy(R3Mesh * mesh);
	
		int GetSelectedVertex(ParamParser * drawParams);
	
	protected:
		R3Mesh * CheckMeshForProblems(R3Mesh * mesh);
	
		VKString m_sSurfaceName;
		int m_iCameraID;
	
		R3Mesh * m_pMesh;
		std::map<VKString, SurfaceSampleSet *> m_mNameToSampleSet;
		std::map<VKString, SurfaceFeature *> m_mNameToFeatureSet;
		std::map<VKString, SurfaceDistance *> m_mNameToDistance;
		std::map<VKString, SurfaceMap *> m_mNameToMapActingAsFrom;
		std::map<VKString, SurfaceMap *> m_mNameToMapActingAsTo;
		AnalysisWindow * m_pWindow;
		bool valid;

		int m_iCachedVertexFrom;
		int m_iCachedConfID;
		VKString m_sCachedNameMapOnSurf;		
		VKString m_CachedMainColorSet;	
		SurfaceFeature * m_pCachedSurfaceFeature;
		SurfaceDistance * m_pCachedSurfaceDistance;
		SurfaceMap * m_pCachedSurfaceMap;
		SurfaceMap * m_pCachedSurfaceOtherMap;
		VKString m_sCachedMapSimilarity;
		VKString m_sCachedMapConfidence;

		std::vector<double> m_CachedPerVertexColors;
	
		double m_fCachedMinX, m_fCachedMinY, m_fCachedMinZ, m_fCachedMaxX, m_fCachedMaxY, m_fCachedMaxZ;
		bool m_bCachedBBox;

		double m_fArea;
		R3MeshSearchTree * m_pMeshSearchTree;
};

/**
 * Dummy surface: 2D plane
 */
class Surface2DPlane : public SampledSurface
{
	public:
		static SampledSurface * m_pInfinitePlanePseudosurface;
		static R3Mesh * m_pInfinitePlanePseudomesh;
};


#endif




