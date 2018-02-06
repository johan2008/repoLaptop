#ifndef __MAP_COARSE_H
#define __MAP_COARSE_H

#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "gaps.h"

#include "VKString.h"
#include "LinAlgMatrixReal.h"

#include "SampledSurface.h"
#include "SurfaceMap.h"
#include "LinAlgComplex.h"

using namespace std;

class MapGeoFeature;
class SurfaceDistance;
class MapConformal;
class FlatSearchNode;

class MapCoarse : public SurfaceMap
{
	public:	
		static SurfaceSampleSet * GetSetFeaturePoints(SampledSurface * M, 
													  const VKString & verticesToFeatPnts);
		
		enum FineToCoarseGeneration
		{
			F2C_MUTUALLY_CLOSEST_NEIGHBORS,
			F2C_LINEAR_ASSIGNMENT,
			F2C_FORWARD_ONLY,
			F2C_BACKWARD_ONLY,
			F2C_FORWARD_BACKWARD,
			F2C_MCN_CONFORMAL_EUCLIDEAN,		// works only if SurfaceMap is conformal
			F2C_MCN_ONLY_ON_SURF2_FORWARD,
			F2C_MCN_ONLY_ON_SURF1_BACKWARD
		};
	
		MapCoarse(SampledSurface * M1, SampledSurface * M2);	
		MapCoarse(SampledSurface * M1, SampledSurface * M2, 
				  const VKString & sampleSet1, const VKString & sampleSet2);
		MapCoarse(SampledSurface * M1, SampledSurface * M2, 
				  SurfaceSampleSet * sampleSet1, SurfaceSampleSet * sampleSet2);
		MapCoarse(SurfaceSampleSet * sampleSet1, SurfaceSampleSet * sampleSet2,
				  SurfaceMap * surfaceMap, FineToCoarseGeneration method,
				  const VKString & distanceMetric1="default", const VKString & distanceMetric2="default");
	
		void AssociateCoarseToRender(MapCoarse * mapCoarse);
		void PreinitializeAllToNULL();
		void InitializeSets(SurfaceSampleSet * sampleSet1, SurfaceSampleSet * sampleSet2);
		void InitializeVoting();

		virtual ~MapCoarse();

		virtual void ClearAllCorrespondences();

		virtual void SetFinalCorrMap(const vector<int> & M1ToM2);
		
		virtual double CastVote(SurfaceMap * lowDimMap, FineToCoarseGeneration voteForMCN,
							  const VKString & confidenceName="none");
	
		virtual void CastVote(double value, int s1, int s2);
		virtual void CastVote(double value, const SurfaceSample & s1, const SurfaceSample & s2);
	
		virtual void AssignGreedyCorr(double takeBestMaps, bool furtherPairsSearch);
	// further pairs search support functions
		virtual bool FPS_NearestFeatureWithinThrehsold(SurfaceSampleSet * searchInSet,
													   double * refFeature,
													   const SurfaceSample & refSample,
													   double ** featuresInSearch,
													   SurfaceDistance * distMetric,
													   double threshold, int numGoodCorrs);
		virtual double FPS_FeatureDistance(double * a, double * b, int vecSize);
	
		virtual void CleanupMap(double takeBestMaps);	
		virtual void FillBestCorr(SurfaceSampleSet ** saveSet1, SurfaceSampleSet ** saveSet2,
								  int numBest) const;
//		virtual MapGeoFeature * GetInterpolatedGeoFeature(double fracBest, const VKString & distMetric1, 
//														  const VKString & distMetric2);
	
	// inherited map functions
		virtual SurfaceSample ForwardMap(const SurfaceSample & s);
		virtual SurfaceSample InverseMap(const SurfaceSample & s);
		virtual const SurfaceSampleSet * GetValidDomain() const;
		virtual const SurfaceSampleSet * GetValidRange() const;
		virtual int GetMappedDomainSize();
		virtual VKString GetSurfaceMapType();

	// Read / Save map in some specific format
		// Funkhouser's *.map file
		virtual bool LoadMapInOneWayVertexIDsFormat(const VKString & filename, bool surf0ToSurf1);
		// Funkhouser's *.cor file
		virtual bool LoadMapInVertexByVertexFormat(const VKString & filename);
		virtual void SaveMapInVertexByVertexFormat(const VKString & filename);
		// dense preceise: vertex to barycentric coordinates
		virtual bool LoadMapDensePreceiseFormat(const VKString & filename);

		virtual void SaveMap(std::ofstream & textStream);
		virtual void LoadMap(std::ifstream & textStream);

		virtual void WriteCorrByRow(std::ofstream & textStream);
		virtual bool LoadCorrByRow(std::ifstream & textStream);
	
		virtual bool LoadCorrFeaturePntSymmetryMap(const VKString & verticesToFeatPnts,
												   const VKString & featurePntsMap);
		virtual bool LoadCorrSameFeaturePntIDs(const VKString & verticesToFeatPnts1, 
											   const VKString & verticesToFeatPnts2);
		virtual bool LoadCorrSameFeaturePntIDsFlipSymmetry(const VKString & verticesToFeatPnts1, 
														   const VKString & verticesToFeatPnts2,
														   const VKString & symInfoFile);
		virtual bool LoadCorrSameVertexIDs();
		static void FillFeatureToVertex(const VKString & vertexIDsOfFeatPnts, 
										std::vector<int> &vertexIDs);
		static bool FillFeatureToFeatureMap(const VKString & corrInfoFile, 
											std::map<int, int> & featureToFeatureMap);
	
	// Rendering details
		virtual void SetCompareToTruthRendering(const VKString & confidenceName, MapCoarse * trueMap);
		virtual VKString GetTruthConfidenceName();
		virtual void SetCorrespondenceColor(double r, double g, double b);
		virtual void Draw(AnalysisWindow * window,
						  ParamParser * params,
						  const VKString & renderingParams="RendererDefault", 
						  const VKString & surfaceName="none");
	
	// static functions
		static FineToCoarseGeneration StrToCoarseningMethod(const VKString & method);	
	
	protected:
	// creating a coarse map: simple forw/backw (caching), processed=MCN
		virtual void CreateSimpleCoarse(SurfaceSampleSet * sampleSetFrom, 
										SurfaceSampleSet * sampleSetTo,
										std::vector<int> &correspondenceMap, 
										const SurfaceSampleSet * sampleSet1, 
										const SurfaceSampleSet * sampleSet2,
										SurfaceMap * surfaceMap, FineToCoarseGeneration method);
		virtual void CreateProcessedCoarse(SurfaceSampleSet * sampleSetFrom, 
										   SurfaceSampleSet * sampleSetTo,
										   std::vector<int> &correspondenceMap, 
										   std::vector<double> &dists, 
										   const SurfaceSampleSet * sampleSet1, 
										   const SurfaceSampleSet * sampleSet2,
										   SurfaceMap * surfaceMap, FineToCoarseGeneration method,
										   SurfaceDistance * dist1, SurfaceDistance * dist2);

		void ConformalEuclideanMCN(const SurfaceSampleSet * set1, const SurfaceSampleSet * set2,
								   MapConformal * map1, 
								   std::vector<int> & sampleMap, std::vector<double> & corrDist);
	
		void AssignMutuallyClosest(const SurfaceSampleSet * set1, const SurfaceSampleSet * set2, 
								   LinAlgMatrixReal & distMatrix1To2, LinAlgMatrixReal & distMatrix2To1,
								   std::vector<int> & sampleMap, std::vector<double> & corrDist);
			
		void FillDistanceMatrix(const SurfaceSampleSet * set1, const SurfaceSampleSet * set2, 
								LinAlgMatrixReal & distMatrix1To2, LinAlgMatrixReal & distMatrix2To1,
								SurfaceDistance * dist1, SurfaceDistance * dist2, SurfaceMap * surfaceMap, 
								FineToCoarseGeneration f2cmethod);
	
		void InitializeOrderOfDomainRendering();
		
		map<int, int> m_mVertexIDToSetID1;
		map<int, int> m_mVertexIDToSetID2;

		vector<int> m_mFinalCorrMapM1ToM2;
		vector<double> m_mFinalCorrValue;
		vector<int> m_OrderOfDomainRendering;

		LinAlgMatrixReal * m_VotedCorrMatrix;

		SurfaceSampleSet * m_pSet1;
		SurfaceSampleSet * m_pSet2;
	
		MapCoarse * m_pTruthMapForDrawing;
		double m_CorrColorR, m_CorrColorG, m_CorrColorB;
		VKString m_sTruthConfidence;
		double m_fEpsilon;
	
		MapCoarse * m_pRenderAnotherCoarse;
	
	// temp stuff for mobius voting. Cleanup in future
		LinAlgComplex * m_pOriginalComplexCoords1;
		LinAlgComplex * m_pOriginalComplexCoords2;
		LinAlgComplex * m_pTransformedComplexCoords1;
		LinAlgComplex * m_pTransformedComplexCoords2;
	
		FlatSearchNode * m_pSearchNodeArray1;
		FlatSearchNode * m_pSearchNodeArray2;	
		int m_iSearchArraySize1;
		int m_iSearchArraySize2;	
};

#endif




