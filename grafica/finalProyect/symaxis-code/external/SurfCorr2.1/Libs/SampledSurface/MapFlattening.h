#ifndef __MAP_FLATTENING_H
#define __MAP_FLATTENING_H

#include "SurfaceMap.h"
#include "PlanarTransform.h"
#include "SurfaceMidEdgeConf.h"

/**
 * This class maps 3D surface to 2D plane (unbounded)
 */
class MapFlattening : public SurfaceMap
{
	public:		
		MapFlattening(SampledSurface * s1);
		MapFlattening(SampledSurface * flattenMe, 
					  PlanarTransform * transform,
					  bool bilateralReflection = false,
					  R3Mesh * flatMesh = NULL,
					  R2Kdtree<FlatSearchNode*> * flatMeshTree = NULL);
	
		virtual ~MapFlattening();
	
		virtual void SetPlanarTransformation(const VKString & coarseSampleSet,
											 std::vector<int> & sampleIDs,
											 std::vector<LinAlgComplex> & samples2DPositions,
											 const VKString & planarTransformClass="SameAsCurrent");
	
	
		virtual void InteractiveUpdateFromCorrs(const VKString & sampleSet,
												const VKString & textureName);
	
		virtual void DrawFlatMesh(ParamParser & drawParams);

	
		/** Maps 3D surface to 2D plane. Note that returned sample should have Z()=0 */
		virtual SurfaceSample ForwardMap(const SurfaceSample & s);
		/** Maps 2D point on a plane to 3  */
		virtual SurfaceSample InverseMap(const SurfaceSample & s);
	
		/** Gets flattened position of a surface (does not apply planar transform) */
		virtual LinAlgComplex Get2DPositionFromSample(const SurfaceSample & s);
		/** Gets point on a surface given its flattening (does not apply planar transform) */
		virtual SurfaceSample GetSampleFrom2DPosition(const LinAlgComplex & s);
	
		/** Maps sample on 3D surface (m_pM1) to m_pFlatMesh */
		virtual SurfaceSample MapSurfaceToFlatMesh(const SurfaceSample & s) = 0;
		/** Maps sample on m_pFlatMesh to 3D surface (m_pM1)  */
		virtual SurfaceSample MapFlatMeshToSurface(const SurfaceSample & s) = 0;
	
		void SetFlatMesh(R3Mesh * flatMesh, R2Kdtree<FlatSearchNode*> * kdTree=NULL);
	
		double GetFlatMeshArea();
		R3Mesh * GetFlatMesh();
		R2Kdtree<FlatSearchNode*> * GetKdTree();
	
	protected:
		bool m_bBilateralReflection;
		double m_fFlatMeshArea;
		FlatSearchNode * m_aSearchNodes;	
		R3Mesh * m_pFlatMesh;
		PlanarTransform * m_pPlanarTransform;
		
		R2Kdtree<FlatSearchNode*> * m_pFlatMeshKdVertexSearch;
};

/**
 * Mid-Edge flattening. Note that this class heavily uses SurfaceMidEdgeConf - 
 * with this class one can use simple SampledSurface).
 */ 
class MapFlatMidEdge : public MapFlattening
{
	public:	
		MapFlatMidEdge(SampledSurface * s1);
		MapFlatMidEdge(SampledSurface * flattenMe, 
					   PlanarTransform * transform,
					   bool bilateralReflection = false,
					   R3Mesh * flatMesh = NULL, 
					   R2Kdtree<FlatSearchNode*> * flatMeshTree = NULL);
	
		MapFlatMidEdge(SampledSurface * flattenMe, 
					   const VKString & findMobisuSampleSet, 
					   std::vector<int> & sampleIDs,
					   std::vector<LinAlgComplex> & samples2DPositions,
					   const VKString & planarTransformClass,
					   bool bilateralReflection = false,
					   R3Mesh * flatMesh = NULL, 
					   R2Kdtree<FlatSearchNode*> * flatMeshTree = NULL);
	
		MapFlatMidEdge(SampledSurface * flattenMe, 
					   PlanarTransform * transform, 
					   SurfaceFeature * agd,
					   bool bilateralReflection = false,
					   R3Mesh * flatMesh = NULL, 
					   R2Kdtree<FlatSearchNode*> * flatMeshTree = NULL);
		
		MapFlatMidEdge(SampledSurface * flattenMe, 
					   const VKString & findMobisuSampleSet, 
					   std::vector<int> & sampleIDs,
					   std::vector<LinAlgComplex> & samples2DPositions,
					   const VKString & planarTransformClass,
					   SurfaceFeature * agd,
					   bool bilateralReflection = false,
					   R3Mesh * flatMesh = NULL, 
					   R2Kdtree<FlatSearchNode*> * flatMeshTree = NULL);
	
		virtual ~MapFlatMidEdge();
	
		virtual VKString GetSurfaceMapType();
		virtual void SaveMap(std::ofstream & textStream);
		virtual void LoadMap(std::ifstream & textStream);

	
		/** Fills m_pFlatMesh and m_pFlatMeshKdVertexSearch using mid-edge flattening */
		virtual void GenerateFlatMesh(int cutTriangle=-1);
		/** Surface to mid-edge mesh mapping */
		virtual SurfaceSample MapSurfaceToFlatMesh(const SurfaceSample & s);
		/** Mid-edge mesh to surface mapping */
		virtual SurfaceSample MapFlatMeshToSurface(const SurfaceSample & s);
	
	protected:
};

/**
 * Inter-surface map via 2D plane (note: this class provides another potential 
 * implementation of MapConformal, architecture of this class is more versatile). 
 */
class MapVia2DPlane : public SurfaceMap
{
	public:
		MapVia2DPlane(SampledSurface * s1, SampledSurface * s2);
		MapVia2DPlane(SampledSurface * surface1, SampledSurface * surface2, 
					  MapFlattening * flattening1, MapFlattening * flattening2);
	
		MapVia2DPlane(SampledSurface * surface1, SampledSurface * surface2, 
					  const VKString & corrSetName, std::vector<int> & correspondences,
					  const VKString & createFlatteningType, 
					  const VKString & planarTransformClass,
					  bool bilateralReflection = false,
					  R3Mesh * flat1=NULL, R2Kdtree<FlatSearchNode*> * tree1=NULL,
					  R3Mesh * flat2=NULL, R2Kdtree<FlatSearchNode*> * tree2=NULL);
	
		virtual ~MapVia2DPlane();
	
		virtual void Draw(AnalysisWindow * window, ParamParser * params,
						  const VKString & renderingParams="RendererDefault", 
						  const VKString & surfaceName="none");
	
		virtual void DrawGenerators(bool drawLines, bool drawSpheres,
									int lineWidth, int sphereScale,
									AnalysisWindow * window, ParamParser * params,
									const VKString & renderingParams="RendererDefault", 
									const VKString & surfaceName="none");
	
		virtual VKString GetSurfaceMapType();
		virtual void SetFlattenings(R3Mesh * flat1, R2Kdtree<FlatSearchNode*> * tree1,
									R3Mesh * flat2, R2Kdtree<FlatSearchNode*> * tree2);
		
		virtual SurfaceSample ForwardMap(const SurfaceSample & s);
		virtual SurfaceSample InverseMap(const SurfaceSample & s);
		
		virtual void SaveMap(std::ofstream & textStream);
		virtual void LoadMap(std::ifstream & textStream);
		
		virtual MapFlattening * GetFlattening(int surfaceID);
	
		virtual void GetGeneratingSet(std::vector<int> ** corrs, SurfaceSampleSet ** sampleSet1, 
									  SurfaceSampleSet ** sampleSet2);
		
	protected:
		MapFlattening * m_pFlattening1;
		MapFlattening * m_pFlattening2;
		VKString m_CorrsSampleSet;
		std::vector<int> m_Corrs;
};



#endif

