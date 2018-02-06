#ifndef __SURFACE_MID_EDGE_CONF_H
#define __SURFACE_MID_EDGE_CONF_H

#include "SampledSurface.h"
#include "FeatureAGD.h"
#include "LinAlgComplex.h"
#include "LinAlgVectorReal.h"
#include "LinAlgMatrixReal.h"
#include "MobiusTransformation.h"

struct FlatSearchNode
{
	FlatSearchNode(int id, const R2Point & pnt);
	FlatSearchNode();
	int m_iID;
	R2Point m_Pnt;
};

R2Point global_FlatNodePosition(FlatSearchNode * currNode, void * data);

class SurfaceMidEdgeConf : public SampledSurface
{
	public:
		friend class MapFlattening;
		friend class MapFlatMidEdge;
	
		SurfaceMidEdgeConf(R3Mesh * mesh, bool complexFlip);
		SurfaceMidEdgeConf(R3Mesh * mesh, int faceToRip, int canID1, int canID2, int canID3, 
				R3Mesh ** onlyFlattened, bool complexFlip);  //Does flattening	
		SurfaceMidEdgeConf(R3Mesh * mesh, R3Mesh * flatMesh, bool complexFlip);
	
		R3Mesh * Flatten(int numSmooth, SurfaceFeature * agdFeature = NULL);

		void Transform(const SurfaceSample & s1, const SurfaceSample & s2, const SurfaceSample & s3,
					const LinAlgComplex & t1, const LinAlgComplex & t2, const LinAlgComplex & t3);
		void Transform(const MobiusTransformation & m);
		void LoadOriginalFlattening();
		MobiusTransformation GetCurrentTransform();
		LinAlgComplex GetConfCoord(const SurfaceSample & s);
		LinAlgComplex GetConfCoordOriginal(const SurfaceSample & s);
		SurfaceSample GetSurfaceSample(const LinAlgComplex & pos);
		virtual bool Reflected();

		virtual void SetParentWindow(AnalysisWindow * window);	
	
		static void GetBestFlatteningParams(R3Mesh * mesh, int & faceToRip, 
											int & v1, int & v2, int & v3, 
											SurfaceFeature * agd);
	
		virtual double GetStandardRadius(ParamParser * params, 
										 const VKString & renderingParams="RendererDefault", 
										 const VKString & surfaceName="none");
	
		virtual void Draw(ParamParser * params, 
						  const VKString & renderingParams="RendererDefault", 
						  const VKString & surfaceName="none");
	
		virtual void DrawConformal(ParamParser * params, 
								   const VKString & renderingParams="RendererDefault", 
								   const VKString & surfaceName="none");
	
		virtual void LoadConformalCamera();
	
		virtual void DrawConformalDoNotLoadCamera(ParamParser * params, 
												  const VKString & renderingParams="RendererDefault", 
												  const VKString & surfaceName="none");
	
	
		virtual R3Point GetDrawablePosition(const SurfaceSample & sample, 
											ParamParser * params,
											const VKString & renderingParams="RendererDefault",
											const VKString & surfaceName="none");
	
	
		virtual VKString GetSurfaceType();
	
	protected:
		void Initialize(R3Mesh * mesh, R3Mesh * flatMesh, bool complexFlip, bool copy);
		
		bool m_bReflection;
		R3Mesh * m_pFlatMesh;
		MobiusTransformation m_CurrentTransform;
		R3Mesh * m_pFlatMeshTransformed;

		R3MeshSearchTree * m_pFlatMeshTransSearch;
		R2Kdtree<FlatSearchNode*> * m_pFlatMeshKdVertexSearch;
		FlatSearchNode * m_VertexArray;
		bool m_bUpdateTransformedMesh;
		void UpdateFlatSearchTree();

		static void AddVerticesToFlatMesh(R3Mesh* originalMesh, R3Mesh * flattenedMesh);
		static R2Kdtree<FlatSearchNode*> * CreateKdTree(R3Mesh * flattenedMesh, FlatSearchNode ** vertexArray);
	
		static SurfaceSample MapMeshToMidEdge(R3Mesh * mesh, R3Mesh * flat, const SurfaceSample & samp);
		static SurfaceSample MapMidEdgeToMesh(R3Mesh * mesh, R3Mesh * flat, const SurfaceSample & samp);	
		static SurfaceSample MapMidEdgeToMeshNearVertex(R3Mesh * mesh, R3Mesh * flat, 
														const SurfaceSample & samp, 
														int localVertexID, int localToVertexID1, 
														int localToVertexID2);
	
		static int GetTriangleID(int v[3], int p[3], R3Mesh * mesh);
		static int OriginalToMidEdgeMap(int v, R3Mesh * mesh, R3Mesh * flat);
		static int OriginalToMidEdgeMap(int v1, int v2, R3Mesh * mesh, R3Mesh * flat);	
		static int MidEdgeOpposingOriginalVertex(int edge, int otherEdge, R3Mesh * mesh, R3Mesh * flat);
		static int MidEdgeAdjacentOriginalVertex(int flatVertexFromEdgeID, int notThisVertex,
											  R3Mesh * mesh, R3Mesh * flat);
		static void MidEdgeToOriginalMap(int flatVertexID, int & edgeID, int & vertexID, R3Mesh * mesh);

		static LinAlgComplex GetConfCoord(const SurfaceSample & s, R3Mesh * mesh, R3Mesh * flat);

		// FLATTENING - IMPLEMENTED BY TOM FUNKHOUSER
		static void TransformToCanonicalFrame(R3Mesh * mesh, R3Mesh * flat, int v1, int v2, int v3);
		static void Transform(MobiusTransformation m, R3Mesh * meshIn, R3Mesh * meshOut);
		static SurfaceSample GetSurfaceSampleSearchTree(const LinAlgComplex & pos, R3Mesh * flatMesh, 
														R3MeshSearchTree * flatMeshSearch, 
														R3Mesh * mesh, bool reflected);
		static SurfaceSample GetSurfaceSampleKdTree(const LinAlgComplex & pos, R3Mesh * flatMesh, 
													R2Kdtree<FlatSearchNode*> * flatMeshSearch,
													R3Mesh * mesh);
		static SurfaceSample GetFlatSurfaceSampleKdTree(const LinAlgComplex & pos, R3Mesh * flatMesh, 
														R2Kdtree<FlatSearchNode*> * flatMeshSearch);
	
		static double Cotangent(const R3Point& position0, const R3Point& position1, const R3Point& position2);
		static LinAlgVectorReal * ComputeXCoordinates(R3Mesh *mesh, R3MeshFace *remove_face);
		static LinAlgVectorReal * ComputeYCoordinates(R3Mesh *mesh, R3MeshFace *remove_face, LinAlgVectorReal &x);
		static R3Mesh * CreateFlattenedMidEdgeMesh(R3Mesh *mesh, R3MeshFace *remove_face, 
						LinAlgVectorReal *x, LinAlgVectorReal *y);
		static R3Mesh * CreateFlattenedMidEdgeMesh(R3Mesh *mesh, int removeFaceID);
		
		int m_iMaxNumConfCoordsToCache;
		std::map<int, LinAlgComplex> m_iVertexToConfCoord; 
};

#endif
