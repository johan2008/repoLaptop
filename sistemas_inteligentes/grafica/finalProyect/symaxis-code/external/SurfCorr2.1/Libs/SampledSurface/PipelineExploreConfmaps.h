#ifndef __PipelineExploreConfMaps_H
#define __PipelineExploreConfMaps_H

#include "AnalysisPipeline.h"
#include "LinAlgMatrixSparseReal.h"
#include "Sortable.h"

using namespace std;

class OneStepSpectralConsistentConfMaps
{
	public:
		bool AddMap(MapConformal* confMap, double eigenvectorValue, int id);
		MapMultiConformal * CreateMultiConformalMap(std::map<VKString, bool> & existingCorrs,
													LinAlgMatrixReal & cachedConfVals,
													SurfaceDistance * dist1, double distNorm,
													SurfaceSampleSet * integrationSet,
													const VKString & multiConfWeights,
													const VKString & interpStr);
		std::map<int, int> m_Corrs12;
		std::map<int, int> m_Corrs21;		
		std::vector<MapConformal*> m_ConfMaps;
		std::vector<int> m_OriginalIDs;
		std::vector<double> m_EigenvectorValue;
	
};

struct SHREC10CorrVertex
{
	R3Point pnt;
	int id;
	int triID;
	double bary[3];
};
R3Point global_SHREC10CorrVertex(SHREC10CorrVertex * vertex, void * data);

class PipelineExploreConfmaps : public AnalysisPipeline
{
	public:	
		PipelineExploreConfmaps(ParamParser & params);
		virtual ~PipelineExploreConfmaps();

		virtual void OneStepSpectralCreateMultiConformalFromIBM(const VKString & algorithmName);
		virtual void OneStepSpectralLoadConformalMaps(const VKString & algorithmName);
		virtual void OneStepSpectralGenerateConformalMaps(const VKString & algorithmName);
		virtual void OneStepSpectralFillIntegratedBlendedMatrix(const VKString & algorithmName);
		virtual void OneStepSpectralFindEigenvalues(const VKString & algorithmName);
		virtual void LoadAndProcessOneStepSpectralMap(const VKString & algorithmName);
	
		virtual void SurfaceCorrespondenceConfMap(const VKString & algorithmName);
		virtual void LoadAndProcessMaps(const VKString & algorithmName);
			// Explore conformal maps, aggregate maps
		virtual void FindCoarseCorrespondences(const VKString & algorithmName);
		virtual void LoadExternalMap(const VKString & algorithmName, 
									 bool & stillNeedCoarse, bool & stillNeedCoarseToFine);
			// use our method, GMDS, or 'none' if input map is already full
		virtual void InterpolateCoarseCorrespondences(const VKString & algorithmName);
		virtual void InterpolateCoarseMultiConfSmartTakeClosest(const VKString & algorithmName);
		virtual int AddLoadedCoarse(const VKString & algorithmName, 
									const VKString & extrapolationType);
	
		virtual void LoadCoarseCorrsFromTruth();		
		virtual void ExploreConformalMaps();
		virtual bool PickNextNplet(int * curr=NULL, int * total=NULL);
		virtual void VoteForCurrentConfMap(SurfaceMidEdgeConf * M1, 
										   SurfaceMidEdgeConf * M2,
										   const VKString & confidenceName,
										   bool mobiusVoting);
		virtual VKString GetCurrentVoteDescriptor();

		virtual void WritePostFinalProcessingResults();
		virtual void WriteFuzzymaps();
		virtual void WriteDenseSHREC10(const VKString & benchMeshFromNull,
									   const VKString & benchMeshToXform,
									   const VKString & outputCorrFile);
		virtual void WriteSparseSHREC10(const VKString & benchMeshFromNull,
										const VKString & benchMeshToXform,
										const VKString & outputCorrFile,
										int numSparseCorrs,
										const VKString & samplesetName,
										const VKString & confidenceName);
		virtual void ReadSHREC10CorrVerticesIntoKDTree(const VKString & benchMeshSHREC10Corr, 
													   int * numVertices,
													   SHREC10CorrVertex ** vertexArray,
													   R3Kdtree<SHREC10CorrVertex*> ** kdtree=NULL);
	
	// Rendering
		virtual void Draw(ParamParser & drawParams);
	//	virtual void DrawClusterArrow(double minX, double maxX, double minY, double maxY, int elementID, bool from);
		virtual void DrawGeneratorSet(SampledSurface * surface, int surfID, char secondChar, ParamParser * params,
									  const VKString & renderingParams="RendererDefault", 
									  const VKString & surfaceName="none");
		virtual void DrawSymmetricGeneratorSet(SampledSurface * surface, ParamParser * params,
											   const VKString & renderingParams="RendererDefault", 
											   const VKString & surfaceName="none");


	
	// additional stuff
		virtual void RescoreCollection(MapScoredCollection * collection, 
									   const VKString & confidenceName);

	// some ad-hoc temp stuff
		virtual void WriteSymmetryFunction(const VKString & outputFunction, 
										   const VKString & relation);
	
		virtual void PrintBlendingMatrixAtCurrentPoint(ParamParser * drawParams);
		virtual void PrintIntegratedBlendingMatrix(ParamParser * drawParams);
		virtual void MapSelectedBasedOnIntegratedBlendedMatrix(ParamParser * drawParams);
		
	protected:
		MapScoredCollection * m_pCollectionOfConformalMaps;
		std::vector<int> * m_CurrentSequence;
		MapCoarse * m_pFinalCoarseMap;
		VKString m_sCoarseSetName;
	
	// Stuff done right before siggraph
		std::vector<SurfaceSample> m_MapsDueToTopEigenvectors;
		std::vector<double> m_EigenValues;
	
	// eigenanalysis of the integrated blending matrix	
	
		LinAlgMatrixSparseReal * m_pIntegratedBlendingMatrix;
		
		LinAlgMatrixReal * m_pEigenvectorsIBM;
		LinAlgVectorReal * m_pEigenvaluesIBM;	
		LinAlgMatrixReal * m_pCachedConfidenceVals;
		std::vector<Sortable> m_MapIDToEigenvalue;
	
	// optional: fuzzy correspondences
		double FuzzyConsistentConfidenceAtSample(SurfaceMap * surfMap, 
												 const SurfaceSample & samp, 
												 SurfaceMapConfidence * confidence);
		double FuzzyConsistentConfidence(SurfaceMap * surfMap,
										 const VKString & algName,
										 const VKString & confidenceName);
	
		virtual void ReadFuzzyCorrespondences(const VKString & filename);
		
		virtual double GetFuzzyValue(int vertex1, int vertex2);
		int * m_VertexToFuzzySetID1;
		int * m_VertexToFuzzyVertexID2;
		std::map<int, double> * m_aFuzzyValues;
};

#endif




