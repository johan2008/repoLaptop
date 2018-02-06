#
# Package name and list of source files.
#

NAME=SampledSurface
CCSRCS=	AnalysisWindow.cpp				\
	AnalysisPipeline.cpp				\
		PipelineExploreConfmaps.cpp		\
		PipelineGeneral.cpp			\
		PipelineExploreConfmapsOneStepSpectral.cpp	\
	AnalysisStats.cpp				\
	PlanarTransform.cpp				\
		MobiusTransformation.cpp 		\
		PlanarTransformLSCM.cpp			\
		PlanarTransformQuasiConformal.cpp	\
		PlanarTransformMVC.cpp			\
	SurfaceSample.cpp 				\
	SurfaceSampleSet.cpp 				\
	SampleSetMapGenerator.cpp			\
	SurfacePerVertexValue.cpp			\
	SurfaceFeature.cpp 				\
		FeatureCurvature.cpp			\
		FeatureAGD.cpp 				\
		FeatureMGD.cpp  			\
		FeatureIRSA.cpp				\
	SampledSurface.cpp 				\
		SurfaceMidEdgeConf.cpp			\
	SurfaceTexture.cpp				\
	SurfaceDistance.cpp 				\
		DistanceGeodesic.cpp 			\
		DistanceLazySubsets.cpp 		\
		DistanceOnTheFly.cpp			\
	SurfaceMap.cpp 					\
		MapFlattening.cpp			\
		MapMultiConformal.cpp			\
		MapScoredCollection.cpp			\
		MapCoarse.cpp 				\
		MapConformal.cpp			\
		MapEuclidean.cpp 			\
		MapGeoFeature.cpp 			\
	SurfaceMapConfidence.cpp			\
		MapConfidenceMutuallyClosest.cpp	\
		MapConfidenceDiscreteArea.cpp		\
		MapConfidenceFaceArea.cpp		\
		MapConfidenceDistanceToGen.cpp		\
		MapConfidenceCompareToTruth.cpp		\
		MapConfidenceMultiConf.cpp		\
	SurfaceMapSimilarity.cpp			\
		MapSimilarityOverlap.cpp		\
		MapSimilarityDistance.cpp		\


PKG_LIBS=BasicRoutines BasicNumerical
EXT_LIBS = gsl CXSparse UFconfig

#
# Set up compiler options, etc.
#
include ../../scripts/Makefile.pkgs






