#include "AnalysisPipeline.h"
#include "AnalysisWindow.h"
#include "GLUTSurfaceWindow.h"
#include "LinAlgMatrixSparseReal.h"
#include "LinAlgMatrixReal.h"
#include "LinAlgVectorReal.h"
#include "svdlib.h"

int main(int argc, char** argv )
{	

	bool valid;
	if (argc<2)
	{
		std::cout<<"Usage: "<<argv[0]<<" renderingSettings.txt [optionalSettings]"<<std::endl;
		exit(0);
	}
	const char ** constArgv=(const char**)argv;
	ParamParser params(argv[1], &(constArgv[2]), argc-2);
	AnalysisPipeline * pipeline = AnalysisPipeline::CreatePipeline(params);
	
	if (params.GetStrValue("Pipeline", "Visualize", valid)=="true" && valid)
	{
		ParamParser drawParams(params.GetStrValue("Pipeline", "VisSettings", valid));
		assert(valid);
		GLUTSurfaceWindow vis(drawParams, pipeline, &argc, argv);
	}
	
	return 0;
}

