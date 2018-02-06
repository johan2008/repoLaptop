#ifndef __GLUT_SURFACE_WINDOW_H
#define __GLUT_SURFACE_WINDOW_H

#include <fstream>
#include <sstream>

#include "gaps.h"

#include "VKString.h"
#include "ParamParser.h"

#include "SampledSurface.h"
#include "SurfaceMap.h"
#include "SurfaceMidEdgeConf.h"
#include "DistanceGeodesic.h"
#include "FeatureAGD.h"
#include "FeatureMGD.h"
#include "MapCoarse.h"
#include "MapConformal.h"
#include "MapEuclidean.h"
#include "MapGeoFeature.h"
#include "SurfaceMap.h"
#include "AnalysisPipeline.h"
#include "AnalysisWindow.h"
#include "PipelineExploreConfmaps.h"
#include "OpenGLHeaders.h"

using namespace std;

void GLUTRedraw();
void GLUTResize(int w, int h);
void GLUTKeyboard(unsigned char key, int x, int y);
void GLUTSpecial(int key, int x, int y);
void GLUTMouse(int button, int state, int x, int y);
void GLUTMotion(int x, int y);

/**
 * GLUT implementation of UI winodw.
 */
class GLUTSurfaceWindow : public AnalysisWindow
{
	public:
		GLUTSurfaceWindow(ParamParser & drawParams, AnalysisPipeline * pipeline,
						  int * argc, char ** argv);

		virtual void RenderText(R3Point pnt, const VKString & text, 
								double r, double g, double b, TextSize size);
		
		virtual void Terminate();
		virtual void Update();
	
		void GLUTRedrawCallback();
		void GLUTResizeCallback(int w, int h);
		void GLUTKeyboardCallback(unsigned char key, int x, int y);
		void GLUTSpecialCallback(int key, int x, int y);
		void GLUTMouseCallback(int button, int state, int x, int y);
		void GLUTMotionCallback(int x, int y);
	
		void KeyboardCallbackCorrespondence(unsigned char key, int x, int y);
		void KeyboardCallbackGeneral(unsigned char key, int x, int y);
		void UpdateFeature();
		void UpdateSampleSet();	
		
		static GLUTSurfaceWindow * s_MainGLUTWindow;

	
	protected:
		int m_GLUTwindowID;
		int m_GLUTmouse[2];
		int m_GLUTbutton[3];
		int m_GLUTmodifiers;

		int m_iFeatureID;
		int m_iSampleSetID;
		int m_iConformalMapID;
};

#endif
