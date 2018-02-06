#include "GLUTSurfaceWindow.h"
#include "PipelineGeneral.h"

GLUTSurfaceWindow * GLUTSurfaceWindow::s_MainGLUTWindow = NULL;

////////// Window memeber functions ///////////// 

GLUTSurfaceWindow::GLUTSurfaceWindow(ParamParser & drawParams, 
									 AnalysisPipeline * pipeline,
									 int * argc, char ** argv)
: AnalysisWindow(drawParams, pipeline)
{
	m_iFeatureID = -1;
	m_iSampleSetID = -1;
	
	s_MainGLUTWindow = this;
	m_GLUTwindowID = 0;
	m_GLUTmouse[0] = 0;		m_GLUTmouse[1] = 0;	
	m_GLUTbutton[0] = 0;	m_GLUTbutton[1] = 0;	m_GLUTbutton[1] = 0;
	m_GLUTmodifiers = 0;
	
	// Open window 
	glutInit(argc, argv);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(Width(), Height());
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_STENCIL
	m_GLUTwindowID = glutCreateWindow("Surface Viewer (GLUT)");
	
	InitializeGL();	
	
	// Initialize GLUT callback functions 
	glutDisplayFunc(GLUTRedraw);
	glutReshapeFunc(GLUTResize);
	glutKeyboardFunc(GLUTKeyboard);
	glutSpecialFunc(GLUTSpecial);
	glutMouseFunc(GLUTMouse);
	glutMotionFunc(GLUTMotion);
	
	glutMainLoop();
}

void GLUTSurfaceWindow::RenderText(R3Point p, const VKString & text, 
								   double r, double g, double b, TextSize size)
{
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	
	// Draw text string s and position p
	glColor3d(r, g, b);
	glRasterPos3d(p[0], p[1], p[2]);

	const char * s = text.c_str();
	if (size==TEXT_SIZE_SMALL)
		while (*s) glutBitmapCharacter(	GLUT_BITMAP_HELVETICA_10 , *(s++));	
	else if (size==TEXT_SIZE_MEDIUM)
		while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *(s++));	
	else if (size==TEXT_SIZE_LARGE)
		while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *(s++));	
}


void GLUTSurfaceWindow::Terminate()
{
	// Destroy window 
	glutDestroyWindow(m_GLUTwindowID);
	
	// Exit
	exit(0);	
}
void GLUTSurfaceWindow::Update()
{
	glutPostRedisplay();
}

///////// GLUT Callback /////////

void GLUTSurfaceWindow::GLUTRedrawCallback()
{
	Draw();
	// Swap buffers 
	glutSwapBuffers();		
}

void GLUTSurfaceWindow::GLUTResizeCallback(int w, int h)
{
	Resize(w, h);
	// Redraw
	glutPostRedisplay();
}

void GLUTSurfaceWindow::GLUTMotionCallback(int x, int y)
{
	// Invert y coordinate
	y = Height() - y;
	
	// Compute mouse movement
	int dx = x - m_GLUTmouse[0];
	int dy = y - m_GLUTmouse[1];
	
	bool shift = (m_GLUTmodifiers & GLUT_ACTIVE_SHIFT) != 0;
	bool ctrl =  (m_GLUTmodifiers & GLUT_ACTIVE_CTRL) != 0;
	bool iCor = m_DrawParams.GetStrValue("RendererDefault", "InteractionMode", valid)=="ManualCorrespondences";
	
	bool done = false;
	if (iCor)
	{
		bool valid;
		VKStringList sampleDescr = m_DrawParams.GetStrValues("RendererDefault", "InteractingSample", valid);
		if (valid && sampleDescr.count()==3)
		{
			VKString typeObject = sampleDescr[0];
			VKString objName = sampleDescr[1];
			int sampleID = sampleDescr[2].toInt(&valid);
			assert(valid);
			SurfaceSampleSet * sampleSet = NULL;
			VKString surfaceName, textureName;
			if (typeObject=="Texture")
			{
				const VKString & setname = objName+s_InteractiveSampleSetName;
				sampleSet = Surface2DPlane::m_pInfinitePlanePseudosurface->GetSampleSet(setname);
				bool inWindow;
				double xWindow, yWindow;
				VKString windowName = InSupplementalWindow(x, y, &inWindow, &xWindow, &yWindow);
				if (inWindow && windowName=="Textures")
					sampleSet->ReplaceSample(sampleID, SurfaceSample(-1, xWindow, yWindow, 0,
																	 Surface2DPlane::m_pInfinitePlanePseudomesh));
				else
					std::cout<<"[WARNING] Moved out of window: "<<windowName.c_str()<<std::endl;
			}
			else if (typeObject=="Surface")
			{
				SampledSurface * surface = PipelineGeneral::s_pGeneralPipeline->GetSurfaceByName(objName);
				assert (surface!=NULL);
				sampleSet = surface->GetSampleSet(s_InteractiveSampleSetName);
				int vertexID = FindSample(x, y, "AllVertices", surfaceName, textureName);
				if (vertexID>=0)
					sampleSet->ReplaceSample(sampleID, SurfaceSample(vertexID, surface->GetMesh()));
			}
			done = true;
			glutPostRedisplay();
			// TODO: move map update to fn
			KeyboardCallbackGeneral('u', 0, 0);
		}
	}
	
	bool inWindow;
	double xWindow, yWindow;
	VKString windowContent = InSupplementalWindow(x, y, &inWindow, &xWindow, &yWindow);
	if (inWindow && !done)
	{
//		if (windowContent=="Clustering")
//			SelectMapInClusteringWindow(xWindow, yWindow, shift);
//		else if (iCor && windowContent=="Textures")
//			SelectAndMovePointInTexturesWindow(xWindow, yWindow);
			
		glutPostRedisplay();
	}
	else if (!done)
	{	
		// World in hand navigation  
		if (m_GLUTbutton[0] && !ctrl && !shift) 
			s_MainGLUTWindow->RotateCurrentMesh(x, y, dx, dy);
		else if (m_GLUTbutton[1] || (m_GLUTbutton[0] && ctrl && !shift)) 
			s_MainGLUTWindow->ScaleCurrentMesh(x, y, dx, dy);
		else if (m_GLUTbutton[2] || (m_GLUTbutton[0] && !ctrl && shift)) 
			s_MainGLUTWindow->TranslateCurrentMesh(x, y, dx, dy);
		else if (m_GLUTbutton[2] || (m_GLUTbutton[0] && shift && ctrl))
			s_MainGLUTWindow->RotateLights(1, x, y, dx, dy);
//		else if (m_GLUTbutton[2] || (m_GLUTbutton[0] && !shift && !ctrl)) 
//			s_MainGLUTWindow->RotateLights(0, x, y, dx, dy);

		
		if (m_GLUTbutton[0] || m_GLUTbutton[1] || m_GLUTbutton[2]) 
			glutPostRedisplay();
	}
	// Remember mouse position 
	m_GLUTmouse[0] = x;
	m_GLUTmouse[1] = y;
}

void GLUTSurfaceWindow::GLUTMouseCallback(int button, int state, int x, int y)
{
	bool valid;
	// Invert y coordinate
	y = Height() - y;
	bool shift = (glutGetModifiers() & GLUT_ACTIVE_SHIFT) != 0;
//	bool ctrl =  (glutGetModifiers() & GLUT_ACTIVE_CTRL) != 0;	
		
	// Remember button state 
	int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
	m_GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;
	
	// Remember modifiers 
	m_GLUTmodifiers = glutGetModifiers();
	
	// Remember mouse position 
	m_GLUTmouse[0] = x;
	m_GLUTmouse[1] = y;
	
	bool iCor = m_DrawParams.GetStrValue("RendererDefault", "InteractionMode", valid)=="ManualCorrespondences";
	
	bool done = false;
	if (iCor)
	{
		if (b==0 && m_GLUTbutton[b]==1)	// press down
		{
			VKString surfaceName, textureName;
			int sampID = FindSample(x, y, s_InteractiveSampleSetName, surfaceName, textureName);
			if (sampID >= 0)
			{
				m_DrawParams.m_ParamModules[m_RendererName]["InteractingSample"].clear();
				if (surfaceName!="")
				{
					m_DrawParams.m_ParamModules[m_RendererName]["InteractingSample"].push_back("Surface");
					m_DrawParams.m_ParamModules[m_RendererName]["InteractingSample"].push_back(surfaceName);
				}
				else
				{
					assert(textureName!="none");
					m_DrawParams.m_ParamModules[m_RendererName]["InteractingSample"].push_back("Texture");
					m_DrawParams.m_ParamModules[m_RendererName]["InteractingSample"].push_back(textureName);
				}
				m_DrawParams.m_ParamModules[m_RendererName]["InteractingSample"].push_back(sampID);
				done = true;
			}
		}
		else if (b==0 && m_GLUTbutton[b]==0)	// release
		{
			m_DrawParams.m_ParamModules[m_RendererName]["InteractingSample"].clear();
			done = true;
		}
	}
	
	// Process mouse button event
//	bool inWindow;
//	double xWindow, yWindow;
//	VKString windowContent = InSupplementalWindow(x, y, &inWindow, &xWindow, &yWindow);
//	if (inWindow && !done)
//	{
//		if (windowContent=="Clustering")		
//			SelectMapInClusteringWindow(xWindow, yWindow, shift);
//	}		
	//else...
	
	
	// Redraw
	glutPostRedisplay();
}

void GLUTSurfaceWindow::GLUTKeyboardCallback(unsigned char key, int x, int y)
{
	// Invert y coordinate
	y = Height() - y;
	
	if (m_DrawParams.GetStrValue("RendererDefault", "Type", valid)=="SurfaceCorrespondence")
		KeyboardCallbackCorrespondence(key, x, y);
	else
		KeyboardCallbackGeneral(key, x, y);
	
	// Remember mouse position 
	m_GLUTmouse[0] = x;
	m_GLUTmouse[1] = y;
	
	// Remember modifiers 
	m_GLUTmodifiers = glutGetModifiers();
	
	// Redraw
	glutPostRedisplay();  	
}

void GLUTSurfaceWindow::KeyboardCallbackGeneral(unsigned char key, int x, int y)
{
	PipelineGeneral * pipeline = PipelineGeneral::s_pGeneralPipeline;
	switch (key) 
	{
		case '`':
			std::cout<<"Selecting Elements. Iterate over..."<<std::endl;
			std::cout<<"\t\t [ - surfaces"<<std::endl;
			std::cout<<"\t\t p - textures"<<std::endl;
			std::cout<<"\t\t . - samples"<<std::endl;
			std::cout<<"\t\t ' - surface maps"<<std::endl;
			std::cout<<"\t\t\t ; - `submaps` (if current map is a collection)"<<std::endl;;
			std::cout<<"\t\t l - surface maps (for similarity)"<<std::endl;
			std::cout<<"\t\t\t k - `submaps` (for similarity) (in collection)"<<std::endl;
			std::cout<<"\t\t ] - features"<<std::endl;
			std::cout<<"\t\t / - map confidence/distortion"<<std::endl;
			std::cout<<"\t\t , - map (dis?)similarites"<<std::endl;
			std::cout<<"\t\t m - distances"<<std::endl;
			std::cout<<"Selecting Signal To Render:"<<std::endl;
			std::cout<<"\t\t ! - Map Confidence"<<std::endl;			
			std::cout<<"\t\t @ - Map XYZ rendering"<<std::endl;
			std::cout<<"\t\t # - Map Far Vertices rendering"<<std::endl;
			std::cout<<"\t\t $ - Map Similarity"<<std::endl;
			std::cout<<"\t\t % - Distance"<<std::endl;
			std::cout<<"\t\t ^ - Texture"<<std::endl;
			std::cout<<"\t\t & - Feature"<<std::endl;
			std::cout<<"\t\t * - Auto"<<std::endl;
			std::cout<<"Other Rendering Options:"<<std::endl;
			std::cout<<"\t\t S - toggle rendering samples vs sample IDs"<<std::endl;
			std::cout<<"Interaction:"<<std::endl;
			std::cout<<"\t\t i - change interaction mode"<<std::endl;
			std::cout<<"\t\t SPACE - select a point on a surface"<<std::endl;
			std::cout<<"\t\t u - update current mesh using interactive correspondences"<<std::endl;
			std::cout<<"Print names:"<<std::endl;
			std::cout<<"\t\t 1 - Confdences"<<std::endl;			
			std::cout<<"\t\t 2 - Maps"<<std::endl;
			std::cout<<"\t\t 3 - Surfaces"<<std::endl;
			std::cout<<"\t\t 4 - Similarities"<<std::endl;
			std::cout<<"\t\t 5 - Distances"<<std::endl;
			std::cout<<"\t\t 6 - Textures"<<std::endl;
			std::cout<<"\t\t 7 - Features"<<std::endl;
			std::cout<<"\t\t 8 - Sample Sets"<<std::endl;
			break;
		case '1':
			PrintList("Map Confidences", pipeline->GetConfidenceNames());
			break;
		case '2':
			PrintList("Maps", pipeline->GetMapNames());
			break;
		case '3':
			PrintList("Surfaces", pipeline->GetSurfaceNames());
			break;
		case '4':
			PrintList("Map Similarities", pipeline->GetSimilarityNames());
			break;
		case '5':
			PrintList("Surface Distance Metrics", pipeline->GetDistanceMetricNames());
			break;
		case '6':
			PrintList("Surfaces", pipeline->GetTextureNames());
			break;
		case '7':
			PrintList("Surface Features", pipeline->GetFeatureNames());
			break;
		case '8':
			PrintList("Surface Sample Sets", pipeline->GetSampleSetNames());
			break;
		case '[':
			GenPipeline_UpdateRenderField("RenderSurface", 1);
			break;
		case '{':
			GenPipeline_UpdateRenderField("RenderSurface", -1);
			break;
		case 'p':
			GenPipeline_UpdateRenderField("RenderTexture", 1);
			break;
		case 'P':
			GenPipeline_UpdateRenderField("RenderTexture", -1);
			break;
		case '.':
			GenPipeline_UpdateRenderField("RenderSampleSet", 1);
			break;
		case '>':
			GenPipeline_UpdateRenderField("RenderSampleSet", -1);
			break;
		case '\'':
			GenPipeline_UpdateRenderField("RenderMap", 1);
			break;
		case '"':
			GenPipeline_UpdateRenderField("RenderMap", -1);
			break;
		case ';':
		case ':':
		{
			VKStringList & subMapID = m_DrawParams.m_ParamModules[m_RendererName]["RenderSubmap"];
			assert(subMapID.count()==1);
			bool ok;
			int intSubmapID = subMapID[0].toInt(&ok);
			assert(ok);
			if (key==';')
				intSubmapID++;
			else
				intSubmapID--;
			subMapID.clear();
			subMapID.push_back(VKString(intSubmapID));
			break;
		}
		case 'l':
			GenPipeline_UpdateRenderField("RenderOtherMap", 1);
			break;
		case 'L':
			GenPipeline_UpdateRenderField("RenderOtherMap", -1);
			break;
		case 'k':
		case 'K':
			std::cout<<"TODO: select surface submaps in collection 2"<<std::endl;
			break;
		case ']':
			GenPipeline_UpdateRenderField("RenderFeatureValue", 1);
			break;
		case '}':
			GenPipeline_UpdateRenderField("RenderFeatureValue", -1);
			break;
		case '/':
			GenPipeline_UpdateRenderField("RenderConfidence", 1);
			break;
		case '?':
			GenPipeline_UpdateRenderField("RenderConfidence", -1);
			break;
		case ',':
			GenPipeline_UpdateRenderField("RenderSimilarity", 1);
			break;
		case '<':
			GenPipeline_UpdateRenderField("RenderSimilarity", -1);
			break;			
		case 'm':
			GenPipeline_UpdateRenderField("RenderDistanceMetric", 1);
			break;
		case 'M':
			GenPipeline_UpdateRenderField("RenderDistanceMetric", -1);
			break;
		case ' ':
		{
			VKString inter = m_DrawParams.GetStrValue("RendererDefault", "InteractionMode", valid);
			if (inter=="none")
				SelectSample(x, y);
			else if (inter=="ManualCorrespondences")
				AddCorrespondence(x, y);
			break;			
		}
		case 'S':
			FlipFlag("RendererDefault", "SurfFlags", "SampleIDs");
			break;
		case 'e':
			FlipFlag(m_RendererName, "MeshFlags", "Edges");
			break;
		case 'f':
			FlipFlag(m_RendererName, "MeshFlags", "Faces");
			break;
		case 'i':
		{
			VKStringList & interactMode = m_DrawParams.m_ParamModules[m_RendererName]["InteractionMode"];
			assert(interactMode.count()==1);
			if (interactMode[0] == "none")
			{
				interactMode.clear();
				interactMode.push_back("ManualCorrespondences");
				
			}
			else
			{
				interactMode.clear();
				interactMode.push_back("none");
			}
			std::cout<<"Interaction Mode: "<<interactMode[0].c_str()<<std::endl;	
			break;
		}
		case 'u':
		{
			if (m_DrawParams.GetStrValue("RendererDefault", "InteractionMode", valid)=="ManualCorrespondences")
			{	
				VKString surfName = m_DrawParams.GetStrValue("RendererDefault", "RenderSurface", valid);
				SampledSurface * surf = PipelineGeneral::s_pGeneralPipeline->GetSurfaceByName(surfName);
				if (surf==NULL)
				{
					std::cout<<"[WARNING] Cannot find surface: "<<surfName.c_str()<<std::endl;
					return;
				}
				
				VKString mapName = m_DrawParams.GetStrValue("RendererDefault", "RenderMap", valid);
				SurfaceMap * surfMap = surf->GetMapActingAsFrom(mapName);
				if (surfMap==NULL)
				{
					std::cout<<"[WARNING] Surface "<<surfName.c_str()<<" does not have map "<<mapName.c_str()<<std::endl;;
					return;
				}
				VKString textureName = m_DrawParams.GetStrValue("RendererDefault", "RenderTexture", valid);
				surfMap->InteractiveUpdateFromCorrs(AnalysisWindow::s_InteractiveSampleSetName,
													textureName);
				//std::cout<<"Updating map: "<<mapName.c_str()<<" for surface "<<surfName.c_str()<<std::endl;
				//std::cout<<"with texturename: "<<textureName.c_str();
				//std::cout<<" samples="<<AnalysisWindow::s_InteractiveSampleSetName.c_str()<<std::endl;
			}
			else
				std::cout<<"[WARNING] Change interaction mode"<<std::endl;
			break;
		}
		case 'c':
			if (m_DrawParams.GetStrValue("RendererDefault", "InteractionMode", valid)=="ManualCorrespondences")
			{
				VKString surfName = m_DrawParams.GetStrValue("RendererDefault", "RenderSurface", valid);
				SampledSurface * surf = PipelineGeneral::s_pGeneralPipeline->GetSurfaceByName(surfName);
				if (surf==NULL)
				{
					std::cout<<"[WARNING] Cannot find surface: "<<surfName.c_str()<<std::endl;
					return;
				}
				VKString textureName = m_DrawParams.GetStrValue("RendererDefault", "RenderTexture", valid);
				SurfaceSampleSet * sampSet = surf->GetSampleSet(AnalysisWindow::s_InteractiveSampleSetName);
				if (sampSet==NULL)
				{
					std::cout<<"[WARNING] Cannot find sample set on surface: ";
					std::cout<<AnalysisWindow::s_InteractiveSampleSetName.c_str()<<std::endl;
					return;
				}
				sampSet->ClearSamples();
				sampSet = Surface2DPlane::m_pInfinitePlanePseudosurface->GetSampleSet(textureName+AnalysisWindow::s_InteractiveSampleSetName);
				if (sampSet==NULL)
				{
					std::cout<<"[WARNING] Cannot find sample set on texture: ";
					std::cout<<AnalysisWindow::s_InteractiveSampleSetName.c_str()<<std::endl;
					return;
				}
				sampSet->ClearSamples();
			}
			else
				std::cout<<"[WARNING] Change interaction mode"<<std::endl;
			break;
		default:
			std::cout<<"[WARNING] Unknown character. Press ` to see help"<<std::endl;
			break;
	}
}

void GLUTSurfaceWindow::KeyboardCallbackCorrespondence(unsigned char key, int x, int y)
{
	// Process keyboard button event 
	switch (key) {
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
			if (key == '0') SetCurrentMesh(-1);
			else 
				SetCurrentMesh(key-'1');
			break;			
		case '`':			
			std::cout<<"Special Chars:"<<std::endl;
			std::cout<<"\tSurface Signals: "<<std::endl;
			std::cout<<"\t\t ! - MapConfidence"<<std::endl;
			std::cout<<"\t\t @ - XYZ"<<std::endl;
			std::cout<<"\t\t # - FarVertices"<<std::endl;
			std::cout<<"\t\t $ - MapDiscrepancy"<<std::endl;
			std::cout<<"\t\t % - Distance"<<std::endl;
			std::cout<<"\t\t Q - MultiMapDistGen"<<std::endl;
			std::cout<<"\t\t W - MultiMapConsistency"<<std::endl;
			std::cout<<"\t\t E - MultiMapConfidence"<<std::endl;
			std::cout<<"\t\t R - MultiMapFinal"<<std::endl;
			std::cout<<"\t\t T - MultiMapTop3Weights"<<std::endl;			
			std::cout<<"\tSingle Map Colors: "<<std::endl;			
			std::cout<<"\t\t + - Consistency"<<std::endl;
			std::cout<<"\t\t - - MapConfidence"<<std::endl;
			std::cout<<"\t\t ) - DistanceToGenerators"<<std::endl;
			std::cout<<"\t\t ( - CombinedWeight"<<std::endl;			
			std::cout<<"\tMore:"<<std::endl;			
			std::cout<<"\t\t C - Advance Conformal Map in Multi-conformal map"<<std::endl;
			break;
		case '+':
		case '_':
		case ')':
		case '(':
		{
			VKStringList & singleCorrColor = m_DrawParams.m_ParamModules[m_RendererName]["SingleCorrColor"];
			if (singleCorrColor.contains("none"))
			{
				singleCorrColor.clear();
				if (key=='+')
					singleCorrColor.push_back("Consistency");
				else if (key=='_')
					singleCorrColor.push_back("MapConfidence");
				else if (key==')')
					singleCorrColor.push_back("DistanceToGenerators");
				else if (key=='(')
					singleCorrColor.push_back("CombinedWeight");
				else
					assert(false);
				
				std::cout<<"Enabled "<<singleCorrColor[0].c_str()<<" Rendering"<<std::endl;				
			}
			else
			{
				singleCorrColor.clear();
				singleCorrColor.push_back("none");
				std::cout<<"Disabled Single Corr Colors"<<std::endl;				
			}
			break;
		}	
		case '!':
		case '@':
		case '#':	
		case '$':
		case '%':	
		case '^':	
		case '&':
		case 'Q':
		case 'W':
		case 'E':
		case 'R':
		case 'T':
		{
			VKStringList & drawSignal = m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"];
			if (drawSignal.contains("none"))
			{
				drawSignal.clear();
				if (key=='!')
					drawSignal.push_back("MapConfidence");
				else if (key=='@')
					drawSignal.push_back("XYZ");
				else if (key=='#')
					drawSignal.push_back("FarVertices");
				else if (key=='$')
					drawSignal.push_back("MapDiscrepancy");
				else if (key=='%')
					drawSignal.push_back("Distance");
				else if (key=='^')
					drawSignal.push_back("ThresholdedDistance");
				else if (key=='&')
					drawSignal.push_back("ErrorColors");
				else if (key=='Q')
					drawSignal.push_back("MultiMapDistGen");
				else if (key=='W')
					drawSignal.push_back("MultiMapConsistency");
				else if (key=='E')
					drawSignal.push_back("MultiMapConfidence");
				else if (key=='R')
					drawSignal.push_back("MultiMapFinal");
				else if (key=='T')
					drawSignal.push_back("MultiMapTop3Weights");
				else
					assert(false);
				
				std::cout<<"Enabled "<<drawSignal[0].c_str()<<" Rendering"<<std::endl;				
			}
			else
			{
				drawSignal.clear();
				drawSignal.push_back("none");
				std::cout<<"Disabled Map On Surface Rendering"<<std::endl;				
			}
			break;
		}
		case 'C':
		{
			int confMap = m_DrawParams.GetIntValue(m_RendererName, "ConfMapInMultiMap", valid);
			assert(valid);
			confMap++;
			m_DrawParams.m_ParamModules[m_RendererName]["ConfMapInMultiMap"].clear();
			m_DrawParams.m_ParamModules[m_RendererName]["ConfMapInMultiMap"].push_back(VKString::number(confMap));
			break;
		}
			//		case 'B':
			//		case 'b':
			//			FlipFlag(m_RendererName, "MeshFlags", "Backfacing");
			//			break;
			
		case 'p': 
			PrintCamera();
			break;			
			//		case 'E':
			//			FlipFlag(m_RendererName, "MeshFlags", "EdgeIDs");
			//			break;			
		case 'e':
			FlipFlag(m_RendererName, "MeshFlags", "Edges");
			break;
			
		case 'F':
			FlipFlag(m_RendererName, "MeshFlags", "FaceIDs");
			break;
			
		case 'f':
			FlipFlag(m_RendererName, "MeshFlags", "Faces");
			break;
			
		case 'V':
			FlipFlag(m_RendererName, "MeshFlags", "VertexIDs");
			break;
		case 'v':
			FlipFlag(m_RendererName, "MeshFlags", "Vertices");
			break;
		case 'S':
			FlipFlag(m_RendererName, "SurfFlags", "SampleIDs");
			break;
		case ';':
		{
			if (m_pPipeline->GetCurrCollection()!=NULL)
				m_pPipeline->GetCurrCollection()->IncreaseSelectedMapID();
			int mapID = m_pPipeline->GetCurrCollection()->GetSelectedMapID();
			std::cout<<"Map #"<<mapID<<" Score="<<m_pPipeline->GetCurrCollection()->GetMapValue(mapID)<<std::endl;
			break;
		}
		case ':':
		{
			if (m_pPipeline->GetCurrCollection()!=NULL)
				m_pPipeline->GetCurrCollection()->DecreaseSelectedMapID();
			int mapID = m_pPipeline->GetCurrCollection()->GetSelectedMapID();
			std::cout<<"Map #"<<mapID<<" Score="<<m_pPipeline->GetCurrCollection()->GetMapValue(mapID)<<std::endl;
			break;
		}
		case '\'':
		{
			VKStringList & renderMaps = m_DrawParams.m_ParamModules[m_RendererName]["RenderMap"];
			if (renderMaps.contains("BestConformal"))
			{
				renderMaps.clear();
				renderMaps.push_back("BestExtrapolated");
			}
			else if (renderMaps.contains("BestExtrapolated"))
			{
				renderMaps.clear();
				//renderMaps.push_back("BenchmarkQuery");
				renderMaps.push_back("Truth");
			}
			else if (renderMaps.contains("BenchmarkQuery"))
			{
				renderMaps.clear();
				renderMaps.push_back("Truth");
			}
			else if (renderMaps.contains("Truth"))
			{
				renderMaps.clear();
				renderMaps.push_back("Collection");
			}
			else if (renderMaps.contains("Collection"))
				renderMaps.clear();
			else 
			{
				renderMaps.clear();
				renderMaps.push_back("BestConformal");
			}
			if (renderMaps.count()==0)
				std::cout<<"Disabled Map Rendering. "<<std::endl;
			else
				std::cout<<"Rendering Map: "<<renderMaps[0].c_str()<<std::endl;
			break;
		}
		case 'm':
			if (m_DrawParams.GetStrValue(m_RendererName, "ConfSurf", valid)=="RenderOriginal")
			{
				VKStringList & tempList = m_DrawParams.m_ParamModules[m_RendererName]["ConfSurf"];
				tempList.clear();
				tempList.push_back("RenderConformal");
				SetCurrentMeshInCollection(1);				
			}
			else
			{
				VKStringList & tempList = m_DrawParams.m_ParamModules[m_RendererName]["ConfSurf"];
				tempList.clear();
				tempList.push_back("RenderOriginal");
				SetCurrentMeshInCollection(0);			
			}
			break;
		case 'M':
			if (m_DrawParams.GetStrValue(m_RendererName, "ConfSurf", valid)=="RenderConformalMap1")
			{
				VKStringList & tempList = m_DrawParams.m_ParamModules[m_RendererName]["ConfSurf"];
				tempList.clear();
				tempList.push_back("RenderConformalMap2");
				SetCurrentMeshInCollection(1);
			}			
			else
			{
				VKStringList & tempList = m_DrawParams.m_ParamModules[m_RendererName]["ConfSurf"];
				tempList.clear();
				tempList.push_back("RenderConformalMap1");
				SetCurrentMeshInCollection(1);
			}			
			break;
		case 'c':
			FlipFlag(m_RendererName, "MapFlags", "Edges");
			break;
		case '[':	// feature browsing
		{
			VKStringList featureNames = m_pPipeline->GetSurface(m_iCurrentMesh)->GetAllFeatureNames();
			if (m_iFeatureID<0)			
				m_iFeatureID = featureNames.count()-1;
			else
				m_iFeatureID--;
			UpdateFeature();
			break;
		}
		case ']':		// feature browsing	
		{
			VKStringList featureNames = m_pPipeline->GetSurface(m_iCurrentMesh)->GetAllFeatureNames();
			if (m_iFeatureID>=featureNames.count()-1)	
				m_iFeatureID = -1;
			else
				m_iFeatureID++;
			UpdateFeature();
			break;			
		}
		case ',':	// sample set browsing
		case '<':
		{
			VKStringList sampleSetNames = m_pPipeline->GetSurface(m_iCurrentMesh)->GetAllSampleSetNames();
			if (m_iSampleSetID<0)			
				m_iSampleSetID = sampleSetNames.count()-1;
			else
				m_iSampleSetID--;
			UpdateSampleSet();
			break;			
		}
		case '.':
		case '>':
		{
			VKStringList sampleSetNames = m_pPipeline->GetSurface(m_iCurrentMesh)->GetAllSampleSetNames();
			if (m_iSampleSetID>=sampleSetNames.count()-1)	
				m_iSampleSetID = -1;
			else
				m_iSampleSetID++;
			UpdateSampleSet();
			break;
		}		
		case 'A':
			//		assert(m_pPipeline->GetPipelineType()=="PipelineExploreConfmaps");
			((PipelineExploreConfmaps*)m_pPipeline)->PrintIntegratedBlendingMatrix(&m_DrawParams);
			break;
		case 'a':
			//		assert(m_pPipeline->GetPipelineType()=="PipelineExploreConfmaps");
			SelectSample(x, y);
			((PipelineExploreConfmaps*)m_pPipeline)->MapSelectedBasedOnIntegratedBlendedMatrix(&m_DrawParams);
			break;
		case 'n':
			//		assert(m_pPipeline->GetPipelineType()=="PipelineExploreConfmaps");
			((PipelineExploreConfmaps*)m_pPipeline)->PickNextNplet();
			break;
			//		case 'b':
			//assert(m_pPipeline->GetPipelineType()=="PipelineExploreConfmaps");			
			//((PipelineExploreConfmaps*)m_pPipeline)->VoteForCurrentConfMap();
			//std::cout<<((PipelineExploreConfmaps*)m_pPipeline)->GetCurrentVoteDescriptor().c_str()<<std::endl;
			//			break;
		case ' ':
			SelectSample(x, y);
			break;
		case 27: // ESCAPE
			Terminate();
			break;
	}
}

void GLUTSurfaceWindow::UpdateFeature()
{
	if (m_iFeatureID<0)
	{
		std::cout<<"Disabled Feature Rendering"<<std::endl;
		m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"].clear();
		m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"].push_back("none");
	}
	else
	{
		m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"].clear();
		m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"].push_back("Feature");
		for (int i=0; i<m_pPipeline->NumSurfaces(); i++)
		{
			VKStringList featureNames = m_pPipeline->GetSurface(i)->GetAllFeatureNames();				
			assert(m_iFeatureID < featureNames.count());			
			if (i==0)
				std::cout<<"Rendering Feature: "<<featureNames[m_iFeatureID].c_str()<<std::endl;
			m_pPipeline->GetSurface(i)->SelectFeature(featureNames[m_iFeatureID]);
		}
	}
}

void GLUTSurfaceWindow::UpdateSampleSet()
{
	if (m_iSampleSetID<0)
	{
		std::cout<<"Disabled Sample Rendering"<<std::endl;
		SetFlag(m_RendererName, "SurfFlags", "DefaultSamples", false);
	}
	else
	{
		SetFlag(m_RendererName, "SurfFlags", "DefaultSamples", true);
		for (int i=0; i<m_pPipeline->NumSurfaces(); i++)
		{
			VKStringList sampleSetNames = m_pPipeline->GetSurface(i)->GetAllSampleSetNames();				
			assert(m_iSampleSetID < sampleSetNames.count());
			if (i==0)
				std::cout<<"Rendering Samples: "<<sampleSetNames[m_iSampleSetID].c_str()<<std::endl;			
			m_pPipeline->GetSurface(i)->SelectSampleSet(sampleSetNames[m_iSampleSetID]);
		}
		//TODO: do the same for symmetric mesh
	}
}


void GLUTSurfaceWindow::GLUTSpecialCallback(int key, int x, int y)
{
	// Invert y coordinate
	y = Height() - y;
	
	// Process keyboard button event 
	
	// Remember mouse position 
	m_GLUTmouse[0] = x;
	m_GLUTmouse[1] = y;
	
	// Remember modifiers 
	m_GLUTmodifiers = glutGetModifiers();
	
	// Redraw
	glutPostRedisplay();
}



///////// STATIC GLUT FUNCTIONS /////////
void GLUTRedraw()
{
	GLUTSurfaceWindow::s_MainGLUTWindow->GLUTRedrawCallback();
}

void GLUTResize(int w, int h)
{
	GLUTSurfaceWindow::s_MainGLUTWindow->GLUTResizeCallback(w, h);
}

void GLUTMotion(int x, int y)
{
	GLUTSurfaceWindow::s_MainGLUTWindow->GLUTMotionCallback(x, y);
}

void GLUTMouse(int button, int state, int x, int y)
{
	GLUTSurfaceWindow::s_MainGLUTWindow->GLUTMouseCallback(button, state, x, y);
}

void GLUTSpecial(int key, int x, int y)
{
	GLUTSurfaceWindow::s_MainGLUTWindow->GLUTSpecialCallback(key, x, y);
}

void GLUTKeyboard(unsigned char key, int x, int y)
{
	GLUTSurfaceWindow::s_MainGLUTWindow->GLUTKeyboardCallback(key, x, y);
}



