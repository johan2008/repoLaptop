#include "AnalysisWindow.h"
#include "AnalysisPipeline.h"
#include "PipelineExploreConfmaps.h"
#include "PipelineGeneral.h"

AnalysisWindow * AnalysisWindow::s_pMainWindow = NULL;

const VKString AnalysisWindow::s_InteractiveSampleSetName = "InteractiveSampleSet";

/**
 * Creates a new UI window.
 * @param drawParams	rendering parameters
 * @param pipeline		actual algorithm for surface analysis
 * @param rendererName	name of a renderer to be used (in drawParams)
 */
AnalysisWindow::AnalysisWindow(ParamParser & drawParams, AnalysisPipeline * pipeline, 
							   const VKString & rendererName)
: m_DrawParams(drawParams)
{
	s_pMainWindow = this;
	m_AvailiableSubWindows.push_back("TopRightWindow");
	m_AvailiableSubWindows.push_back("TopLeftWindow");
	m_AvailiableSubWindows.push_back("BottomRightWindow");
	m_AvailiableSubWindows.push_back("BottomLeftWindow");
	
	m_iCurrentCameraCollectionID = 0;
	m_iCurrentCameraMeshID = 0;
	
	m_iCurrentMesh = 0;
	m_iCurrentSubMesh = 0;
	m_pPipeline = pipeline;
	m_RendererName = rendererName;	
	m_pPipeline->SetWindow(this);
	
	m_MulticolorSequence.push_back(RNRgb(0, 0, 1));
	m_MulticolorSequence.push_back(RNRgb(.6, .6, 1));
	m_MulticolorSequence.push_back(RNRgb(0, 1, 0));
	m_MulticolorSequence.push_back(RNRgb(1, 1, 0));
	m_MulticolorSequence.push_back(RNRgb(1, 0, 0));	

//	std::cout<<"Cluster to analyze"<<std::endl;
//	VKString numClusters = VKString::number(m_pPipeline->ClustersToAnalyze());
//	std::cout<<"Cluster to analyze: "<<numClusters.c_str()<<std::endl;	
//	m_DrawParams.m_ParamModules["RendererDefault"]["VisCollectItems"].push_back(numClusters);
}

/**
 * Draws an openGL context. Calls to the AnalysisPipeline object. Creates snapshots.
 */
void AnalysisWindow::Draw()
{
	// Clear window 
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	bool takingShot=false;
	VKStringList restoreVal;
	VKString screenshotName="";
	VKString colors;
	VKStringList mapFlags;
	VKString maxNumPnts;
	if (m_DrawParams.GetStrValues(m_RendererName, "ShotTypes", valid).contains("InitialEmptyShot"))
		screenshotName = "InitialEmptyShot";
	else if (m_DrawParams.GetStrValues(m_RendererName, "ShotTypes", valid).contains("BestExtrapolated"))
		screenshotName="BestExtrapolated";
	else if (m_DrawParams.GetStrValues(m_RendererName, "ShotTypes", valid).contains("BenchmarkQuery"))
		screenshotName="BenchmarkQuery";
	else if (m_DrawParams.GetStrValues(m_RendererName, "ShotTypes", valid).contains("BestConformal"))
		screenshotName="BestConformal";
	else if (m_DrawParams.GetStrValues(m_RendererName, "ShotTypes", valid).contains("FinalColors"))
		screenshotName="FinalColors";
	else if (m_DrawParams.GetStrValues(m_RendererName, "ShotTypes", valid).contains("ErrorColors"))
		screenshotName="ErrorColors";
	else if (m_DrawParams.GetStrValue(m_RendererName, "ShotOnly", valid)=="true")
		Terminate();
	
	if (screenshotName!="")
	{
		takingShot = true;
		m_DrawParams.m_ParamModules[m_RendererName]["ShotTypes"].remove(screenshotName);
		restoreVal = m_DrawParams.GetStrValues(m_RendererName, "RenderMap", valid);
		m_DrawParams.m_ParamModules[m_RendererName]["RenderMap"].clear();
		m_DrawParams.m_ParamModules[m_RendererName]["RenderMap"].push_back(screenshotName);
		if (screenshotName=="FinalColors")
		{
			m_DrawParams.m_ParamModules[m_RendererName]["RenderMap"].clear();
			m_DrawParams.m_ParamModules[m_RendererName]["RenderMap"].push_back("BestExtrapolated");
			
			colors=m_DrawParams.GetStrValue(m_RendererName, "DrawMapColorsOnSurf", valid);
			m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"].clear();
			m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"].push_back("XYZ");
			mapFlags = m_DrawParams.GetStrValues(m_RendererName, "MapFlags", valid);
			
			m_DrawParams.m_ParamModules[m_RendererName]["MapFlags"].clear();
		}
		else if (screenshotName=="BestExtrapolated")
		{
			maxNumPnts=m_DrawParams.GetStrValue(m_RendererName, "MaxNumMapColors", valid);
			m_DrawParams.m_ParamModules[m_RendererName]["MaxNumMapColors"].clear();
			m_DrawParams.m_ParamModules[m_RendererName]["MaxNumMapColors"].push_back("60");
		}
		else if (screenshotName=="ErrorColors")
		{
			m_DrawParams.m_ParamModules[m_RendererName]["RenderMap"].clear();
			m_DrawParams.m_ParamModules[m_RendererName]["RenderMap"].push_back("BestExtrapolated");
			
			colors=m_DrawParams.GetStrValue(m_RendererName, "DrawMapColorsOnSurf", valid);
			m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"].clear();
			m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"].push_back("ErrorColors");
			mapFlags = m_DrawParams.GetStrValues(m_RendererName, "MapFlags", valid);			
		}
	}
	
	m_pPipeline->Draw(m_DrawParams);
	
	if (takingShot)
	{
		m_DrawParams.m_ParamModules[m_RendererName]["RenderMap"] = restoreVal;
		m_pPipeline->WriteScreenshot(CaptureImage(), screenshotName);
		if (screenshotName=="FinalColors")
		{
			m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"].clear();
			m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"].push_back(colors);
			m_DrawParams.m_ParamModules[m_RendererName]["MapFlags"] = mapFlags;
		}
		else if (screenshotName=="BestExtrapolated")
		{
			m_DrawParams.m_ParamModules[m_RendererName]["MaxNumMapColors"].clear();
			m_DrawParams.m_ParamModules[m_RendererName]["MaxNumMapColors"].push_back(maxNumPnts);
		}
		else if (screenshotName=="ErrorColors")
		{
			m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"].clear();
			m_DrawParams.m_ParamModules[m_RendererName]["DrawMapColorsOnSurf"].push_back(colors);
			m_DrawParams.m_ParamModules[m_RendererName]["MapFlags"] = mapFlags;
		}
		
		Update();	
	}
}

/**
 * Loads standard lights at some fixed positions
 */
void AnalysisWindow::SetLights(R3Mesh * mesh)
{
	 glEnable (GL_LIGHT0);
	 glEnable (GL_LIGHT1);	
//	static GLfloat specular[] = {1.0, 1.0, 1.0}; 
	static GLfloat specular[] = {0, 0, 0}; 	
	static GLfloat ambient[] = {0, 0, 0}; 	
	static GLfloat diffuse[] = {.8, .8, .8};
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
	
	// mesh center, radius
	R3Point cntr = mesh->Centroid();
	double rad = mesh->AverageRadius() * 20;	

	assert(m_MeshToLight1.find(mesh)!=m_MeshToLight1.end());
	assert(m_MeshToLight2.find(mesh)!=m_MeshToLight2.end());
	double theta1 = m_MeshToLight1[mesh][0];
	double phi1 = m_MeshToLight1[mesh][1];	
	double theta2 = m_MeshToLight2[mesh][0];
	double phi2 = m_MeshToLight2[mesh][1];	
	
	R3Vector pSphere1(cos(theta1)*sin(phi1), sin(theta1)*sin(phi1), cos(theta1));
	R3Vector pSphere2(cos(theta2)*sin(phi2), sin(theta2)*sin(phi2), cos(theta2));	
	
	R3Point light0Pos = cntr + rad * pSphere1;
	R3Point light1Pos = cntr + rad * pSphere2;	
	// position on a sphere around mesh
	
	// Set lights
//	static GLfloat light0_position[] = { 3.0, 4.0, 5.0, 0.0 };
	static GLfloat light0_position[3];
	light0_position[0] = light0Pos.X();
	light0_position[1] = light0Pos.Y();
	light0_position[2] = light0Pos.Z();	
	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
//	static GLfloat light1_position[] = { -3.0, -2.0, -3.0, 0.0 };
	static GLfloat light1_position[3];
	light1_position[0] = light1Pos.X();
	light1_position[1] = light1Pos.Y();
	light1_position[2] = light1Pos.Z();	
	glLightfv(GL_LIGHT1, GL_POSITION, light1_position);	
}

/**
 * Initializes openGL buffer & standard settings
 */
void AnalysisWindow::InitializeGL()
{
	// Initialize background color 
	glClearColor(200.0/255.0, 200.0/255.0, 200.0/255.0, 1.0);
	
	// Initialize lights
	static GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	static GLfloat light0_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
	glEnable(GL_LIGHT0);
	static GLfloat light1_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
	glEnable(GL_LIGHT1);
	glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHTING); 
	
	// Initialize graphics modes  
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);	
	
	// initialize alpha-blending
	glEnable (GL_BLEND); 
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	
}

/**
 * Captures image from the openGL buffer
 */
R2Image AnalysisWindow::CaptureImage()
{
    R2Image image(Width(), Height(), 3);
    image.Capture();
	return image;
}

/**
 * Width
 */
int AnalysisWindow::Width()
{
	return 800;
//	return 1600;
}

/**
 * Height
 */
int AnalysisWindow::Height()
{
	return 600;
//	return 1200;
}

/**
 * sets a mesh in collection (for manipulation)
 * @param collectionID	id of mesh collection (mesh + flattened mesh of a human)
 * @param meshID		id of mesh in collection, type of mesh (e.g. flattened mesh)
 */
void AnalysisWindow::SetCurrentMesh(int collectionID, int meshID)
{
	m_iCurrentMesh = collectionID;
	m_iCurrentSubMesh = meshID;
}

/**
 * Set mesh type
 * @param meshID	id of mesh in collection, type of mesh (e.g. flattened mesh)
 */
void AnalysisWindow::SetCurrentMeshInCollection(int meshID)
{
	m_iCurrentSubMesh = meshID;
}

/**
 * Set collection
 * @param collectionID	id of mesh collection (mesh + flattened mesh of a human)
 */
void AnalysisWindow::SetCurrentCollection(int collectionID)
{
	m_iCurrentMesh = collectionID;
}

/**
 * Print current camera settings: originXYZ towardsXYZ upXYZ FOV
 */
void AnalysisWindow::PrintCamera()
{
	// Print camera
	if (m_iCurrentMesh>=0 && m_iCurrentMesh<(int)m_ViewerCollections.size()
		&& m_iCurrentSubMesh>=0 && m_iCurrentSubMesh<(int)m_ViewerCollections.size())
	{
		const R3Camera& camera = m_ViewerCollections[m_iCurrentMesh][m_iCurrentSubMesh]->Camera();
		printf("#camera  %g %g %g  %g %g %g  %g %g %g  %g \n",
		   camera.Origin().X(), camera.Origin().Y(), camera.Origin().Z(),
		   camera.Towards().X(), camera.Towards().Y(), camera.Towards().Z(),
		   camera.Up().X(), camera.Up().Y(), camera.Up().Z(),
		   camera.YFOV());
	}
}

/**
 * Translate mesh 
 * @param collectionID
 * @param meshID			mesh identification (maps to viewer)
 * @param x		
 * @param y					coordinates of a mouse (screen space)
 * @param dx		
 * @param dy				displacement of a mouse (screen space)
 */
void AnalysisWindow::TranslateMesh(int collectionID, int meshID,  int x, int y, int dx, int dy)
{
	if (collectionID>=0 && collectionID<(int)m_ViewerCollections.size()
		&& meshID>=0 && meshID<(int)m_ViewerCollections[collectionID].size())
	{
		R3Point origin = m_MeshCollections[collectionID][meshID]->BBox().Centroid();
		m_ViewerCollections[collectionID][meshID]->TranslateWorld(1.0, origin, x, y, dx, dy);
	}
}

/**
 * Rotate mesh 
 * @param collectionID
 * @param meshID			mesh identification (maps to viewer)
 * @param x		
 * @param y					coordinates of a mouse (screen space)
 * @param dx		
 * @param dy				displacement of a mouse (screen space)
 */
void AnalysisWindow::RotateMesh(int collectionID, int meshID, int x, int y, int dx, int dy)
{
	if (collectionID>=0 && collectionID<(int)m_ViewerCollections.size()
		&& meshID>=0 && meshID<(int)m_ViewerCollections[collectionID].size())
	{
		R3Point origin = m_MeshCollections[collectionID][meshID]->BBox().Centroid();
		m_ViewerCollections[collectionID][meshID]->RotateWorld(1.0, origin, x, y, dx, dy);
	}
}

/**
 * Scale mesh 
 * @param collectionID
 * @param meshID			mesh identification (maps to viewer)
 * @param x		
 * @param y					coordinates of a mouse (screen space)
 * @param dx		
 * @param dy				displacement of a mouse (screen space)
 */
void AnalysisWindow::ScaleMesh(int collectionID, int meshID, int x, int y, int dx, int dy)
{
	if (collectionID>=0 && collectionID<(int)m_ViewerCollections.size()
		&& meshID>=0 && meshID<(int)m_ViewerCollections[collectionID].size())
	{
		R3Point origin = m_MeshCollections[collectionID][meshID]->BBox().Centroid();
		m_ViewerCollections[collectionID][meshID]->ScaleWorld(1.0, origin, x, y, dx, dy);
	}		
}


/**
 * Scale mesh 
 * @param collectionID
 * @param meshID			mesh identification (maps to viewer)
 * @param scale			scale factor
 */void AnalysisWindow::ScaleMesh(int collectionID, int meshID, 
								  double scale)
{
	if (collectionID>=0 && collectionID<(int)m_ViewerCollections.size()
		&& meshID>=0 && meshID<(int)m_ViewerCollections[collectionID].size())
	{
		R3Point origin = m_MeshCollections[collectionID][meshID]->BBox().Centroid();
		m_ViewerCollections[collectionID][meshID]->ScaleWorld(origin, scale);
	}			
}



/**
 * Translate mesh (currently set for manipulation) 
 * @param x		
 * @param y					coordinates of a mouse (screen space)
 * @param dx		
 * @param dy				displacement of a mouse (screen space)
 */
void AnalysisWindow::TranslateCurrentMesh(int x, int y, int dx, int dy)
{
	TranslateMesh(m_iCurrentMesh, m_iCurrentSubMesh, x, y, dx, dy);
}

/**
 * Rotate mesh (currently set for manipulation) 
 * @param x		
 * @param y					coordinates of a mouse (screen space)
 * @param dx		
 * @param dy				displacement of a mouse (screen space)
 */
void AnalysisWindow::RotateCurrentMesh(int x, int y, int dx, int dy)
{
	RotateMesh(m_iCurrentMesh, m_iCurrentSubMesh, x, y, dx, dy);
}

void AnalysisWindow::RotateLights(int lightID, int x, int y, int dx, int dy)
{
	if (m_iCurrentMesh>=0 && m_iCurrentMesh<(int)m_ViewerCollections.size()
		&& m_iCurrentSubMesh>=0 && m_iCurrentSubMesh<(int)m_ViewerCollections[m_iCurrentMesh].size())
	{
		double * angles=NULL;
		if (lightID==0)
			angles = m_MeshToLight1[m_MeshCollections[m_iCurrentMesh][m_iCurrentSubMesh]];
		else
			angles = m_MeshToLight2[m_MeshCollections[m_iCurrentMesh][m_iCurrentSubMesh]];
		assert(angles!=NULL);
			
		angles[0] += ((double)dx / (double)Width()) * 3.14;
		angles[1] += ((double)dy / (double)Height()) * 3.14;		
	}
	
}


/**
 * Scale mesh (currently set for manipulation) 
 * @param x		
 * @param y					coordinates of a mouse (screen space)
 * @param dx		
 * @param dy				displacement of a mouse (screen space)
 */
void AnalysisWindow::ScaleCurrentMesh(int x, int y, int dx, int dy)
{
	ScaleMesh(m_iCurrentMesh, m_iCurrentSubMesh, x, y, dx, dy);
}

/**
 * Add mesh (creates viewer, allows to select mesh)
 * @param mesh
 * @return collectionID of the mesh (meshID=0 in the collection)
 */
int AnalysisWindow::AddMesh(R3Mesh * mesh)
{
	double pi = 3.14159265;
	m_MeshToLight1[mesh] = new double[2];
	m_MeshToLight2[mesh] = new double[2];	
	m_MeshToLight1[mesh][0] = 0;
	m_MeshToLight1[mesh][1] = pi / 2.;	
	m_MeshToLight2[mesh][0] = pi;
	m_MeshToLight2[mesh][1] = -pi / 2.;	
	
	m_MeshCollections.push_back(std::vector<R3Mesh*>());
	m_MeshCollections[m_MeshCollections.size()-1].push_back(mesh);
	
	m_ViewerCollections.push_back(std::vector<R3Viewer *>());
	m_ViewerCollections[m_ViewerCollections.size()-1].push_back(CreateBirdsEyeViewer(mesh));
	PositionInSubWindow(m_ViewerCollections.size()-1);
	return m_ViewerCollections.size()-1;
}

/**
 * Add mesh (creates viewer, allows to select mesh) with specific collection id
 * @param mesh
 * @param collectionID	where to add mesh. If invalid - adds as mesh 0, with next collection
 * @return collectionID of the mesh (meshID=0 in the collection)
 */
int AnalysisWindow::AddMesh(R3Mesh * mesh, int collectionID)
{
	if (collectionID<0 || collectionID>=(int)m_ViewerCollections.size())
		return AddMesh(mesh);
	
	m_MeshCollections[collectionID].push_back(mesh);
	m_ViewerCollections[collectionID].push_back(CreateBirdsEyeViewer(mesh));
	
	PositionInSubWindow(collectionID);
	
	return collectionID;
}

void AnalysisWindow::PositionInSubWindow(int collectionID)
{
	std::vector<int> rowsCols=m_DrawParams.GetIntValues(m_RendererName, "RowsCols", valid);
	int numSurf = m_DrawParams.GetIntValue(m_RendererName, "NumSurfaces", valid);
	if (numSurf==2 && valid)
	{
		if (collectionID==0)
		{
			ScaleMesh(collectionID, 0, 0.6);
			TranslateMesh(collectionID, 0,  0, 0, (int)(-Width()/2.5), Height()/5);
		}
		else if (collectionID==1)
		{
			ScaleMesh(collectionID, 0, 0.75);			
			TranslateMesh(collectionID, 0,  0, 0, (int)(Width()/2.5), -Height()/5);	
		}
	}
	else if (rowsCols.size()==2)
	{
		int rows = rowsCols[1];
		int cols = rowsCols[0];	//NOTE: parameters switched
		double scaleFactor = vkMin(1. / rows, 1. / cols);
		ScaleMesh(collectionID, 0, scaleFactor);
		int row = collectionID / cols;
		int col = collectionID % cols;
		double offsetX = -(double)Width() + (double)Width()/rows;
		double offsetY = ((double)Height() - (double)Height()/cols)*.6;
		int transX = (int)(offsetX+1.5*(double)Width()*(double)row/(double)rows);
		int transY = (int)(offsetY-1.5*(double)Height()*(double)col/(double)cols);
		TranslateMesh(collectionID, 0,  0, 0, transX, transY);
	}
}

/**
 * Set camera for a mesh
 * @param collectionID
 * @param meshID			mesh identification
 * @param ox
 * @param oy
 * @param oz				origin
 * @param tx
 * @param ty
 * @param tz				towards
 * @param ux
 * @param uy
 * @param uz				up
 */
void AnalysisWindow::SetMeshCamera(int collectionID, int meshID, 
								   double ox, double oy, double oz, 
								   double tx, double ty, double tz,
								   double ux, double uy, double uz)
{
	assert(collectionID>=0 && collectionID<(int)m_ViewerCollections.size()
		   && meshID>=0 && meshID < (int)m_ViewerCollections[collectionID].size() );
	
	m_ViewerCollections[collectionID][meshID]->ResetCamera(R3Point(ox, oy, oz), 
														   R3Vector(tx, ty, tz),
														   R3Vector(ux, uy, uz));
	PositionInSubWindow(collectionID);
}

/**
 * Creates a bird-eye viewer for a mesh
 * @param mesh
 * @return the new viewer
 */
R3Viewer * AnalysisWindow::CreateBirdsEyeViewer(R3Mesh *mesh)
{
	R3Vector initial_camera_towards(-0.57735, -0.57735, -0.57735);
	R3Vector initial_camera_up(-0.57735, 0.57735, 0.5773);
	
    // Setup camera view looking down the Z axis
    R3Box bbox = mesh->BBox();
    assert(!bbox.IsEmpty());
    RNLength r = bbox.DiagonalRadius();
    assert((r > 0.0) && RNIsFinite(r));
    R3Point initial_camera_origin = bbox.Centroid() - initial_camera_towards * (2.5 * r);
    R3Camera camera(initial_camera_origin, initial_camera_towards, initial_camera_up, 0.4, 0.4, 0.01 * r, 100.0 * r);
    R2Viewport viewport(0, 0, Width(), Height());
    return new R3Viewer(camera, viewport);
}

/**
 * Load camera for a specified mesh
 */
void AnalysisWindow::LoadCamera(int collectionID, int meshID)
{
	m_iCurrentCameraCollectionID = collectionID;
	m_iCurrentCameraMeshID = meshID;
	assert(collectionID>=0 && collectionID<(int)m_ViewerCollections.size()
		   && meshID>=0 && meshID < (int)m_ViewerCollections[collectionID].size() );

	m_ViewerCollections[collectionID][meshID]->Camera().Load();
	//PrintCamera();
}

/**
 * Get point's position in current coordinate system
 */
R3Point AnalysisWindow::PositionInMyCoords(int collectionID, int meshID, const R3Point & pnt)
{
	assert(collectionID>=0 && collectionID<(int)m_ViewerCollections.size()
		   && meshID>=0 && meshID < (int)m_ViewerCollections[collectionID].size() );
	
	R4Matrix Minv = m_ViewerCollections[m_iCurrentCameraCollectionID][m_iCurrentCameraMeshID]->Camera().CoordSystem().Matrix();
	R4Matrix Mcurr = m_ViewerCollections[collectionID][meshID]->Camera().CoordSystem().InverseMatrix();
	return Minv * (Mcurr * pnt);
}

/**
 * Maps a value to some visible color: blue(small) -> light-blue -> green -> yellow -> red(high)
 * @param value
 * @param minValue
 * @param maxValue		range
 * @param r
 * @param g
 * @param b				output color
 */
void AnalysisWindow::MapScalarToColor(double value, double minValue, double maxValue, 
									  double & r, double & g, double & b)
{
	double scaledVal = (value-minValue) / (maxValue-minValue);	// scale to [0, 1]
	// linear greyscale
//	r = scaledVal;
//	g = scaledVal;
//	b = scaledVal;	
	
	// linear blue -> red
	int nPrev = (int)floor(scaledVal * ((int)m_MulticolorSequence.size()-1));
	int nNext = (int)ceil(scaledVal * ((int)m_MulticolorSequence.size()-1));
	double t = 1;
	if (nNext!=nPrev)
		t = ((double)nNext - (scaledVal * (double)m_MulticolorSequence.size())) / (nNext-nPrev);
	
	r = m_MulticolorSequence[nPrev].R() * t + m_MulticolorSequence[nNext].R() * (1-t);
	g = m_MulticolorSequence[nPrev].G() * t + m_MulticolorSequence[nNext].G() * (1-t);
	b = m_MulticolorSequence[nPrev].B() * t + m_MulticolorSequence[nNext].B() * (1-t);	
}

/**
 * maps integer to a color (colors very different from each other, but fewer than MapIntToColor)
 */
void AnalysisWindow::MapIntToRadicalColor(int i, double & r, double & g, double & b)
{
	switch(i)
	{
		case 0:
			r=0; g=1; b=0;
			break;
		case 1:
			r=0; g=0; b=1;
			break;
		case 2:
			r=1; g=1; b=0;
			break;
		case 3:
			r=0; g=1; b=1;
			break;
		case 4:
			r=1; g=0; b=1;
			break;
		case 5:
			r=0; g=0; b=0;
			break;
		case 6:
			r=0; g=.5; b=0;
			break;
		case 7:
			r=0; g=0; b=.5;
			break;
		case 8:
			r=.5; g=.5; b=0;
			break;
		case 9:
			r=0; g=.5; b=.5;
			break;
		default:
			r=0; g=0; b=0;
	}
}

/**
 * Maps integer to a color
 */
void AnalysisWindow::MapIntToColor(int i, double & r, double & g, double & b)
{
	r = (31 + 48 * ((3*i) % 5)) / 256.0;
	g = (31 + 32 * ((5*i) % 7)) / 256.0;
	b = (31 + 19 * ((7*i) % 11)) / 256.0;	
}

/**
 * Standard surface color - grayish, with some alpha
 */
void AnalysisWindow::StandardSurfaceColor(double & r, double & g, double & b, double & a)
{
	r = .8;
	g = .8;
	b = .8;
	a = GetAlpha();
}

/**
 * maps integer to a color
 */
void AnalysisWindow::MapIntToColor(int i, GLfloat material[4])
{
	double r, g, b;
	MapIntToColor(i, r, g, b);
	material[0] = r;
	material[1] = g;	
	material[2] = b;
	material[3] = 1;	
}

/**
 * Standard surface color - grayish with some alpha
 */
void AnalysisWindow::StandardSurfaceColor(GLfloat material[4])
{
	double r, g, b, a;
	StandardSurfaceColor(r, g, b, a);
	material[0] = r;
	material[1] = g;
	material[2] = b;
	material[3] = GetAlpha();	
}

/**
 * Flip rendering flag
 * @param renderer			name of a standard renderer
 * @param flagsCollName		collection of flags (flag in collection -> it is set to true, false o.w.) False flags can be prefixed with '#'
 * @param flagName			name of a flag to be flipped
 */
void AnalysisWindow::FlipFlag(const VKString & renderer, 
							  const VKString & flagsCollName,
							  const VKString & flagName)
{
	VKStringList flagList = m_DrawParams.GetStrValues(renderer, flagsCollName, valid);
	if (!valid)
		std::cout<<"[WARNING] AnalysisWindow.cpp failed to find flag to flip: "<<renderer.c_str()<<"::"<<flagsCollName.c_str()<<std::endl;

	for (int i=0; i<flagList.count(); i++)
	{
		if (flagList[i]==flagName)
		{
			m_DrawParams.m_ParamModules[renderer][flagsCollName].remove(i);
			return;
		}
	}
	
	m_DrawParams.m_ParamModules[renderer][flagsCollName].push_back(flagName);
}

/**
 * Flip rendering flag
 * @param renderer			name of a standard renderer
 * @param flagsCollName		collection of flags (flag in collection -> it is set to true, false o.w.) False flags can be prefixed with '#'
 * @param flagName			name of a flag 
 * @param setTo				value to set the flag to 
 */
void AnalysisWindow::SetFlag(const VKString & renderer, 
							 const VKString & flagsCollName,
							 const VKString & flagName, bool setTo)
{
	VKStringList flagList = m_DrawParams.GetStrValues(renderer, flagsCollName, valid);
	if (!valid)
		std::cout<<"[WARNING] AnalysisWindow.cpp failed to find flag to set: "<<renderer.c_str()<<"::"<<flagsCollName.c_str()<<std::endl;
	
	for (int i=0; i<flagList.count(); i++)
	{
		if (flagList[i]==flagName)
		{
			if (!setTo)
				m_DrawParams.m_ParamModules[renderer][flagsCollName].remove(i);
			return;
		}
	}
	
	if (setTo)
		m_DrawParams.m_ParamModules[renderer][flagsCollName].push_back(flagName);
	
}

/**
 * @return alpha value of a default surface
 */
double AnalysisWindow::GetAlpha()
{
	return m_DrawParams.GetDoubleValue(m_RendererName, "Alpha", valid);
}

/**
 * Checks value of a renderer flag for a default renderer
 * @param flagsCollName		collection is flag supposed to be in
 * @param flagName			flag name 
 * @param modelDetails		model can have some flags re-flipped, so this is link to a specific model additional flag flips
 * @param flipsCollection	flip collection specific to the model
 * @return true iff flag is set
 */
bool AnalysisWindow::IsFlagSetDefault(const VKString & flagsCollName, 
							   const VKString & flagName,
							   const VKString & modelDetails,
							   const VKString & flipsCollection)
{
	return IsFlagSet(m_RendererName, flagsCollName, flagName, modelDetails, flipsCollection);
}

/**
 * Checks value of a renderer flag for a default renderer
 * @param renderer			rendering module
 * @param flagsCollName		collection is flag supposed to be in
 * @param flagName			flag name 
 * @param modelDetails		model can have some flags re-flipped, so this is link to a specific model additional flag flips
 * @param flipsCollection	flip collection specific to the model
 * @return true iff flag is set
 */
bool AnalysisWindow::IsFlagSet(const VKString & renderer, 
							   const VKString & flagsCollName, 
							   const VKString & flagName,
							   const VKString & modelDetails,
							   const VKString & flipsCollection)
{
	bool setFlag = m_DrawParams.GetStrValues(renderer, flagsCollName, valid).contains(flagName);
	if (!valid)
		std::cout<<"[WARNING] AnalysisWindow.cpp Could not find flags: "<<renderer.c_str()<<"::"<<flagsCollName.c_str()<<std::endl;
	if (modelDetails!="none")
	{
		bool flipFlag = m_DrawParams.GetStrValues(modelDetails, flipsCollection, valid).contains(flagName);
		if (!valid)
			std::cout<<"[WARNING] AnalysisWindow.cpp Failed to check flag flip: "<<modelDetails.c_str()<<"::"<<flipsCollection.c_str()<<std::endl;
		if (flipFlag)
			setFlag = !setFlag;
	}
	
	return setFlag;
}

/**
 * called upon resize (children need to make sure ot implement / register calls to resize)
 */
void AnalysisWindow::Resize(int w, int h)
{
//	std::cout<<"[WARNING] AnalysisWindow.cpp Resize is not implemented!"<<std::endl;
}

void AnalysisWindow::PrintList(const VKString & listName, const VKStringList & stringList)
{
	std::cout<<listName.c_str()<<": "<<std::flush;
	for (int i=0; i<stringList.count(); i++)
		std::cout<<((VKStringList&)stringList)[i].c_str()<<" "<<std::flush;
	std::cout<<std::endl;
}

void AnalysisWindow::AddCorrespondence(int x, int y)
{
	VKString interactSmpSet=s_InteractiveSampleSetName;
	bool inWindow;
	double xWindow, yWindow;
	VKString windowContent = InSupplementalWindow(x, y, &inWindow, &xWindow, &yWindow);
	if (windowContent=="Textures" && inWindow)
	{
		// add corrs on texture
		VKString textureName = m_DrawParams.GetStrValue("RendererDefault", "RenderTexture", valid);
		if (textureName!="none")
		{
			SurfaceSampleSet *smpSet = Surface2DPlane::m_pInfinitePlanePseudosurface->GetSampleSet(textureName+interactSmpSet);
			if (smpSet==NULL)
			{
				smpSet = new SurfaceSampleSet();
				Surface2DPlane::m_pInfinitePlanePseudosurface->AddSampleSet(textureName+interactSmpSet, smpSet);
			}
			SurfaceSample samp(-1, xWindow, yWindow, 0, Surface2DPlane::m_pInfinitePlanePseudomesh);
			smpSet->AddSample(samp);

			std::cout<<"Adding Interactive Sample to Texture "<<textureName.c_str()<<": ";
			std::cout<<"["<<xWindow<<", "<<yWindow<<"] (N="<<smpSet->NumSamples()<<")"<<std::endl;;
		}
		else
		{
			std::cout<<"[WARNING] Select Texture (press ` for help)"<<std::endl;
		}
	}
	else
	{
		// add corrs on surface
		VKString surfaceName = m_DrawParams.GetStrValue("RendererDefault", "RenderSurface", valid);
		if (surfaceName!="none")
		{
			assert(PipelineGeneral::s_pGeneralPipeline!=NULL);
			SampledSurface * surf = PipelineGeneral::s_pGeneralPipeline->GetSurfaceByName(surfaceName);
			assert(surf!=NULL);
			SurfaceSampleSet * smpSet = surf->GetSampleSet(interactSmpSet);
			if (smpSet==NULL)
			{
				smpSet = new SurfaceSampleSet();
				surf->AddSampleSet(interactSmpSet, smpSet);
			}
			
			R3Ray screenRay = m_ViewerCollections[m_iCurrentMesh][m_iCurrentSubMesh]->WorldRay(x, y);
			double radius = surf->GetStandardRadius(&m_DrawParams);
			double t = -1;
			int selectedVertex = -1;
			for (int i=0; i<surf->GetMesh()->NVertices(); i++)
			{
				R3Point pos = surf->GetDrawablePosition(SurfaceSample(i, surf->GetMesh()), &m_DrawParams);
				double candT = CheckIntersection(pos, screenRay, radius);
				if ((t<0 || t>candT) && candT>=0 )
				{
					selectedVertex = i;
					t = candT;
				}
			}
			
			if (selectedVertex!=-1)
			{
				smpSet->AddSample(SurfaceSample(selectedVertex, surf->GetMesh()));
				std::cout<<"Adding Interactive Sample to Surface "<<surfaceName.c_str()<<": ";
				std::cout<<selectedVertex<<" (N="<<smpSet->NumSamples()<<")"<<std::endl;;
			}
		}
		else
		{
			std::cout<<"[WARNING] Select Surface (press ` for help)"<<std::endl;
		}
	}
}

int AnalysisWindow::FindSample(int x, int y, const VKString & interactSmpSet,
							   VKString & surfaceName, VKString & textureName)
{
	surfaceName = "";
	textureName = "";
	assert(m_pPipeline->GetPipelineType()=="PipelineGeneral");
	bool inWindow;
	double xWindow, yWindow;
	VKString windowContent = InSupplementalWindow(x, y, &inWindow, &xWindow, &yWindow);
	if (windowContent=="Textures" && inWindow)
	{
		double radius = .02;
		// add corrs on texture
		textureName = m_DrawParams.GetStrValue("RendererDefault", "RenderTexture", valid);
		if (textureName!="none")
		{
			SurfaceSampleSet *smpSet = Surface2DPlane::m_pInfinitePlanePseudosurface->GetSampleSet(textureName+interactSmpSet);
			if (smpSet==NULL)
				return -1;
			double bestDist = 10000;
			int retVal = -1;
			R3Point p(xWindow, yWindow, 0);
			for (int i=0; i<smpSet->NumSamples(); i++)
			{
				R3Point sp = smpSet->GetSample(i).GetPosition();
				double dist = (sp-p).Length();
				if (dist < radius && dist < bestDist)
				{
					retVal = i;
					bestDist = dist;
				}
			}
			return retVal;
		}
		else
		{
			std::cout<<"[WARNING] Select Texture (press ` for help)"<<std::endl;
		}
	}
	else
	{
		double t = -1;
		int selectedSample = -1;
		VKStringList surfaceNames = ((PipelineGeneral*)m_pPipeline)->GetSurfaceNames();
		for (int surfID=0; surfID<surfaceNames.count(); surfID++)
		{
			SampledSurface * surface = PipelineGeneral::s_pGeneralPipeline->GetSurfaceByName(surfaceNames[surfID]);
			assert(surface!=NULL);
			assert(m_MeshCollections[surfID][0]==surface->GetMesh());
			R3Ray screenRay = m_ViewerCollections[surfID][0]->WorldRay(x, y);
			surface->LoadCamera();
			double radius = surface->GetStandardRadius(&m_DrawParams);
			SurfaceSampleSet * surfSet = surface->GetSampleSet(interactSmpSet);
			if (surfSet==NULL)
			{
				std::cout<<"[WARNING] Could not find set (FindSample): "<<interactSmpSet.c_str()<<std::endl;
				return - 1;
			}
			
			for (int i=0; i<surfSet->NumSamples(); i++)
			{
				double candT = CheckIntersection(surface->GetDrawablePosition(surfSet->GetSample(i), &m_DrawParams), 
												 screenRay, radius);
				if ((t<0||t>candT) && candT>=0)
				{
					selectedSample = i;
					t = candT;
					surfaceName = surfaceNames[surfID];
				}
			}
		}
		
		return selectedSample;
	}
	return -1;
}

/**
 * If ray is shot through [x, y] position on the screen, select sample or a vertex on a mesh
 * saved in m_DrawParams.. SelectedVertex or SelectedSample
 */
void AnalysisWindow::SelectSample(int x, int y)
{
	if (m_iCurrentMesh<0 || m_iCurrentMesh>=(int)m_ViewerCollections.size()
		|| m_iCurrentSubMesh<0 || m_iCurrentSubMesh>=(int)m_ViewerCollections[m_iCurrentMesh].size())
		return;
	
	R3Ray screenRay = m_ViewerCollections[m_iCurrentMesh][m_iCurrentSubMesh]->WorldRay(x, y);
	SampledSurface * surface = m_pPipeline->GetSurface(m_iCurrentMesh);
	if (surface==NULL)
		return;
	surface->LoadCamera();	
	double radius = surface->GetStandardRadius(&m_DrawParams);
	
	int selectedSample=-1;
	int selectedVertex=-1;

	double t = -1;
	
	VKString setName = m_DrawParams.GetStrValue(m_RendererName, "RenderSampleSet", valid);
	if (!valid)
	{
		if (m_DrawParams.GetStrValues(m_RendererName, "SurfFlags", valid).contains("DefaultSamples"))
			setName="default";
		else
			setName="none";
	}
	if (setName!="none")	// rendering samples
	{
		SurfaceSampleSet * set = surface->GetSampleSet(setName);;
		if (set==NULL)
		{
			std::cout<<"[WARNING] Could not find set (SelectSample)"<<setName.c_str()<<std::endl;
			return;
		}
		
		for (int i=0; i<set->NumSamples(); i++)
		{
			double candT = CheckIntersection(surface->GetDrawablePosition(set->GetSample(i), &m_DrawParams), 
											 screenRay, radius);
			if ((t<0||t>candT) && candT>=0)
			{
				selectedSample = i;
				t = candT;
			}
		}
		
		if (selectedSample!=-1)
			selectedVertex = set->GetSample(selectedSample).NearestVertex();
	}
	else	// go over all vertices
	{
		for (int i=0; i<surface->GetMesh()->NVertices(); i++)
		{
			R3Point pos = surface->GetDrawablePosition(SurfaceSample(i, surface->GetMesh()), &m_DrawParams);
			double candT = CheckIntersection(pos, screenRay, radius);
			if ((t<0 || t>candT) && candT>=0 )
			{
				//selectedVertex = i;
				selectedVertex = i;
				t = candT;
			}
		}
		
	}
	
	std::cout<<"Selected: sample="<<selectedSample<<", vertex="<<selectedVertex<<std::endl;
	
	m_DrawParams.m_ParamModules[m_RendererName]["SelectedVertex"].clear();
	m_DrawParams.m_ParamModules[m_RendererName]["SelectedVertex"].push_back(m_iCurrentMesh);
	m_DrawParams.m_ParamModules[m_RendererName]["SelectedVertex"].push_back(selectedVertex);

	m_DrawParams.m_ParamModules[m_RendererName]["SelectedSample"].clear();
	m_DrawParams.m_ParamModules[m_RendererName]["SelectedSample"].push_back(m_iCurrentMesh);	
	m_DrawParams.m_ParamModules[m_RendererName]["SelectedSample"].push_back(selectedSample);
}


/**
 * sphere/ray intersection for selecting samples. 
 */
double AnalysisWindow::CheckIntersection(const R3Point & samplePos, 
										 const R3Ray & ray, double radius)
{	
	R3Sphere sphere(samplePos, radius);
	double x0 = ray.Start().X();
	double y0 = ray.Start().Y();
	double z0 = ray.Start().Z();	
	
	double xd = ray.Vector().X();
	double yd = ray.Vector().Y();
	double zd = ray.Vector().Z();	
	
	double xc = sphere.Center().X();
	double yc = sphere.Center().Y();
	double zc = sphere.Center().Z();	
	
	double Sr = sphere.Radius();
	
	double B = 2 * (xd*(x0-xc) + yd*(y0-yc) + zd*(z0-zc));
	double C = (x0-xc)*(x0-xc) + (y0-yc)*(y0-yc) + (z0-zc)*(z0-zc) - Sr*Sr;
	
	if (B*B-4*C < 0)
		return -1;
	else
		return (-B-sqrt(B*B-4*C)) / 2.;
}


void AnalysisWindow::DrawInAdditionalWindowBegin(bool top, bool leftSide)
{	
	VKString windowName;
	if (top)
		windowName="Top";
	else
		windowName="Bottom";
	
	if (leftSide)
		windowName+="Left";
	else
		windowName+="Right";
	windowName+="Window";
	
	VKStringList fracs = m_DrawParams.GetStrValues(m_RendererName, windowName, valid);
	assert(valid && fracs.count()==3);
	double width = fracs[1].toDouble();
	double height = fracs[2].toDouble();
	
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(-0.0, 1.0, -0.0, 1.0, 0., 10.);
	double transX, transY;
	transY = top ? (1-height) : 0;
	transX = leftSide ? 0 : (1-width);
	glTranslated(transX, transY, 0);
	glScaled(width, height, 1);
	
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	
	glDisable(GL_DEPTH_TEST);
	
	glDisable(GL_LIGHTING);
	glEnable (GL_BLEND); 
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	glBegin(GL_QUADS);
	glColor4d(1, 1, 1, .9);
	glVertex2d(0, 0);
	glVertex2d(1, 0);
	glVertex2d(1, 1);
	glVertex2d(0, 1);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glColor3d(0, .2, .5);
	glVertex2d(0, 0);
	glVertex2d(1, 0);
	glVertex2d(1, 1);
	glVertex2d(0, 1);
	glVertex2d(0, 0);	
	glEnd();
	
	glEnable(GL_DEPTH_TEST);

}

void AnalysisWindow::DrawInAdditionalWindowEnd()
{
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

VKString AnalysisWindow::InSupplementalWindow(int x, int y, bool * inWindow, double * xWindow, double * yWindow)
{	
	*inWindow = true;
	for (int i=0; i<(int)m_AvailiableSubWindows.size(); i++)
	{
		VKStringList windowDims = m_DrawParams.GetStrValues(m_RendererName, m_AvailiableSubWindows[i], valid);
		if (!valid)
		{
			*inWindow = false;
			return "";
		}
		assert(windowDims.count()==3);		
		
		double xRel = (double)x / (double)Width();
		double yRel = (double)y / (double)Height();	
		
		double fracX = windowDims[1].toDouble();
		double fracY = windowDims[2].toDouble();
		
		double minX=0, maxX=0, minY=0, maxY=0;
		
		if ((m_AvailiableSubWindows[i]=="TopRightWindow"))
		{
			minX = 1-fracX;
			maxX = 1;
			minY = 1-fracY;
			maxY = 1;
		}
		else if ((m_AvailiableSubWindows[i]=="TopLeftWindow"))
		{
			minX = 0;
			maxX = fracX;
			minY = 1-fracY;
			maxY = 1;
		}
		else if ((m_AvailiableSubWindows[i]=="BottomRightWindow"))
		{
			minX = 1-fracX;
			maxX = 1;
			minY = 0;
			maxY = fracX;
		}
		else if ((m_AvailiableSubWindows[i]=="BottomLeftWindow"))
		{
			minX = 0;
			maxX = fracX;
			minY = 0;
			maxY = fracX;
		}
		
		if (xRel >= minX && xRel<=maxX && yRel>=minY && yRel<=maxY)
		{
			*xWindow = (xRel - minX) / (maxX - minX);
			*yWindow = (yRel - minY) / (maxY - minY);			
			return windowDims[0];
		}
			
	}
	*inWindow = false;
	return "";
}

void AnalysisWindow::GenPipeline_UpdateRenderField(const VKString & name, int increase)
{
	assert(m_pPipeline->GetPipelineType()=="PipelineGeneral");
	VKStringList fields;
	if (name=="RenderSampleSet")
		fields = ((PipelineGeneral*)m_pPipeline)->GetSampleSetNames();
	else if (name=="RenderDistanceMetric")
		fields = ((PipelineGeneral*)m_pPipeline)->GetDistanceMetricNames();
	else if (name=="RenderFeatureValue")
		fields = ((PipelineGeneral*)m_pPipeline)->GetFeatureNames();
	else if (name=="RenderMap")
		fields = ((PipelineGeneral*)m_pPipeline)->GetMapNames();
	else if (name=="RenderOtherMap")
		fields = ((PipelineGeneral*)m_pPipeline)->GetMapNames();
	else if (name=="RenderConfidence")
		fields = ((PipelineGeneral*)m_pPipeline)->GetConfidenceNames();
	else if (name=="RenderSimilarity")
		fields = ((PipelineGeneral*)m_pPipeline)->GetSimilarityNames();
	else if (name=="RenderSurface")
		fields = ((PipelineGeneral*)m_pPipeline)->GetSurfaceNames();
	else if (name=="RenderTexture")
		fields = ((PipelineGeneral*)m_pPipeline)->GetTextureNames();
	else
	{
		std::cout<<"[ERROR] Cannot find renderable: "<<name.c_str()<<std::endl;
		assert(false);
	}
	
	int counterID=-1;
	if (m_mRenderFieldNameToSelectedID.find(name)==m_mRenderFieldNameToSelectedID.end())
		m_mRenderFieldNameToSelectedID[name] = -1;
	counterID = m_mRenderFieldNameToSelectedID[name];
	counterID += increase;
	
	if (counterID <= -2)
		counterID = fields.count()-1;
	else if (counterID >= fields.count())
		counterID = -1;

	VKString counterStrValue;
	if (counterID==-1)
		counterStrValue = "none";
	else
		counterStrValue = fields[counterID];
	
	std::cout<<name.c_str()<<": "<<counterStrValue.c_str()<<std::endl;
	m_DrawParams.m_ParamModules[m_RendererName][name].clear();
	m_DrawParams.m_ParamModules[m_RendererName][name].push_back(counterStrValue);
	m_mRenderFieldNameToSelectedID[name] = counterID;
}



