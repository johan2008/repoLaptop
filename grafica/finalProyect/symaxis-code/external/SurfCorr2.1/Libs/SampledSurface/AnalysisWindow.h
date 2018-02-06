#ifndef __ANALYSIS_WINDOW_H
#define __ANALYSIS_WINDOW_H

#include "ParamParser.h"
#include "SurfaceSample.h"
#include "gaps.h"

class AnalysisPipeline;

/**
 * UI Window. Should have different implementations (e.g. GLUT-based, QT-based, etc...)
 */
class AnalysisWindow
{
	public:
		enum TextSize
		{
			TEXT_SIZE_SMALL,
			TEXT_SIZE_MEDIUM,
			TEXT_SIZE_LARGE,
		};
		AnalysisWindow(ParamParser & drawParams, AnalysisPipeline * pipeline, 
					   const VKString & rendererName = "RendererDefault");
		virtual ~AnalysisWindow(){}
		virtual void RenderText(R3Point pnt, const VKString & text, 
								double r, double g, double b, TextSize size=TEXT_SIZE_MEDIUM) = 0;
		virtual void Terminate() = 0;
		virtual void Update() = 0;
		virtual void Draw();
		virtual void InitializeGL();
	
		virtual void FlipFlag(const VKString & renderer, 
							  const VKString & flagsCollName,
							  const VKString & flagName);

		virtual void SetFlag(const VKString & renderer, 
							 const VKString & flagsCollName,
							 const VKString & flagName, bool setTo);	
	
		virtual bool IsFlagSet(const VKString & renderer, 
							   const VKString & flagsCollName, 
							   const VKString & flagName,
							   const VKString & modelDetails = "none",
							   const VKString & flipsCollection = "none");
		virtual bool IsFlagSetDefault(const VKString & flagsCollName, 
							   const VKString & flagName,
							   const VKString & modelDetails = "none",
							   const VKString & flipsCollection = "none");
	
		virtual void SetCurrentMesh(int collectionID, int meshID = 0);
		virtual void SetCurrentMeshInCollection(int meshID);
		virtual void SetCurrentCollection(int collectionID);	
		virtual void PrintCamera();
		virtual void TranslateCurrentMesh(int x, int y, int dx, int dy);
		virtual void RotateCurrentMesh(int x, int y, int dx, int dy);
		virtual void ScaleCurrentMesh(int x, int y, int dx, int dy);	
	
		virtual void ScaleMesh(int collectionID, int meshID, double scale);	
		
		virtual R2Image CaptureImage();	
		virtual void Resize(int w, int h);
		virtual int Width();
		virtual int Height();
		virtual int AddMesh(R3Mesh * mesh);
		virtual int AddMesh(R3Mesh * mesh, int collectionID);
		virtual void SetMeshCamera(int collectionID, int meshID, 
								   double ox, double oy, double oz, 
								   double tx, double ty, double tz,
								   double ux, double uy, double uz);
		virtual void PositionInSubWindow(int collectionID);
		
		virtual void LoadCamera(int collectionID, int meshID = 0);
		virtual R3Point PositionInMyCoords(int collectionID, int meshID, const R3Point & pnt);
		virtual void SetLights(R3Mesh * mesh);
		virtual R3Viewer * CreateBirdsEyeViewer(R3Mesh *mesh);
	
		virtual void MapScalarToColor(double value, double minValue, double maxValue, 
									  double & r, double & g, double & b);
		virtual void MapIntToColor(int i, double & r, double & g, double & b);
		virtual void StandardSurfaceColor(double & r, double & g, double & b, double & a);

		virtual void MapIntToColor(int i, GLfloat material[4]);
		virtual void StandardSurfaceColor(GLfloat material[4]);
	
		virtual void MapIntToRadicalColor(int i, double & r, double & g, double & b);
	
		virtual VKString InSupplementalWindow(int x, int y, bool * inWindow, double * xWindow, double * yWindow);
	
		virtual void AddCorrespondence(int x, int y);
		virtual int FindSample(int x, int y, const VKString & interactiveSet, 
							   VKString & surfaceName, VKString & textureName);
		virtual void SelectSample(int x, int y);
		virtual double CheckIntersection(const R3Point & samplePos, 
										 const R3Ray & screenRay, double radius);
		
		virtual double GetAlpha();

		virtual void DrawInAdditionalWindowBegin(bool top, bool leftSide);
		virtual void DrawInAdditionalWindowEnd();
		
		virtual void RotateLights(int lightID, int x, int y, int dx, int dy);
		
		virtual void PrintList(const VKString & listName, const VKStringList & stringList);
	
	
		virtual void GenPipeline_UpdateRenderField(const VKString & name, int increase);
	
		static const VKString s_InteractiveSampleSetName;	
		static AnalysisWindow * s_pMainWindow;
	
	protected:
		std::map<VKString, int> m_mRenderFieldNameToSelectedID;
	
		virtual void TranslateMesh(int collectionID, int meshID, int x, int y, int dx, int dy);
		virtual void RotateMesh(int collectionID, int meshID, int x, int y, int dx, int dy);
		virtual void ScaleMesh(int collectionID, int meshID, int x, int y, int dx, int dy);	
		
		bool valid;
		int m_iCurrentMesh;
		int m_iCurrentSubMesh;	//e.g. original vs flattened
		std::vector<std::vector<R3Mesh *> > m_MeshCollections;
		std::vector<std::vector<R3Viewer *> > m_ViewerCollections;
		AnalysisPipeline * m_pPipeline;
		ParamParser & m_DrawParams;
		VKString m_RendererName;
		std::vector<RNRgb> m_MulticolorSequence;
		
		int m_iCurrentCameraCollectionID;
		int m_iCurrentCameraMeshID;
	
		std::vector<VKString> m_AvailiableSubWindows;
		
		std::map<R3Mesh *, double * > m_MeshToLight1;
		std::map<R3Mesh *, double * > m_MeshToLight2;
};

#endif



