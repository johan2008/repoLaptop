#include "MapFlattening.h"

////////////////// GENERAL FLATTENING PARENT CLASS ////////////////
MapFlattening::MapFlattening(SampledSurface * s1)
: SurfaceMap(s1, Surface2DPlane::m_pInfinitePlanePseudosurface)
{
	m_fFlatMeshArea = -1;
	m_pFlatMesh = NULL;
	m_pFlatMeshKdVertexSearch = NULL;
	m_pPlanarTransform = NULL;
	m_aSearchNodes = NULL;
	m_bBilateralReflection = false;
}

MapFlattening::MapFlattening(SampledSurface * flattenMe, 
							 PlanarTransform * transform,
							 bool bilateralReflection,
							 R3Mesh * flatMesh,
							 R2Kdtree<FlatSearchNode*> * flatMeshTree)
: SurfaceMap(flattenMe, Surface2DPlane::m_pInfinitePlanePseudosurface)
{
	m_fFlatMeshArea = -1;
	m_pFlatMesh = flatMesh;
	m_pFlatMeshKdVertexSearch = flatMeshTree;
	m_aSearchNodes = NULL;
	m_pPlanarTransform = transform;
	m_bBilateralReflection = bilateralReflection;
}

void MapFlattening::SetPlanarTransformation(const VKString & coarseSampleSet,
									 std::vector<int> & sampleIDs,
									 std::vector<LinAlgComplex> & samples2DPositions,
									 const VKString & planarTransformClass)
{	
	if (planarTransformClass=="SameAsCurrent")
		assert(m_pPlanarTransform!=NULL);
	else if (planarTransformClass=="MobiusTransformation")
		m_pPlanarTransform = new MobiusTransformation();
	else
	{
		std::cout<<"[ERROR] Unknown planar transformation - "<<planarTransformClass.c_str()<<std::endl;
		assert(false);
	}
	assert(m_pPlanarTransform!=NULL);
	
	SurfaceSampleSet * sampleSet = m_pM1->GetSampleSet(coarseSampleSet);
	assert(sampleSet!=NULL);
	std::vector<LinAlgComplex> z;
	for (int i=0; i<(int)sampleIDs.size(); i++)
	{
		SurfaceSample samp = sampleSet->GetSample(sampleIDs[i]);
		z.push_back(Get2DPositionFromSample(samp));
	}
	
	assert(samples2DPositions.size()==z.size());
	m_pPlanarTransform->FindTransformation(z, samples2DPositions);
}

MapFlattening::~MapFlattening()
{
	if (m_aSearchNodes!=NULL)
		delete [] m_aSearchNodes;
	delete m_pPlanarTransform;
}

void MapFlattening::InteractiveUpdateFromCorrs(const VKString & sampleSetName,
											   const VKString & textureName)
{
	SurfaceSampleSet * sampleSetSurf = m_pM1->GetSampleSet(sampleSetName);
	if (sampleSetSurf==NULL)
		return;
	SurfaceSampleSet * sampleSetTx = m_pM2->GetSampleSet(textureName+sampleSetName);
	if (sampleSetTx==NULL)
		return;
	if (sampleSetSurf->NumSamples()!=sampleSetTx->NumSamples())
		return;

	std::vector<LinAlgComplex> z;
	for (int i=0; i<sampleSetSurf->NumSamples(); i++)
		z.push_back(Get2DPositionFromSample(sampleSetSurf->GetSample(i)));

	std::vector<LinAlgComplex> w;
	for (int i=0; i<sampleSetTx->NumSamples(); i++)	
	{
		R3Point fp = sampleSetTx->GetSample(i).GetPosition();
		w.push_back(LinAlgComplex(fp.X(), fp.Y()));
	}
	m_pPlanarTransform->FindTransformation(z, w);
//	for (int i=0; i<3; i++)
//		std::cout<<"Z["<<i<<"]="<<z[i].r<<" + "<<z[i].i<<"i"<<std::endl;
//	for (int i=0; i<3; i++)
//		std::cout<<"W["<<i<<"]="<<w[i].r<<" + "<<w[i].i<<"i"<<std::endl;
//	for (int i=0; i<3; i++)
//	{
//		LinAlgComplex tz = m_pPlanarTransform->Transform(z[i]);
//		std::cout<<"T(Z["<<i<<"])="<<tz.r<<" + "<<tz.i<<"i"<<std::endl;
//	}
	
	GetSurface(0)->ClearCachedSurfaceColors();
}

void MapFlattening::DrawFlatMesh(ParamParser & drawParams)
{
	assert(m_pFlatMesh!=NULL);
	SurfaceSampleSet * sampleSet = GetSurface(0)->GetRenderableSampleSet(&drawParams);

	double minX=0, maxX=0, minY=0, maxY=0;
	std::vector<LinAlgComplex> fPnts;
	std::vector<double> allXs;
	std::vector<double> allYs;
	for (int i=0; i<m_pFlatMesh->NVertices(); i++)
	{
		R3Point fpnt = m_pFlatMesh->VertexPosition(m_pFlatMesh->Vertex(i));
		allXs.push_back(fpnt.X());
		allYs.push_back(fpnt.Y());
		LinAlgComplex fp(fpnt.X(), fpnt.Y());
	}
	
	std::sort(allXs.begin(), allXs.end());
	std::sort(allYs.begin(), allYs.end());
	double throwPercent=.05;
	assert(allXs.size()==allYs.size());
	int lastElement = (int)((double)allXs.size() * (1.-throwPercent));
	int firstElement = (int)((double)allXs.size() * throwPercent);
	minX = allXs[firstElement];
	minY = allYs[firstElement];
	maxX = allXs[lastElement];
	maxY = allYs[lastElement];
	
	if (sampleSet!=NULL)
	{
		for (int i=0; i<sampleSet->NumSamples(); i++)
		{
			//R3Point fp = ForwardMap(sampleSet->GetSample(i)).GetPosition();
			LinAlgComplex fp = Get2DPositionFromSample(sampleSet->GetSample(i));
			//std::cout<<"X["<<i<<"] = "<<fp.r<<" , "<<fp.i<<std::endl;
			fPnts.push_back(fp);
		}
	
		//std::cout<<"X in ["<<minX<<", "<<maxX<<"]"<<std::endl;
		//std::cout<<"Y in ["<<minY<<", "<<maxY<<"]"<<std::endl;
		double width = maxX - minX;
		double height = maxY - minY;
		
		static GLfloat material[4];	
		material[3] = 1.;

		glEnable(GL_LIGHTING);
		double radius = .02;
			
		for (int i=0; i<(int)fPnts.size(); i++)
		{
			AnalysisWindow::s_pMainWindow->MapIntToColor(i, material);
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material); 
			double x = (fPnts[i].r-minX)/width;
			double y = (fPnts[i].i-minY)/height;
			//std::cout<<"sphere ["<<x<<", "<<y<<"]"<<std::endl;
			R3Sphere(R3Point(x, y, -.9), radius).Draw();
		}
	}
}

void MapFlattening::SetFlatMesh(R3Mesh * flatMesh, R2Kdtree<FlatSearchNode*> * kdTree)
{
	m_pFlatMesh = flatMesh;
	m_pFlatMeshKdVertexSearch = kdTree;
	assert(m_pFlatMesh!=NULL);
	if (m_pFlatMeshKdVertexSearch==NULL)
		m_pFlatMeshKdVertexSearch = SurfaceMidEdgeConf::CreateKdTree(m_pFlatMesh, &m_aSearchNodes);
	assert(m_pFlatMeshKdVertexSearch!=NULL);
}

SurfaceSample MapFlattening::ForwardMap(const SurfaceSample & s)
{
	assert(s.CheckMeshIsSame(m_pM1->GetMesh()));
	// map to flat
	LinAlgComplex pz = m_pPlanarTransform->Transform(Get2DPositionFromSample(s));
	return SurfaceSample(-1, pz.r, pz.i, 0, m_pM2->GetMesh());
}

SurfaceSample MapFlattening::InverseMap(const SurfaceSample & s)
{
	assert(s.CheckMeshIsSame(m_pM2->GetMesh()));
	R3Point pos = s.GetPosition();
	LinAlgComplex pz = m_pPlanarTransform->TransformInv(LinAlgComplex(pos.X(), pos.Y()));
	return GetSampleFrom2DPosition(pz);
}

LinAlgComplex MapFlattening::Get2DPositionFromSample(const SurfaceSample & s)
{
	assert(m_pFlatMesh!=NULL);
	assert(m_pFlatMeshKdVertexSearch!=NULL);
	
	SurfaceSample flatSamp = MapSurfaceToFlatMesh(s);
	R3Point flatPnt = flatSamp.GetPosition();
	R3Point p1 = m_pFlatMesh->VertexPosition(m_pFlatMesh->Vertex(0));
	R3Point p2 = m_pFlatMesh->VertexPosition(m_pFlatMesh->Vertex(1));
	int n = m_pFlatMesh->NVertices();
	p1 = m_pFlatMesh->VertexPosition(m_pFlatMesh->Vertex(n-1));
	p2 = m_pFlatMesh->VertexPosition(m_pFlatMesh->Vertex(n-2));
	assert(flatPnt.Z()==0);
	LinAlgComplex retVal = LinAlgComplex(flatPnt.X(), flatPnt.Y());
	if (m_bBilateralReflection)
		retVal.i *= -1.;
	return retVal;
}

SurfaceSample MapFlattening::GetSampleFrom2DPosition(const LinAlgComplex & s)
{
	assert(m_pFlatMesh!=NULL);
	assert(m_pFlatMeshKdVertexSearch!=NULL);
	
	LinAlgComplex pos2D = s;
	if (m_bBilateralReflection)
		pos2D.i *= -1.;
	// map s to flat mesh
	SurfaceSample fs = SurfaceMidEdgeConf::GetFlatSurfaceSampleKdTree(pos2D, m_pFlatMesh, 
																	  m_pFlatMeshKdVertexSearch);
	return MapFlatMeshToSurface(fs);
}

double MapFlattening::GetFlatMeshArea()
{
	if (m_fFlatMeshArea==-1)
	{
		m_fFlatMeshArea = 0;
		for (int i=0; i<m_pFlatMesh->NFaces(); i++)
			m_fFlatMeshArea += m_pFlatMesh->FaceArea(m_pFlatMesh->Face(i));
	}
	return m_fFlatMeshArea;
}

R3Mesh * MapFlattening::GetFlatMesh()
{
	return m_pFlatMesh;
}

R2Kdtree<FlatSearchNode*> * MapFlattening::GetKdTree()
{
	return m_pFlatMeshKdVertexSearch;
}


/////////////////////  Mid-Edge flattening /////////////////////////
MapFlatMidEdge::MapFlatMidEdge(SampledSurface * s1)
: MapFlattening(s1)
{
}

MapFlatMidEdge::MapFlatMidEdge(SampledSurface * flattenMe, 
							   PlanarTransform * transform, 
							   bool bilateralReflection,
							   R3Mesh * flatMesh, 
							   R2Kdtree<FlatSearchNode*> * flatMeshTree)
: MapFlattening(flattenMe, transform, bilateralReflection, flatMesh, flatMeshTree)
{
	if (m_pFlatMesh==NULL)
		GenerateFlatMesh();
	SetFlatMesh(m_pFlatMesh, m_pFlatMeshKdVertexSearch);
	assert(m_pFlatMesh!=NULL && m_pFlatMeshKdVertexSearch!=NULL);
}


MapFlatMidEdge::MapFlatMidEdge(SampledSurface * flattenMe, 
							   const VKString & findMobisuSampleSet, 
							   std::vector<int> & sampleIDs,
							   std::vector<LinAlgComplex> & samples2DPositions,
							   const VKString & planarTransformClass,
							   bool bilateralReflection,
							   R3Mesh * flatMesh, 
							   R2Kdtree<FlatSearchNode*> * flatMeshTree)
: MapFlattening(flattenMe, NULL, bilateralReflection, flatMesh, flatMeshTree)
{
	if (m_pFlatMesh==NULL)
		GenerateFlatMesh();
	SetFlatMesh(m_pFlatMesh, m_pFlatMeshKdVertexSearch);
	assert(m_pFlatMesh!=NULL && m_pFlatMeshKdVertexSearch!=NULL);
	SetPlanarTransformation(findMobisuSampleSet, sampleIDs, samples2DPositions, planarTransformClass);
}

MapFlatMidEdge::MapFlatMidEdge(SampledSurface * flattenMe, 
							   PlanarTransform * transform, 
							   SurfaceFeature * agd,
							   bool bilateralReflection,
							   R3Mesh * flatMesh, 
							   R2Kdtree<FlatSearchNode*> * flatMeshTree)
: MapFlattening(flattenMe, transform, bilateralReflection, flatMesh, flatMeshTree)
{
	if (m_pFlatMesh==NULL)
	{
		int faceToRip, v1, v2, v3;
		SurfaceMidEdgeConf::GetBestFlatteningParams(m_pM1->GetMesh(), faceToRip, v1, v2, v3, agd); 
		GenerateFlatMesh(faceToRip);
	}
	SetFlatMesh(m_pFlatMesh, m_pFlatMeshKdVertexSearch);
	assert(m_pFlatMesh!=NULL && m_pFlatMeshKdVertexSearch!=NULL);
}

MapFlatMidEdge::MapFlatMidEdge(SampledSurface * flattenMe, 
							   const VKString & findMobisuSampleSet, 
							   std::vector<int> & sampleIDs,
							   std::vector<LinAlgComplex> & samples2DPositions,
							   const VKString & planarTransformClass,
							   SurfaceFeature * agd, bool bilateralReflection, R3Mesh * flatMesh, 
							   R2Kdtree<FlatSearchNode*> * flatMeshTree)
: MapFlattening(flattenMe, NULL, bilateralReflection, flatMesh, flatMeshTree)
{
	if (m_pFlatMesh==NULL)
	{
		int faceToRip, v1, v2, v3;
		SurfaceMidEdgeConf::GetBestFlatteningParams(m_pM1->GetMesh(), faceToRip, v1, v2, v3, agd); 
		GenerateFlatMesh(faceToRip);
	}
	SetFlatMesh(m_pFlatMesh, m_pFlatMeshKdVertexSearch);
	assert(m_pFlatMesh!=NULL && m_pFlatMeshKdVertexSearch!=NULL);
	SetPlanarTransformation(findMobisuSampleSet, sampleIDs, samples2DPositions, planarTransformClass);
}

MapFlatMidEdge::~MapFlatMidEdge()
{
}

void MapFlatMidEdge::GenerateFlatMesh(int cutTriangle)
{	
	if (cutTriangle==-1)
		cutTriangle = 0;
	m_pFlatMesh = SurfaceMidEdgeConf::CreateFlattenedMidEdgeMesh(m_pM1->GetMesh(), cutTriangle);
	SurfaceMidEdgeConf::AddVerticesToFlatMesh(m_pM1->GetMesh(), m_pFlatMesh);
	m_pFlatMeshKdVertexSearch = SurfaceMidEdgeConf::CreateKdTree(m_pFlatMesh, &m_aSearchNodes);

}

SurfaceSample MapFlatMidEdge::MapSurfaceToFlatMesh(const SurfaceSample & s)
{
	return SurfaceMidEdgeConf::MapMeshToMidEdge(m_pM1->GetMesh(), m_pFlatMesh, s);
}

SurfaceSample MapFlatMidEdge::MapFlatMeshToSurface(const SurfaceSample & s)
{
	return SurfaceMidEdgeConf::MapMidEdgeToMesh(m_pM1->GetMesh(), m_pFlatMesh, s);	
}

VKString MapFlatMidEdge::GetSurfaceMapType()
{
	return "MapFlatMidEdge";
}

void MapFlatMidEdge::SaveMap(std::ofstream & textStream)
{
	textStream<<"Format MapFlatMidEdge\n";
	textStream<<"Reflection "<<m_bBilateralReflection<<"\n";
	m_pPlanarTransform->SaveTransformation(textStream);
}

void MapFlatMidEdge::LoadMap(std::ifstream & textStream)
{
	std::string tempStr;
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="Format");
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="MapFlatMidEdge");	
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="Reflection");
	textStream>>m_bBilateralReflection;
	
	m_pPlanarTransform = PlanarTransform::CreateTransform(textStream);
}


//////////////////// MapVia2DPlane ////////////////////
MapVia2DPlane::MapVia2DPlane(SampledSurface * s1, SampledSurface * s2)
: SurfaceMap(s1, s2)
{
}

MapVia2DPlane::MapVia2DPlane(SampledSurface * surface1, SampledSurface * surface2, 
							 MapFlattening * flattening1, MapFlattening * flattening2)
: SurfaceMap(surface1, surface2)
{
	m_pFlattening1 = flattening1;
	m_pFlattening2 = flattening2;
}

MapVia2DPlane::MapVia2DPlane(SampledSurface * surface1, SampledSurface * surface2, 
							 const VKString & corrSetName, std::vector<int> & corrs,
							 const VKString & createFlatteningType,
							 const VKString & planarTransformClass,
							 bool bilateralReflection,
							 R3Mesh * flat1, R2Kdtree<FlatSearchNode*> * tree1,
							 R3Mesh * flat2, R2Kdtree<FlatSearchNode*> * tree2)
: SurfaceMap(surface1, surface2)
{
	m_CorrsSampleSet = corrSetName;
	m_Corrs = corrs;
	if (createFlatteningType=="MapFlatMidEdge")
	{
		SurfaceSampleSet * set1 = surface1->GetSampleSet(corrSetName);
		SurfaceSampleSet * set2 = surface2->GetSampleSet(corrSetName);
		assert(set1!=NULL && set2!=NULL);
		
		PlanarTransform * identity = NULL;
		if (planarTransformClass=="MobiusTransformation")
			identity = new MobiusTransformation();
		else
		{
			std::cout<<"[ERROR] Unknown Planar Transform (MapVia2DPlane): ";
			std::cout<<planarTransformClass.c_str()<<std::endl;
			assert(false);
		}
		
		m_pFlattening1 = new MapFlatMidEdge (m_pM1, identity, bilateralReflection, flat1, tree1);
		std::vector<LinAlgComplex> w;
		std::vector<int> z;
		
		for (int i=0; i<(int)corrs.size(); i+=2)
		{
			SurfaceSample sampOnPlane = m_pFlattening1->ForwardMap(set1->GetSample(corrs[i+0]));
			R3Point flatPnt = sampOnPlane.GetPosition();
			w.push_back(LinAlgComplex(flatPnt.X(), flatPnt.Y()));
			z.push_back(corrs[i+1]);
		}

		m_pFlattening2 = new MapFlatMidEdge(m_pM2, corrSetName, z, w,
											planarTransformClass, false, flat2, tree2);
	}
	else
	{
		assert(false);
	}
}

MapVia2DPlane::~MapVia2DPlane()
{
}

VKString MapVia2DPlane::GetSurfaceMapType()
{
	return "MapVia2DPlane";
}

MapFlattening * MapVia2DPlane::GetFlattening(int surfaceID)
{
	if (surfaceID==0)
		return m_pFlattening1;
	else if (surfaceID==1)
		return m_pFlattening2;

	assert(false);
}

SurfaceSample MapVia2DPlane::ForwardMap(const SurfaceSample & s)
{
	return m_pFlattening2->InverseMap(m_pFlattening1->ForwardMap(s));
}

SurfaceSample MapVia2DPlane::InverseMap(const SurfaceSample & s)
{
	return m_pFlattening1->InverseMap(m_pFlattening2->ForwardMap(s));
}

void MapVia2DPlane::SetFlattenings(R3Mesh * flat1, R2Kdtree<FlatSearchNode*> * tree1,
								   R3Mesh * flat2, R2Kdtree<FlatSearchNode*> * tree2)
{
	assert(m_pFlattening1!=NULL && m_pFlattening2!=NULL);
	m_pFlattening1->SetFlatMesh(flat1, tree1);
	m_pFlattening2->SetFlatMesh(flat2, tree2);
}

void MapVia2DPlane::GetGeneratingSet(std::vector<int> ** corrs, SurfaceSampleSet ** sampleSet1, 
									 SurfaceSampleSet ** sampleSet2)
{
	*corrs = &m_Corrs;
	*sampleSet1 = GetSurface(0)->GetSampleSet(m_CorrsSampleSet);
	*sampleSet2 = GetSurface(1)->GetSampleSet(m_CorrsSampleSet);
}

void MapVia2DPlane::SaveMap(std::ofstream & textStream)
{
	textStream<<"Format MapVia2DPlane\n";
	
	// save generators (if any)
	textStream<<m_Corrs.size()<<"\n";
	if (m_Corrs.size()!=0)
	{
		for (int i=0; i<(int)m_Corrs.size(); i++)
			textStream<<m_Corrs[i]<<" ";
		textStream<<"\n";
		textStream<<m_CorrsSampleSet.c_str()<<"\n";
	}
	
	m_pFlattening1->SaveMap(textStream);
	m_pFlattening2->SaveMap(textStream);
}

void MapVia2DPlane::LoadMap(std::ifstream & textStream)
{
	std::string tempStr;
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="Format");
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="MapVia2DPlane");	
	// save generators (if any)
	int numCorrespondences=0;
	textStream>>numCorrespondences;
	
	if (numCorrespondences!=0)
	{
		for (int i=0; i<numCorrespondences; i++)
		{
			int tempVal;
			textStream>>tempVal;
			m_Corrs.push_back(tempVal);
		}
		m_CorrsSampleSet.readToken(textStream);
	}
	
	
	m_pFlattening1 = (MapFlattening*) SurfaceMap::CreateMap(textStream, m_pM1, NULL);
	m_pFlattening2 = (MapFlattening*) SurfaceMap::CreateMap(textStream, m_pM2, NULL);
}

void MapVia2DPlane::Draw(AnalysisWindow * window, ParamParser * params,
						 const VKString & renderingParams, 
						 const VKString & surfaceName)
{
	DrawCorrespondenceForSelected(params, renderingParams, surfaceName, true);	
	DrawGenerators(true, true, 1., 2., window, params, renderingParams, surfaceName);
}

void MapVia2DPlane::DrawGenerators(bool drawLines, bool drawSpheres,
								   int lineWidth, int sphereScale,
								   AnalysisWindow * window, ParamParser * params,
								   const VKString & renderingParams,
								   const VKString & surfaceName)
{
	SurfaceSampleSet * S1 = GetSurface(0)->GetSampleSet(m_CorrsSampleSet);
	SurfaceSampleSet * S2 = GetSurface(1)->GetSampleSet(m_CorrsSampleSet);	
	if (S1==NULL || S2==NULL)
	{
		std::cout<<"[ERROR] MapConformal::DrawGenerators - cannot find sample set: ";
		std::cout<<m_CorrsSampleSet.c_str()<<std::endl;
		assert(S1!=NULL && S2!=NULL);
	}
	
	for (int i=0; i<(int)m_Corrs.size(); i+=2)
	{
		if (drawSpheres)
		{
			double radius1 = GetSurface(0)->GetStandardRadius(params, renderingParams, surfaceName) * sphereScale;
			double radius2 = GetSurface(1)->GetStandardRadius(params, renderingParams, surfaceName) * sphereScale;
			
			glEnable(GL_LIGHTING);	
			static GLfloat material[4];	
			double r1, g1, b1;
			window->MapIntToRadicalColor(i/2, r1, g1, b1);
			material[0] = r1;		material[1] = g1;		material[2] = b1;		material[3] = 1;
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material);
			
			GetSurface(0)->LoadCamera();
			R3Point p1 = GetSurface(0)->GetDrawablePosition(S1->GetSample(m_Corrs[i+0]), 
															params, renderingParams, surfaceName);				
			
			if (GetSurface(0)!=GetSurface(1))
			{
				R3Sphere(p1, radius1).Draw();
				
				GetSurface(1)->LoadCamera();				
				R3Point p2 = GetSurface(1)->GetDrawablePosition(S2->GetSample(m_Corrs[i+1]), 
																params, renderingParams, surfaceName);				
				R3Sphere(p2, radius2).Draw();
			}
			else
			{
				if (m_Corrs[i+0]==m_Corrs[i+1])	// stationary
					R3Sphere(p1, radius1*1.5).Draw();
				else
					R3Sphere(p1, radius1).Draw();
			}
		}
		
		if (drawLines)
		{
			GetSurface(0)->LoadCamera();	
			R3Point p1 = GetSurface(0)->GetDrawablePosition(S1->GetSample(m_Corrs[i+0]), 
															params, renderingParams, surfaceName);	
			R3Point p2 = GetSurface(1)->GetDrawablePosition(S2->GetSample(m_Corrs[i+1]), 
															params, renderingParams, surfaceName);				
			
			glLineWidth(lineWidth);		
			glDisable(GL_LIGHTING);	
			glColor3d(.8, 0, 0);
			
			glBegin(GL_LINES);
			glVertex3d(p1.X(), p1.Y(), p1.Z());	
			glVertex3d(p2.X(), p2.Y(), p2.Z());
			glEnd();			
		}
	}
}

