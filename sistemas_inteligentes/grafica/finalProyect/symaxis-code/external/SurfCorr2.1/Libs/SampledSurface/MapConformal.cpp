#include "MapConformal.h"
#include "MapCoarse.h"
#include "SurfaceDistance.h"
#include "AnalysisStats.h"

MapConformal::MapConformal(SurfaceMidEdgeConf * M1, SurfaceMidEdgeConf * M2, 
						   const std::vector<int> & correspondences,
						   const VKString & genSetName )
: SurfaceMap(M1, M2)
{
	m_fDrawOnlyGenR = -1;
	m_fDrawOnlyGenG = -1;
	m_fDrawOnlyGenB = -1;	
	
	//InitializeCache();
	m_pConfM1 = M1;
	m_pConfM2 = M2;
	
	// store generating correspondences
	m_GenSetName = genSetName;
	m_GeneratorCorrespondences = correspondences;
	UpdateGeneratingCorrsChanged();
}

MapConformal::MapConformal(SurfaceMidEdgeConf * M1, SurfaceMidEdgeConf * M2, 
						   const MobiusTransformation & m1, const MobiusTransformation & m2)
: SurfaceMap(M1, M2)
{	
	m_fDrawOnlyGenR = -1;
	m_fDrawOnlyGenG = -1;
	m_fDrawOnlyGenB = -1;	
	
	//InitializeCache();
	m_pConfM1 = M1;
	m_pConfM2 = M2;
	m_MapTransform1 = m1;
	m_MapTransform2 = m2;
}

void MapConformal::UpdateGeneratingCorrsChanged()
{
	if (IsCacheEnabled())
		ClearCache();
	ClearAllConfidences();
	ClearAllSimilarities();
	
	SurfaceSampleSet * S1 = GetSurface(0)->GetSampleSet(m_GenSetName);
	SurfaceSampleSet * S2 = GetSurface(1)->GetSampleSet(m_GenSetName);	
	assert(S1!=NULL);
	assert(S2!=NULL);
	
	// find transformations
	assert(m_GeneratorCorrespondences.size()>=6);
	//	std::cout<<"\nGetting Conf Coords: "<<std::endl;
	//	std::cout<<"\t"<<correspondences[0]<<" -> "<<correspondences[1]<<std::endl;
	//	std::cout<<"\t"<<correspondences[2]<<" -> "<<correspondences[3]<<std::endl;
	//	std::cout<<"\t"<<correspondences[4]<<" -> "<<correspondences[5]<<std::endl;	
	LinAlgComplex m1z1 = m_pConfM1->GetConfCoordOriginal(S1->GetSample(m_GeneratorCorrespondences[0]));
	LinAlgComplex m2z1 = m_pConfM2->GetConfCoordOriginal(S2->GetSample(m_GeneratorCorrespondences[1]));
	LinAlgComplex m1z2 = m_pConfM1->GetConfCoordOriginal(S1->GetSample(m_GeneratorCorrespondences[2]));
	LinAlgComplex m2z2 = m_pConfM2->GetConfCoordOriginal(S2->GetSample(m_GeneratorCorrespondences[3]));
	LinAlgComplex m1z3 = m_pConfM1->GetConfCoordOriginal(S1->GetSample(m_GeneratorCorrespondences[4]));	
	LinAlgComplex m2z3 = m_pConfM2->GetConfCoordOriginal(S2->GetSample(m_GeneratorCorrespondences[5]));
	//	std::cout<<"Getting Conf Coords - DONE!"<<std::endl;	
	
	if (m_GeneratorCorrespondences.size()==6)	// align both to canonical frame
	{
		LinAlgComplex y1(1, 0);
		LinAlgComplex y2(cos(2 * 3.14159265 / 3), sin(2 * 3.14159265 / 3));
		LinAlgComplex y3(cos(4 * 3.14159265 / 3), sin(4 * 3.14159265 / 3));
		
		m_MapTransform1.ExactFit(m1z1, m1z2, m1z3, y1, y2, y3);
		m_MapTransform2.ExactFit(m2z1, m2z2, m2z3, y1, y2, y3);		
	}
	else if (m_GeneratorCorrespondences.size()==8)	
		// align second in canonical frame (approx), then first to 2nd
	{
		LinAlgComplex m1z4 = m_pConfM1->GetConfCoordOriginal(S1->GetSample(m_GeneratorCorrespondences[6]));	
		LinAlgComplex m2z4 = m_pConfM2->GetConfCoordOriginal(S2->GetSample(m_GeneratorCorrespondences[7]));	
		
		std::vector<LinAlgComplex> z1;		
		std::vector<LinAlgComplex> z2;
		std::vector<LinAlgComplex> w1;		
		std::vector<LinAlgComplex> w2;
		z2.push_back(m2z1);		z2.push_back(m2z4);		z2.push_back(m2z2);		z2.push_back(m2z3);
		for (int i=0; i<4; i++)
		{
			double alpha = 3.14159265 / 4 + i*3.14159265/2;
			w2.push_back(LinAlgComplex(cos(alpha), sin(alpha)));
		}
		m_MapTransform2.LeastSquaresFit(z2, w2);
		m_pConfM2->Transform(m_MapTransform2);
		
		z1.push_back(m1z1);		z1.push_back(m1z2);		z1.push_back(m1z3);		z1.push_back(m1z4);		
		w1.push_back(m_pConfM2->GetConfCoordOriginal(S2->GetSample(m_GeneratorCorrespondences[1])));
		w1.push_back(m_pConfM2->GetConfCoordOriginal(S2->GetSample(m_GeneratorCorrespondences[3])));
		w1.push_back(m_pConfM2->GetConfCoordOriginal(S2->GetSample(m_GeneratorCorrespondences[5])));
		w1.push_back(m_pConfM2->GetConfCoordOriginal(S2->GetSample(m_GeneratorCorrespondences[7])));		
		
		m_MapTransform1.LeastSquaresFit(z1, w1);
	}
	else	// find a mobius to align first with the second
	{
		std::vector<LinAlgComplex> z1;
		std::vector<LinAlgComplex> w1;		
		
		for (int i=0; i<(int)m_GeneratorCorrespondences.size(); i+=2)
		{
			LinAlgComplex m1z = m_pConfM1->GetConfCoordOriginal(S1->GetSample(m_GeneratorCorrespondences[i+0]));
			LinAlgComplex m2w = m_pConfM1->GetConfCoordOriginal(S2->GetSample(m_GeneratorCorrespondences[i+1]));	
			z1.push_back(m1z);
			w1.push_back(m2w);
		}
		m_MapTransform1.LeastSquaresFit(z1, w1);
	}
}

SurfaceSample MapConformal::ForwardMap(const SurfaceSample & s)
{
	assert(!s.Invalid());
//	std::cout<<"ConformapMap::ForwardMap: "<<std::endl;
//	std::cout<<"\ts1=["<<s.VertID(0)<<", "<<s.VertID(1)<<", "<<s.VertID(2)<<"]"<<std::flush;
//	std::cout<<"\tweights=["<<s.B(0)<<", "<<s.B(1)<<", "<<s.B(2)<<"]"<<std::endl;		
	
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("MapConformal_ForwardMap");
	if (IsCacheEnabled())
	{
		bool valid;
		SurfaceSample retVal = GetCachedForw(s, valid);
		if (valid)
			return retVal;
	}

	SurfaceSample ret = MapSample(s, m_MapTransform1, m_MapTransform2, m_pConfM1, m_pConfM2);
//	std::cout<<"\ts2=["<<ret.VertID(0)<<", "<<ret.VertID(1)<<", "<<ret.VertID(2)<<"] "<<std::flush;
//	std::cout<<"\tweights=["<<ret.B(0)<<", "<<ret.B(1)<<", "<<ret.B(2)<<"]"<<std::endl;
		
	if (IsCacheEnabled() && !IsCacheLockedForSave())
		AddSampleToForwMap(s, ret);
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("MapConformal_ForwardMap");
	return ret;
}


SurfaceSample MapConformal::InverseMap(const SurfaceSample & s)
{
	assert(!s.Invalid());	
//	std::cout<<"ConformapMap::InverseMap: "<<std::endl;
//	std::cout<<"\ts1=["<<s.VertID(0)<<", "<<s.VertID(1)<<", "<<s.VertID(2)<<"]"<<std::flush;
//	std::cout<<"\tweights=["<<s.B(0)<<", "<<s.B(1)<<", "<<s.B(2)<<"]"<<std::endl;		
	
//	AnalysisStats::m_GlobalStats.m_Timing.startedProcess("MapConformal_InverseMap");	
	if (IsCacheEnabled())
	{
		bool valid;
		SurfaceSample retVal = GetCachedBack(s, valid);
		if (valid)
			return retVal;
	}
	
	SurfaceSample ret = MapSample(s, m_MapTransform1, m_MapTransform2, m_pConfM1, m_pConfM2);
	assert(!ret.Invalid());
	
//	std::cout<<"\ts2=["<<ret.VertID(0)<<", "<<ret.VertID(1)<<", "<<ret.VertID(2)<<"] "<<std::flush;
//	std::cout<<"\tweights=["<<ret.B(0)<<", "<<ret.B(1)<<", "<<ret.B(2)<<"]"<<std::endl;	
	
	if (IsCacheEnabled() && !IsCacheLockedForSave())
		AddSampleToBackMap(s, ret);
//	AnalysisStats::m_GlobalStats.m_Timing.finishedProcess("MapConformal_InverseMap");
	return ret;
}

SurfaceSample MapConformal::MapSample(const SurfaceSample & s, 
									  const MobiusTransformation & t1, 
									  const MobiusTransformation & t2,
									  SurfaceMidEdgeConf * surf1, 
									  SurfaceMidEdgeConf * surf2)
{	
	assert(!s.Invalid());	
//	std::cout<<"\n///////////// MAPPING CONFORMAL //////////"<<std::endl;
	
	bool forward = s.CheckMeshIsSame(surf1->GetMesh());
	if (!forward)
		assert(s.CheckMeshIsSame(surf2->GetMesh()));
	
//	std::cout<<"Mobius Transforms:"<<std::endl;
//	t1.Print();
//	std::cout<<std::endl;
//	t2.Print();
	
	if (t1!=surf1->GetCurrentTransform())
		surf1->Transform(t1);
	
	if (t2!=surf2->GetCurrentTransform())
		surf2->Transform(t2);
	
//	std::cout<<"\nForward="<<forward<<std::endl;

	LinAlgComplex confCoord = forward ? surf1->GetConfCoord(s) : surf2->GetConfCoord(s);
//	std::cout<<"ConfCoord="<<confCoord.r<<" + "<<confCoord.i<<" i"<<std::endl;
	
	SurfaceSample ret = forward ? surf2->GetSurfaceSample(confCoord) : surf1->GetSurfaceSample(confCoord);
//	std::cout<<"RetVertex = "<<ret.NearestVertex()<<std::endl;
	
	assert(!ret.Invalid());
	return ret;
}

void MapConformal::GetGeneratingCorrespondences(std::vector<int> ** correspondences,
												SurfaceSampleSet ** sampleSet1,
												SurfaceSampleSet ** sampleSet2)
{
	*correspondences = &m_GeneratorCorrespondences;
	*sampleSet1 = GetSurface(0)->GetSampleSet(m_GenSetName);
	*sampleSet2 = GetSurface(1)->GetSampleSet(m_GenSetName);
	assert(*sampleSet1!=NULL);
	assert(*sampleSet2!=NULL);	
}

VKString MapConformal::GetSurfaceMapType()
{
	return "MapConformal";
}

const MobiusTransformation & MapConformal::GetMobiusForSurface(int surfaceID) const
{
	if (surfaceID==0)
		return m_MapTransform1;
	else if (surfaceID==1)
		return m_MapTransform2;
	else
		assert(false);
}

void MapConformal::SetDrawOnlyGenerators(double r, double g, double b)
{
	m_fDrawOnlyGenR = r;
	m_fDrawOnlyGenG = g;
	m_fDrawOnlyGenB = b;	
}

void MapConformal::DrawGenerators(AnalysisWindow * window, ParamParser * params, 
								  const VKString & renderingParams, 
								  const VKString & surfaceName, double linewidth, 
								  double r, double g, double b, bool lines, bool pnts)
{
	SurfaceSampleSet * S1 = GetSurface(0)->GetSampleSet(m_GenSetName);
	SurfaceSampleSet * S2 = GetSurface(1)->GetSampleSet(m_GenSetName);	
	if (S1==NULL || S2==NULL)
	{
		std::cout<<"[ERROR] MapConformal::DrawGenerators - cannot find sample set: ";
		std::cout<<m_GenSetName.c_str()<<std::endl;
		assert(S1!=NULL && S2!=NULL);
	}
	
	for (int i=0; i<(int)m_GeneratorCorrespondences.size(); i+=2)
	{
		if (pnts)
		{
			double radius1 = GetSurface(0)->GetStandardRadius(params, renderingParams, surfaceName) * 3;
			double radius2 = GetSurface(1)->GetStandardRadius(params, renderingParams, surfaceName) * 3;

			glEnable(GL_LIGHTING);	
			static GLfloat material[4];	
			double r1, g1, b1;
			window->MapIntToRadicalColor(i/2, r1, g1, b1);
			material[0] = r1;		material[1] = g1;		material[2] = b1;		material[3] = 1;
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material);
			
			GetSurface(0)->LoadCamera();
			R3Point p1 = GetSurface(0)->GetDrawablePosition(S1->GetSample(m_GeneratorCorrespondences[i+0]), 
															params, renderingParams, surfaceName);				
			R3Sphere(p1, radius1).Draw();
			GetSurface(1)->LoadCamera();				
			R3Point p2 = GetSurface(1)->GetDrawablePosition(S2->GetSample(m_GeneratorCorrespondences[i+1]), 
															params, renderingParams, surfaceName);				
			R3Sphere(p2, radius2).Draw();
		}
		
		if (lines)
		{
			GetSurface(0)->LoadCamera();	
			R3Point p1 = GetSurface(0)->GetDrawablePosition(S1->GetSample(m_GeneratorCorrespondences[i+0]), 
													params, renderingParams, surfaceName);	
			R3Point p2 = GetSurface(1)->GetDrawablePosition(S2->GetSample(m_GeneratorCorrespondences[i+1]), 
													params, renderingParams, surfaceName);				
			
			glLineWidth(linewidth);		
			glDisable(GL_LIGHTING);	
			if (r==-1) r = .8;
			if (g==-1) g = 0.;
			if (b==-1) b = 0.;
			glColor3d(r, g, b);
			
			glBegin(GL_LINES);
			glVertex3d(p1.X(), p1.Y(), p1.Z());	
			glVertex3d(p2.X(), p2.Y(), p2.Z());
			glEnd();			
		}
	}
}

void MapConformal::Draw(AnalysisWindow * window, ParamParser * params,
						const VKString & renderingParams, const VKString & surfaceName)
{
	DrawCorrespondenceForSelected(params, renderingParams, surfaceName, true);	
	bool valid;
		// TODO: can do colors as well
	if (params->GetStrValue(renderingParams, "ConfSurf", valid)=="RenderOriginal")
	{
		bool renderCoarse = params->GetStrValues(renderingParams, "RenderCorrCoarse", valid).contains(GetSurfaceMapType());

		if (renderCoarse && GetCurrentCoarseMap()!=NULL && m_fDrawOnlyGenR==-1)
			GetCurrentCoarseMap()->Draw(window, params, renderingParams, surfaceName);
		
		DrawGenerators(window, params, renderingParams, surfaceName,
					   (m_fDrawOnlyGenB==1.) ? 1. : 2., m_fDrawOnlyGenR, m_fDrawOnlyGenG, m_fDrawOnlyGenB);
			// rendering for paper figure: similarity
//		DrawGenerators(window, params, renderingParams, surfaceName,
//					   (m_fDrawOnlyGenB==1.) ? 2. : 6., m_fDrawOnlyGenR, m_fDrawOnlyGenG, m_fDrawOnlyGenB, 
//					   true, false);
	}
	else	// rendering in flattened
	{
		if (params->GetStrValue(renderingParams, "ConfSurf", valid)=="RenderConformalMap1")
		{
			if (m_MapTransform1!=m_pConfM1->GetCurrentTransform())
				m_pConfM1->Transform(m_MapTransform1);
			((SurfaceMidEdgeConf*)GetSurface(0))->DrawConformal(params, renderingParams);
		}
		
		if (params->GetStrValue(renderingParams, "ConfSurf", valid)=="RenderConformalMap2")	
		{
			if (m_MapTransform2!=m_pConfM2->GetCurrentTransform())
				m_pConfM2->Transform(m_MapTransform2);
			((SurfaceMidEdgeConf*)GetSurface(0))->LoadConformalCamera();
			((SurfaceMidEdgeConf*)GetSurface(1))->DrawConformalDoNotLoadCamera(params, renderingParams);		
		}
	}
}

void MapConformal::SaveMap(std::ofstream & textStream)
{
	textStream<<"Format MapConformal\n";
	// save generators (if any)
	textStream<<m_GeneratorCorrespondences.size()<<"\n";
	if (m_GeneratorCorrespondences.size()!=0)
	{
		for (int i=0; i<(int)m_GeneratorCorrespondences.size(); i++)
			textStream<<m_GeneratorCorrespondences[i]<<" ";
		textStream<<"\n";
		textStream<<m_GenSetName.c_str()<<"\n";
	}
	
	// save mobius transformations
	m_MapTransform1.SaveTransformation(textStream);
	m_MapTransform2.SaveTransformation(textStream);
	
	//WriteCache(textStream);
}

void MapConformal::LoadMap(std::ifstream & textStream)
{	
	std::string tempStr;
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="Format");
	textStream>>tempStr;	assert(VKString(tempStr.c_str())=="MapConformal");	
	// save generators (if any)
	int numCorrespondences=0;
	textStream>>numCorrespondences;
//	m_pGeneratorSampleSet1 = new SurfaceSampleSet();
//	m_pGeneratorSampleSet2 = new SurfaceSampleSet();
	
	if (numCorrespondences!=0)
	{
		for (int i=0; i<numCorrespondences; i++)
		{
			int tempVal;
			textStream>>tempVal;
			m_GeneratorCorrespondences.push_back(tempVal);
		}
		m_GenSetName.readToken(textStream);
	}
	
	// save mobius transformations
	m_MapTransform1.LoadTransformation(textStream);
	m_MapTransform2.LoadTransformation(textStream);
	
	//ReadCache(textStream);
}




