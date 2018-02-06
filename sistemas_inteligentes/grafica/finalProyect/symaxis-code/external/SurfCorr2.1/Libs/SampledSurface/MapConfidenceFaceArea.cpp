#include "MapConfidenceFaceArea.h"
#include "MapFlattening.h"

MapConfidenceFaceArea::MapConfidenceFaceArea(SurfaceMap * correspondingMap, ValueType valType, 
											 const VKString & sampFrom, const VKString & sampTo)
: SurfaceMapConfidence(correspondingMap, true, false, -1, sampFrom, sampTo, true, false)
{
	m_ValueType = valType;
}

MapConfidenceFaceArea::~MapConfidenceFaceArea()
{
}

VKString MapConfidenceFaceArea::GetMapConfidenceType()
{
	return "MapConfidenceFaceArea";
}

SurfaceMapConfidence * MapConfidenceFaceArea::Clone(SurfaceMap * currMap)
{
	return new MapConfidenceFaceArea(currMap, m_ValueType, m_SampleSet1, m_SampleSet2);
}

double MapConfidenceFaceArea::CalculateConfidenceAtSample(const SurfaceSample & sample)
{
	bool atVertex;
	int nearestVertex = sample.NearestVertex(&atVertex);
	assert(atVertex);
	bool forward = sample.CheckMeshIsSame(m_pSurfaceMap->GetSurface(0)->GetMesh());
	
	double conf = GetPerVertexConfidence(m_ValueType, forward, nearestVertex, m_pSurfaceMap);
	return conf;
}

double MapConfidenceFaceArea::GetPerVertexConfidence(ValueType tp, bool forward, 
													 int vertexID, SurfaceMap * map1)
{
//	if (tp==TRIAREA_VAL_DIFF || tp==TRIAREA_VAL_MIN_OVER_MAX)
		return 1. - GetPerVertexError(tp, forward, vertexID, map1);
	// e ^ -val?
	assert(false);
	return -1;
}

double MapConfidenceFaceArea::GetPerVertexError(ValueType tp, bool forward, 
												int vertexID, SurfaceMap * map1)
{
	R3Mesh * mesh = (forward) ? (map1->GetSurface(0)->GetMesh()) : (map1->GetSurface(1)->GetMesh());
	int numFaces = mesh->VertexValence(mesh->Vertex(vertexID));	
	double cumErr = 0;
	double norm = 0;
//	double cumA1 = 0;
//	double cumA2 = 0;	
	for (int i=0; i<numFaces; i++)
	{
		R3MeshEdge * edge =  mesh->EdgeOnVertex(mesh->Vertex(vertexID), i);
		assert(edge!=NULL);
		R3MeshFace * face = mesh->FaceOnVertex(mesh->Vertex(vertexID), edge);
		if (face!=NULL)
		{
			int triID = mesh->FaceID(face);
			
			double a1, a2;
			double val = 0;
			if (tp==TRIAREA_VAL_DIFF || tp==TRIAREA_VAL_MIN_OVER_MAX)
				val = GetPerTriangleErrorArea(tp, forward, triID, map1, &a1, &a2);
			else if (tp==TRISTRETCH_SV_RATIO || tp==TRISTRETCH_SV_L2STRETCH || tp==TRISTRETCH_SV_AREA)
				val = GetPerTriangleErrorStretch(tp, forward, triID, map1, &a1, &a2);
			else
				assert(false);	// for now
			
			if (val>=0)
			{
//				cumA1 += a1;
//				cumA2 += a2;
				double weight = forward ? a1 : a2;
				cumErr += val * weight;		
				norm += weight;
			}
		}
	}
		
	if (norm>0)
	{
		if (!(cumErr/norm>=0 && cumErr/norm<=1))
		{
			std::cout<<"[ERROR] CumErr="<<cumErr<<" norm="<<norm<<std::endl;
			assert(cumErr/norm>=0 && cumErr/norm<=1);
		}
		
		//std::cout<<"cumErr="<<cumErr<<" norm="<<norm<<" val="<<(cumErr / norm)<<std::endl;
		return cumErr / norm;
	}
	//std::cout<<"norm="<<norm<<std::endl;
//	if (cumA1 > 0 && cumA2 > 0)
//		return 1. - vkMin(cumA1, cumA2) / vkMax(cumA1, cumA2);
	return 1;
}

void MapConfidenceFaceArea::FillSamples(bool forward, int triangleID, SurfaceMap * map1,
										double * A1save, double * A2save,
										SurfaceSample & s11, SurfaceSample & s12, SurfaceSample & s13,
										SurfaceSample & s21, SurfaceSample & s22, SurfaceSample & s23,
										double * Anorm1, double * Anorm2)
{
	R3Mesh * mesh1 = map1->GetSurface(0)->GetMesh();
	R3Mesh * mesh2 = map1->GetSurface(1)->GetMesh();
	*Anorm1 = map1->GetSurface(0)->Area();
	*Anorm2 = map1->GetSurface(1)->Area();
	
	if (map1->GetSurface(1)==Surface2DPlane::m_pInfinitePlanePseudosurface
		&& map1->GetSurfaceMapType()=="MapFlatMidEdge")
	{
		*Anorm2 = ((MapFlattening*)map1)->GetFlatMeshArea();
		//*Anorm2 = 1;
	}
		
	if (forward)
	{
		assert(triangleID < mesh1->NFaces());
		s11 = SurfaceSample(mesh1->VertexID(mesh1->VertexOnFace(mesh1->Face(triangleID), 0)), mesh1);
		s12 = SurfaceSample(mesh1->VertexID(mesh1->VertexOnFace(mesh1->Face(triangleID), 1)), mesh1);
		s13 = SurfaceSample(mesh1->VertexID(mesh1->VertexOnFace(mesh1->Face(triangleID), 2)), mesh1);
		s21 = SurfaceSample(map1->ForwardMap(s11));
		s22 = SurfaceSample(map1->ForwardMap(s12));
		s23 = SurfaceSample(map1->ForwardMap(s13));
		*A1save = mesh1->FaceArea(mesh1->Face(triangleID)) / (*Anorm1);
	}
	else
	{
		assert(triangleID < mesh2->NFaces());		
		s21 = SurfaceSample(mesh2->VertexID(mesh2->VertexOnFace(mesh2->Face(triangleID), 0)), mesh2);
		s22 = SurfaceSample(mesh2->VertexID(mesh2->VertexOnFace(mesh2->Face(triangleID), 1)), mesh2);
		s23 = SurfaceSample(mesh2->VertexID(mesh2->VertexOnFace(mesh2->Face(triangleID), 2)), mesh2);
		s11 = SurfaceSample(map1->InverseMap(s21));
		s12 = SurfaceSample(map1->InverseMap(s22));
		s13 = SurfaceSample(map1->InverseMap(s23));
		*A2save = mesh2->FaceArea(mesh2->Face(triangleID)) / (*Anorm2);		
	}
}

double MapConfidenceFaceArea::GetPerTriangleErrorStretch(ValueType tp, bool forward, 
														 int triangleID, SurfaceMap * map1,
														 double * A1save, double * A2save)
{
	*A1save = map1->GetSurface(0)->Area();
	*A2save = map1->GetSurface(1)->Area();
	double A1 = sqrt(*A1save);
	double A2 = sqrt(*A2save);

	SurfaceSample s11, s12, s13, s21, s22, s23;	
	double Anorm1, Anorm2;
	FillSamples(forward, triangleID, map1, A1save, A2save, s11, s12, s13, s21, s22, s23, &Anorm1, &Anorm2);
	R3Point p11, p12, p13, p21, p22, p23;
	p11 = s11.GetPosition()/A1;		p12 = s12.GetPosition()/A1;		p13 = s13.GetPosition()/A1;
	p21 = s21.GetPosition()/A2;		p22 = s22.GetPosition()/A2;		p23 = s23.GetPosition()/A2;
	if (p11==p12 || p11==p13 || p12==p13 || p21==p22 || p21==p23 || p22==p23)	// area of a triangle is 0
		return 1.;	// max error;
	
	R3Vector u = p12-p11;
	R3Vector x = p22-p21;
	
	R3Vector v = p13-p11;
	R3Vector y = p23-p21;
	
	R3Vector e1 = u;		e1.Normalize();
	R3Vector b1 = x;		b1.Normalize();
	
	R3Vector e2 = v - v.Dot(e1) * e1;		e2.Normalize();
	R3Vector b2 = y - y.Dot(b1) * b1;		b2.Normalize();
		
	LinAlgMatrixReal B(2, 2), E(2, 2), Einv(2,2);
	B(0,0) = x.Dot(b1);
	B(0,1) = y.Dot(b1);
	B(1,0) = x.Dot(b2);
	B(1,1) = y.Dot(b2);
	
	E(0,0) = u.Dot(e1);
	E(0,1) = v.Dot(e1);
	E(1,0) = u.Dot(e2);
	E(1,1) = v.Dot(e2);
	
//	E.PrintMatlab();	std::cout<<std::endl;
//	B.PrintMatlab();	std::cout<<std::endl;
	//E.Inverse();
	double determ = E(0,0) * E(1,1) - E(0,1)*E(1,0);
	if (determ==0)
		return 1.;
	Einv(0,0) = E(1,1)/determ;
	Einv(0,1) = -E(0,1)/determ;
	Einv(1,0) = -E(1,0)/determ;
	Einv(1,1) = E(0,0)/determ;
	LinAlgMatrixReal abc(2,2);
	B.Multiply(Einv, abc);
	// do svd on matrix abc
	LinAlgMatrixReal U(2,2);
	LinAlgVectorReal S(2);
	LinAlgMatrixReal Vt(2,2);

	abc.SVDDecomposeDestructive(U, S, Vt);
	if (S(0) < S(1))
	{
		std::cout<<"[ERROR] Singular values not sorted: "<<S(0)<<" < "<<S(1)<<std::endl;
		assert(S(0) >= S(1));
	}
	double sigma2 = pow(1., 2.);
	double distortion = 0;
	if (tp==TRISTRETCH_SV_RATIO)	// only measures 'conformal error'
		distortion = ((S(0) / S(1)) - 1) / (1 + S(0) / S(1));
	else if (tp==TRISTRETCH_SV_L2STRETCH)
		distortion = S(0)*S(0) + S(1)*S(1) + 1. / (S(0)*S(0)) + 1. / (S(1)*S(1));
	else if (tp==TRISTRETCH_SV_AREA)
		distortion = (S(0)*S(1) + 1/(S(0) * S(1))) / 2.;
	else
		assert(false);
	return 1. - 1. / pow(distortion, 1.);
}

double MapConfidenceFaceArea::GetPerTriangleErrorArea(ValueType tp, bool forward, 
													  int triangleID, SurfaceMap * map1,
													  double * A1save, double * A2save)
{
	R3Mesh * mesh1 = map1->GetSurface(0)->GetMesh();
	R3Mesh * mesh2 = map1->GetSurface(1)->GetMesh();
	double A1=-1;
	double A2=-1;
	double Anorm1, Anorm2;
	SurfaceSample s11, s12, s13, s21, s22, s23;
		
	FillSamples(forward, triangleID, map1, &A1, &A2, s11, s12, s13, s21, s22, s23, &Anorm1, &Anorm2);
	
	if ((A2==-1 && !s21.Invalid() && !s22.Invalid() && !s23.Invalid())
		||(A1==-1 && !s11.Invalid() && !s12.Invalid() && !s13.Invalid()))
	{
		if (A2==-1)
		{
			if (s21.GetPosition()==s22.GetPosition() || s22.GetPosition()==s23.GetPosition()
				|| s21.GetPosition()==s23.GetPosition())
			{
				A2 = 0;
			}
			else
			{				
				R3TriangleVertex v1(s21.GetPosition());
				R3TriangleVertex v2(s22.GetPosition());
				R3TriangleVertex v3(s23.GetPosition());
				A2 = R3Triangle(&v1, &v2, &v3).Area() / Anorm2;
			}
		}
		else if (A1==-1)
		{
			if (s11.GetPosition()==s12.GetPosition() || s12.GetPosition()==s13.GetPosition()
				|| s11.GetPosition()==s13.GetPosition())
				A1 = 0;
			else
			{				
				R3TriangleVertex v1(s11.GetPosition());
				R3TriangleVertex v2(s12.GetPosition());
				R3TriangleVertex v3(s13.GetPosition());
				A1 = R3Triangle(&v1, &v2, &v3).Area() / Anorm1;			
			}
		}
		if (A1save!=NULL)
			*A1save = A1;
		if (A2save!=NULL)
			*A2save = A2;
		double dist=0;

		if (vkMax(A1,A2)==0)
			return 1;
		if (tp==TRIAREA_VAL_MIN_OVER_MAX)
		{
			if (A1==A2)
				dist = 0;
			else
				dist = 1 - vkMin(A1, A2) / vkMax(A1, A2);
		}
		else if (tp==TRIAREA_VAL_DIFF)
		{
			assert(false);
			dist = A1-A2;
		}
		else
			assert(false);
		
		assert(dist>=0 && dist<=1);
		return dist * dist;
	}
	
//	std::cout<<"Triangle = "<<triangleID<<" invalid"<<std::endl;
	return 1;
}
