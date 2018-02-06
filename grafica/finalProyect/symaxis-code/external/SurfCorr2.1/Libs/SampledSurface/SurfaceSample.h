
#include "gaps.h"
#include "VKString.h"

#ifndef __SURFACE_SAMPLE_H
#define __SURFACE_SAMPLE_H

class SurfaceSample
{
	public:
		friend class SurfaceSampleSet;
	
		SurfaceSample();
		SurfaceSample(const SurfaceSample & other);
		SurfaceSample(R3MeshIntersection & intersection, R3Point query, R3Mesh * mesh,
					  const VKString & type="none");
		SurfaceSample(int triangleID, double b1, double b2, double b3, R3Mesh * mesh,  const VKString & type="none");
		SurfaceSample(int triangleID, R3Point posOnTriangle, R3Mesh * mesh,  const VKString & type="none");
		SurfaceSample(int vertexID, R3Mesh * mesh,  const VKString & type="none");

		virtual ~SurfaceSample(){}

		bool operator == (const SurfaceSample & other) const;
		void operator = (const SurfaceSample & other);

		bool Invalid() const;
		int TriID() const;
		double B(int coordID) const;
		int VertID(int coordID) const;
		int NearestVertex(bool * atVertex=NULL) const;
		double Interpolate(double b1, double b2, double b3) const;
	
		bool CheckMeshIsSame(const R3Mesh * mesh) const;

		R3Point GetPosition() const;
		R3Vector Normal() const;

		VKString GetSampleType() const;
		void SetSampleType(const VKString & type);
 
	protected:
		bool InitializedSpecialMesh(R3Mesh * mesh, double b1, double b2, double b3);
	
		int m_iTriangleID;
		double m_fBarycentricCoords[3]; 
		R3Mesh * m_pMesh;
		VKString m_sType;
};

#endif
