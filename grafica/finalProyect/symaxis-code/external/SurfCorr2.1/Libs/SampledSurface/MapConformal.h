#ifndef __MAP_CONFORMAL_H
#define __MAP_CONFORMAL_H

#include <fstream>
#include <sstream>

#include "gaps.h"

#include "VKString.h"

#include "SampledSurface.h"
#include "SurfaceMap.h"
#include "SurfaceMidEdgeConf.h"

using namespace std;

class MapCoarse;
class MapMultiConformal;

class MapConformal : public SurfaceMap
{
	public:
		friend class MapMultiConformal;
	
		MapConformal(SurfaceMidEdgeConf * M1, SurfaceMidEdgeConf * M2, 
					 const std::vector<int> & correspondences,
					 const VKString & genSetName);
		
		MapConformal(SurfaceMidEdgeConf * M1, SurfaceMidEdgeConf * M2, 
					 const MobiusTransformation & m1, const MobiusTransformation & m2);
	
		void UpdateGeneratingCorrsChanged();

		virtual SurfaceSample ForwardMap(const SurfaceSample & s);
		virtual SurfaceSample InverseMap(const SurfaceSample & s);

		virtual VKString GetSurfaceMapType();
	
		virtual void SetDrawOnlyGenerators(double r, double g, double b);
		virtual void DrawGenerators(AnalysisWindow * window,
									ParamParser * params, const VKString & renderingParams, 
									const VKString & surfaceName, double linewidth, 
									double r, double g, double b, 
									bool lines = true, bool pnts = true);
		virtual void Draw(AnalysisWindow * window, ParamParser * params,
						  const VKString & renderingParams="RendererDefault", 
						  const VKString & surfaceName="none");
	
		virtual void SaveMap(std::ofstream & textStream);
		virtual void LoadMap(std::ifstream & textStream);
	
		void GetGeneratingCorrespondences(std::vector<int> ** correspondences,
										  SurfaceSampleSet ** sampleSet1,
										  SurfaceSampleSet ** sampleSet2);
		const MobiusTransformation & GetMobiusForSurface(int surfaceID) const; 
	
	protected:
		static SurfaceSample MapSample(const SurfaceSample & s, 
									   const MobiusTransformation & t1, 
									   const MobiusTransformation & t2,
									   SurfaceMidEdgeConf * surf1, 
									   SurfaceMidEdgeConf * surf2);
		
		SurfaceMidEdgeConf * m_pConfM1;
		SurfaceMidEdgeConf * m_pConfM2;
		MobiusTransformation m_MapTransform1;
		MobiusTransformation m_MapTransform2;

		VKString m_GenSetName;
		std::vector<int> m_GeneratorCorrespondences;
	
		double m_fDrawOnlyGenR, m_fDrawOnlyGenG, m_fDrawOnlyGenB; 
};

#endif


