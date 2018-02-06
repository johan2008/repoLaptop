#include "PlanarTransform.h"
#include "VKString.h"
#include <fstream>
#include <sstream>
#include "MobiusTransformation.h"
#include "SurfaceMidEdgeConf.h"
#include "PlanarTransformMVC.h"
#include "PlanarTransformLSCM.h"

PlanarTransform::PlanarTransform(R3Mesh * flatMesh,
								 R2Kdtree<FlatSearchNode*> * flatMeshSearch)
{
	m_pFlatMesh = flatMesh;
	m_pFlatMeshSearch = flatMeshSearch;
}

PlanarTransform * PlanarTransform::CreateTransform(std::ifstream & textStream, R3Mesh * flatMesh)
{
	std::streampos pos=textStream.tellg();
	
	std::string tempStr;
	textStream>>tempStr;	VKString transformName(tempStr.c_str());
	textStream.seekg(pos);
	
	PlanarTransform * returnXform=NULL;
	if (transformName=="MobiusTransformationEncoded:")
		returnXform = new MobiusTransformation();
//	else if (transformName=="PlanarTransformLSCM")
//		returnXform = new PlanarTransformLSCM(flatMesh); 
	else if (transformName=="PlanarTransformMVC")
		returnXform = new PlanarTransformMVC();
	assert(returnXform!=NULL);
	returnXform->LoadTransformation(textStream);
	return returnXform;
}
