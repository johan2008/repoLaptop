#include <assert.h>
#include <iostream>

#include "MapCoarse.h"
#include "SampledSurface.h"
#include "DistanceGeodesic.h"
#include "MapGeoFeature.h"
#include "MeshProcessor.h"
#include "gaps.h"
#include "TimeProfiler.h"


void PrintUsage(const char * progname);
void MapFromTOSCAVert(const VKString & inputDirVert1, const VKString & vertFile1,
					  const VKString & inputDirMesh1, const VKString & meshName1,
					  const VKString & inputDirVert2, const VKString & vertFile2,
					  const VKString & inputDirMesh2, const VKString & meshName2,					  
					  const VKString & outputDir);
void FixNoNumFacesOff(const VKString & inputDirOff, const VKString & inputMeshName, 
					  const VKString & outputDir);
void PairwiseGeodesicDistances(const VKString & inputMeshName, 
							   const VKString & outputDistances);
void ConvertCorrsBronToFunk(const VKString & simplifiedMeshFrom,
							const VKString & simplifiedMeshTo,
							const VKString & trianglesFrom,
							const VKString & baryCoordsFrom,
							const VKString & trianglesTo,
							const VKString & baryCoordsTo,
							const VKString & originalMeshFrom,
							const VKString & originalMeshTo,
							const VKString & outputCorFile);

void ConvertMeshCorrs(const VKString & modifiedMeshName1, const VKString & modifiedMeshName2,
					  const VKString & origMeshName1, const VKString & origMeshName2,
					  const VKString & simpleCorrFile, const VKString & outputOrigMeshname);

void ConvertDefDrivCorsZhang08(const VKString & smfFileFrom, const VKString & smfFileTo,
							   const VKString & corrZhangFile, const VKString & origFromFile,
							   const VKString & origToFile, const VKString & funkCoarseCorFile);

void Coarse2FineByProjection(const VKString & meshCoarse1, const VKString & meshCoarse2,
							 const VKString & fullCoarseMap,
							 const VKString & meshFine1, const VKString & meshFine2, 
							 const VKString & fullFineMap);

int main (int argc, const char ** argv)
{
	if (argc<1)
		PrintUsage(argv[0]);
	else if (strcmp(argv[1], "CheckForConformalProblems")==0)
	{
		R3Mesh mesh;
		mesh.ReadFile(argv[2]);
			
		int numConn = MeshProcessor::NumConnectedComponents(&mesh);
		int numIsoVert = MeshProcessor::NumIsolatedVertices(&mesh);
		
	//	std::cout<<"Examining: "<<argv[2]<<std::endl;
		if (numConn>1)
			std::cout<<"\t[ERROR] Connected Components = "<<numConn<<std::endl;
		if (numIsoVert>0)
			std::cout<<"\t[ERROR] Isolated Vertices = "<<numIsoVert<<std::endl;		
	}
	else if (strcmp(argv[1], "ConvertDefDrivCorsZhang08")==0)
	{
		ConvertDefDrivCorsZhang08(VKString(argv[2]), VKString(argv[3]), 
								  VKString(argv[4]), VKString(argv[5]), 
								  VKString(argv[6]), VKString(argv[7]));
	}
	else if (strcmp(argv[1], "ConvertMeshCorrs")==0)
	{
		ConvertMeshCorrs(VKString(argv[2]), VKString(argv[3]), 
						 VKString(argv[4]), VKString(argv[5]), 
						 VKString(argv[6]), VKString(argv[7]));		
	}
	else if (strcmp(argv[1], "Coarse2FineByProjection")==0)
	{
		Coarse2FineByProjection(VKString(argv[2]), VKString(argv[3]), 
								VKString(argv[4]), VKString(argv[5]), 
								VKString(argv[6]), VKString(argv[7]));		
	}
	else if (strcmp(argv[1], "ConvertCorrsBronToFunk")==0)
	{
		ConvertCorrsBronToFunk(VKString(argv[2]), VKString(argv[3]), VKString(argv[4]), 
							   VKString(argv[5]), VKString(argv[6]), VKString(argv[7]),
							   VKString(argv[8]), VKString(argv[9]), VKString(argv[10]));
	}
	else if (strcmp(argv[1], "PairwiseGeodesicDistances")==0)
	{
		PairwiseGeodesicDistances(VKString(argv[2]), VKString(argv[3]));
	}
	else if (strcmp(argv[1], "SymCoarseToFine")==0)
	{
		assert(false);
//		assert(argc==5 || argc==6);
//		R3Mesh symMesh;
//		symMesh.ReadFile(argv[2]);
//		SampledSurface surf(&symMesh);
//		DistanceGeodesic geoDist(&surf);		
//		std::cout<<"Creating Distance Metric"<<std::endl;		
//		geoDist.PrecomputeDistances();
//		surf.AddDistanceMetric("default", &geoDist);
//		
//		MapCoarse mapCoarse(&surf, &surf);
//
//		if (argc==5)
//			((SurfaceMap&)mapCoarse).LoadMap(argv[3]);
//		else
//			mapCoarse.LoadCorrFeaturePntSymmetryMap(argv[3], argv[4]);
//
//		int offset=0;
//		if (argc==6)
//			offset++;
//		
//		MapGeoFeature * mapFine = mapCoarse.GetInterpolatedGeoFeature(1, "default", "default");
//		std::cout<<"Finding Dense Correspondences"<<std::endl;
//		
//		MapCoarse * resultingFineMap = new MapCoarse(new SurfaceSampleSet(mapFine->GetSurface(0)->GetMesh()),
//													 new SurfaceSampleSet(mapFine->GetSurface(1)->GetMesh()),
//													 mapFine, MapCoarse::F2C_FORWARD_ONLY);
//		std::cout<<"Writing Dense Correspondences: "<<argv[4+offset]<<std::endl;		
//		((SurfaceMap*)resultingFineMap)->SaveMap(argv[4+offset]);
	}
	else if (strcmp(argv[1], "MeshTypeCast")==0)
	{
		assert(argc==7);
		VKString dirIn(argv[2]);
		VKString fileName(argv[3]);
		VKString extIn(argv[4]);
		VKString extOut(argv[5]);
		VKString dirOut(argv[6]);

		R3Mesh mesh;
		if (extIn!="TOSCA")
		{
			mesh.ReadFile((dirIn+fileName+extIn).c_str());
			mesh.WriteFile((dirOut+fileName+extOut).c_str());
		}
		else
		{
			std::vector<double> vertexCoords;
			std::vector<int> faceIDs;
			double x, y, z;
			ifstream verticesStream((dirIn+fileName+".vert").c_str());
			assert(verticesStream.is_open());
			while(!verticesStream.eof())
			{
				verticesStream>>x>>y>>z;
				vertexCoords.push_back(x);
				vertexCoords.push_back(y);
				vertexCoords.push_back(z);				
				mesh.CreateVertex(R3Point(x, y, z));
			}
			
			ifstream trianglesStream((dirIn+fileName+".tri").c_str());
			assert(trianglesStream.is_open());
			int v1, v2, v3;
			while(!trianglesStream.eof())
			{
				trianglesStream>>v1>>v2>>v3;
				faceIDs.push_back(v1-1);
				faceIDs.push_back(v2-1);
				faceIDs.push_back(v3-1);				
				mesh.CreateFace(mesh.Vertex(v1-1), mesh.Vertex(v2-1), mesh.Vertex(v3-1));
			}
			
			if (extOut==".off")
			{
				std::ofstream textStream((dirOut+fileName+extOut).c_str());
				textStream<<"OFF\n";
				textStream<<vertexCoords.size()/3<<" "<<faceIDs.size()/3<<" 0\n";
				for (int i=0; i<(int)vertexCoords.size()/3; i++)
					textStream<<vertexCoords[3*i+0]<<" "<<vertexCoords[3*i+1]<<" "<<vertexCoords[3*i+2]<<"\n";
				for (int i=0; i<(int)faceIDs.size()/3; i++)
					textStream<<"3 "<<faceIDs[3*i+0]<<" "<<faceIDs[3*i+1]<<" "<<faceIDs[3*i+2]<<"\n";
				textStream.close();
			}
			else 
				mesh.WriteFile((dirOut+fileName+extOut).c_str());
		}
		
	}
	else if (strcmp(argv[1], "SmoothMesh")==0)
	{
		assert(argc==6);
		VKString inputDir(argv[2]);
		VKString meshName(argv[3]);
		int numIters = VKString(argv[4]).toInt();
		VKString outputDir(argv[5]);		
		MeshProcessor proc(new R3Mesh());
		proc.m_pMesh->ReadFile((inputDir+meshName).c_str());
		proc.Smooth(numIters);
		proc.m_pMesh->WriteFile((outputDir+meshName).c_str());
	}
	else if (strcmp(argv[1], "FixNoNumFacesOff")==0)
	{
		assert(argc==5);
		VKString inputDirOff(argv[2]);
		VKString inputMeshName(argv[3]);		
		VKString outputDir(argv[4]);	
		FixNoNumFacesOff(inputDirOff, inputMeshName, outputDir);
	}
	else if (strcmp(argv[1], "MapFromTOSCAVert")==0)
	{
		assert(argc==11);
		VKString inputDirVert1(argv[2]);
		VKString vertFile1(argv[3]);
		VKString inputDirMesh1(argv[4]);		
		VKString meshName1(argv[5]);
		VKString inputDirVert2(argv[6]);
		VKString vertFile2(argv[7]);
		VKString inputDirMesh2(argv[8]);		
		VKString meshName2(argv[9]);
		VKString outDir(argv[10]);	
		MapFromTOSCAVert(inputDirVert1, vertFile1, inputDirMesh1, meshName1, 
						 inputDirVert2, vertFile2, inputDirMesh2, meshName2, 
						 outDir);
	}	
	else if (strcmp(argv[1], "DelaunayTriangEdgeFlip")==0)
	{
		assert(argc==4);
		R3Mesh inMesh;
		inMesh.ReadFile(argv[2]);
		MeshProcessor processor(&inMesh);
		processor.DelaunayTriangByEdgeFlip();
		processor.m_pMesh->WriteFile(argv[3]);
	}	
	else if (strcmp(argv[1], "SingleConnectedComponent")==0)
	{
		assert(false);
		assert(argc==5);
		VKString inputDir(argv[2]);
		VKString meshName(argv[3]);
		VKString outputDir(argv[4]);
		MeshProcessor proc(new R3Mesh());
		proc.m_pMesh->ReadFile((inputDir+meshName).c_str());
		proc.SingleConnectedComponent();
		proc.m_pMesh->WriteFile((outputDir+meshName).c_str());
	}
	else if (strcmp(argv[1], "ConsistentFaceOrientation")==0)
	{
		assert(false);
		assert(argc==5);		
		VKString inputDir(argv[2]);
		VKString meshName(argv[3]);
		VKString outputDir(argv[4]);
		MeshProcessor proc(NULL);
		proc.ConsistentFaceOrientation((inputDir + meshName).c_str());
		proc.m_pMesh->WriteFile((outputDir + meshName).c_str());
	}
	else
	{
		PrintUsage(argv[0]);
	}
}  

void PrintUsage(const char * progname)
{
	std::cout<<"\nUsages (1):\n\t"<<progname<<" SymCoarseToFine meshNameIn [coarseMapNameIn OR coarseMapFeatPntsToVrtx coarseMapFeatPntToPnt] fineMapNameOut\n"<<std::endl;	
	std::cout<<"       (2):\n\t"<<progname<<" MeshTypeCast dirIn fileName extIn extOut dirOut\n"<<std::endl;		
	std::cout<<"       (3):\n\t"<<progname<<" SmoothMesh inputDir meshName smoothIter outputDir\n"<<std::endl;		
	std::cout<<"       (4):\n\t"<<progname<<" MapFromTOSCAVert inputDirVert1 vertFile1 inputDirMesh1 meshName1 ";
	std::cout<<"inputDirVert2 vertFile2 inputDirMesh2 meshName2 outputDir\n"<<std::endl;
	std::cout<<"       (5):\n\t"<<progname<<" DelaunayTriangEdgeFlip meshIn meshOut\n"<<std::endl;		
	std::cout<<"       (6):\n\t"<<progname<<" Coarse2FineByProjection meshCoarse1 meshCoarse2 fullCoarseMap meshFine1 meshFine2 fullFineMap\n"<<std::endl;			
	
//	std::cout<<"       (5):\n\t"<<progname<<" SingleConnectedComponent inputDir meshName outputDir\n"<<std::endl;			
//	std::cout<<"       (6):\n\t"<<progname<<" ConsistentFaceOrientation inputDir meshName outputDir\n"<<std::endl;				
}

struct TriVertex
{
	TriVertex(int id, double x, double y, double z)
	:id(id), pos(x, y, z)
	{	}
	
	int id;
	R3Point pos;
};

R3Point TriVertexPositionGlobalFn(TriVertex * v, void *)
{
	return v->pos;
}

void MapFromTOSCAVert(const VKString & inputDirVert1, const VKString & vertFile1,
					  const VKString & inputDirMesh1, const VKString & meshName1,
					  const VKString & inputDirVert2, const VKString & vertFile2,
					  const VKString & inputDirMesh2, const VKString & meshName2,
					  const VKString & outputDir)
{
	// load vertices from into kd tree
	RNArray<TriVertex*> vertexCoords1;
	double x, y, z;
	ifstream verticesStream1((inputDirVert1+vertFile1+".vert").c_str());
	assert(verticesStream1.is_open());
	while(!verticesStream1.eof())
	{
		verticesStream1>>x>>y>>z;
		vertexCoords1.InsertTail(new TriVertex(vertexCoords1.NEntries(), x, y, z));
	}
	R3Kdtree<TriVertex*> kdTree1(vertexCoords1, &TriVertexPositionGlobalFn, NULL);
	
	// load vertices 2
	RNArray<TriVertex*> vertexCoords2;
	ifstream verticesStream2((inputDirVert2+vertFile2+".vert").c_str());
	assert(verticesStream2.is_open());
	while(!verticesStream2.eof())
	{
		verticesStream2>>x>>y>>z;
		vertexCoords2.InsertTail(new TriVertex(vertexCoords2.NEntries(), x, y, z));
	}
	
	R3Mesh mesh1;
	mesh1.ReadFile((inputDirMesh1+meshName1).c_str());

	R3Mesh mesh2;
	mesh2.ReadFile((inputDirMesh2+meshName2).c_str());

	//SampledSurface surf1(&mesh1);
	//SampledSurface surf2(&mesh2);	
	
	R3MeshSearchTree searchTreeOnM2(&mesh2);
	
	std::cout<<"Mapping files: "<<mesh1.NVertices()<<std::endl;
	for (int i=0; i<mesh1.NVertices(); i++)
	{
		int vertexCloudID = kdTree1.FindClosest(mesh1.VertexPosition(mesh1.Vertex(i)))->id;
		R3MeshIntersection closest;
		searchTreeOnM2.FindClosest(vertexCoords2[vertexCloudID]->pos, closest);		
		
		assert(closest.type!=R3_MESH_NULL_TYPE);
		
		R3Point bary = mesh2.FaceBarycentric(closest.face, closest.point);
	
		int otherVertexID=-1;
		if (bary.X() > bary.Y() && bary.X() > bary.Z())
			otherVertexID = mesh2.VertexID(mesh2.VertexOnFace(closest.face, 0));
		else if (bary.Y() > bary.X() && bary.Y() > bary.Z())
			otherVertexID = mesh2.VertexID(mesh2.VertexOnFace(closest.face, 1));
		else if (bary.Z() > bary.X() && bary.Z() > bary.Y())
			otherVertexID = mesh2.VertexID(mesh2.VertexOnFace(closest.face, 2));
		else
			assert(false);
		
		if (i!=otherVertexID)
			std::cout<<"Map["<<i<<"] -> "<<otherVertexID<<std::endl;
	}
}

void FixNoNumFacesOff(const VKString & inputDirOff, const VKString & inputMeshName, 
					  const VKString & outputDir)
{
	std::ofstream outputStream((outputDir+inputMeshName).c_str());
	std::ifstream textStream((inputDirOff+inputMeshName).c_str());
	assert(textStream.is_open());
	std::string tempStr;
	textStream>>tempStr;
	assert(strcmp(tempStr.c_str(), "OFF")==0);
	outputStream<<"OFF\n";
	int vertices, faces, edges;
	textStream>>vertices>>faces>>edges;
	outputStream<<vertices<<" "<<faces<<" "<<edges<<"\n";
	double x, y, z;
	for (int i=0; i<vertices; i++)
	{
		textStream>>x>>y>>z;
		outputStream<<x<<" "<<y<<" "<<z<<"\n";
	}
	
	int v1, v2, v3;
	for (int i=0; i<faces; i++)
	{
		textStream>>v1>>v2>>v3;
		outputStream<<3<<" "<<(v1-1)<<" "<<(v2-1)<<" "<<(v3-1)<<"\n";
	}
}

void ConvertMeshCorrs(const VKString & modifiedMeshName1, const VKString & modifiedMeshName2,
					  const VKString & origMeshName1, const VKString & origMeshName2,
					  const VKString & modifiedCorrFile, const VKString & outputOrigCorrFile)
{
	R3Mesh modifiedMesh1;
	R3Mesh modifiedMesh2;	
	R3Mesh origMesh1;
	R3Mesh origMesh2;	
	
	modifiedMesh1.ReadFile(modifiedMeshName1.c_str());
	modifiedMesh2.ReadFile(modifiedMeshName2.c_str());
	origMesh1.ReadFile(origMeshName1.c_str());
	origMesh2.ReadFile(origMeshName2.c_str());
	
	std::ofstream textStreamOut(outputOrigCorrFile.c_str());	
	assert(textStreamOut.is_open());
	std::ifstream textStreamIn(modifiedCorrFile.c_str());
	assert(textStreamIn.is_open());
	while(!textStreamIn.eof())
	{
		int v1, v2;
		textStreamIn>>v1>>v2;
		if (textStreamIn.eof())
			break;
		R3Point p1 = modifiedMesh1.VertexPosition(modifiedMesh1.Vertex(v1));
		R3Point p2 = modifiedMesh2.VertexPosition(modifiedMesh2.Vertex(v2));
		int v1new=-1, v2new=-1;
		double distv1=-1, distv2=-1;
		for (int i=0; i<origMesh1.NVertices(); i++)
		{
			R3Point p1new = origMesh1.VertexPosition(origMesh1.Vertex(i));
			double dist = (p1new-p1).Length();
			if (v1new==-1 || distv1 > dist)
			{
				distv1 = dist;
				v1new = i;
			}
		}
		for (int i=0; i<origMesh2.NVertices(); i++)
		{
			R3Point p2new = origMesh2.VertexPosition(origMesh2.Vertex(i));
			double dist = (p2new-p2).Length();
			if (v2new==-1 || distv2 > dist)
			{
				distv2 = dist;
				v2new = i;
			}
		}	
		textStreamOut<<v1new<<"\n"<<v2new<<"\n";
		std::cout<<"corrs: "<<v1<<"->"<<v2<<". Changed to: "<<v1new<<"->"<<v2new<<std::endl;		
	}
	textStreamOut.close();
	
}

void PairwiseGeodesicDistances(const VKString & inputMeshName, 
							   const VKString & outputDistances)
{
	std::cout<<"Precomputing pairwise distances"<<std::endl;
	std::cout<<"\tIN="<<inputMeshName.c_str()<<std::endl;	
	std::cout<<"\tOUT="<<outputDistances.c_str()<<std::endl;		
	R3Mesh mesh;
	mesh.ReadFile(inputMeshName.c_str());
	
	double area = 0;
	for (int i=0; i<mesh.NFaces(); i++)
		area += mesh.FaceArea(mesh.Face(i));
	double sqrtArea = sqrt(area);
	
	DistanceGeodesic::GeoVertexData * tempData = DistanceGeodesic::InitializeFunkPrecomputation(&mesh, true);	
	FILE * binFile = fopen(outputDistances.c_str(), "wb");
	assert(binFile!=NULL);
	double * perVertexDists = new double[mesh.NVertices()];
	std::cout<<"Iterating over vertices: "<<mesh.NVertices()<<std::endl;
	TimeProfiler profiler;
	for (int i=0; i<mesh.NVertices(); i++)
	{
		if (i%1000==0)
			profiler.WriteProgress("GeoDistnaces", i, mesh.NVertices());
		DistanceGeodesic::PrecomputeFillFunkhouser(tempData, &mesh, i, perVertexDists);
		for (int i=0; i<mesh.NVertices(); i++)
			perVertexDists[i] /= sqrtArea;
		fwrite (perVertexDists, sizeof(double), mesh.NVertices(), binFile );		
	}
	
	std::cout<<std::endl;
	
	fclose(binFile);
}

void ConvertCorrsBronToFunk(const VKString & simplifiedMeshFrom,
							const VKString & simplifiedMeshTo,
							const VKString & trianglesFrom,
							const VKString & baryCoordsFrom,
							const VKString & trianglesTo,
							const VKString & baryCoordsTo,
							const VKString & originalMeshFrom,
							const VKString & originalMeshTo,
							const VKString & outputCorFile)
{
	// load meshes
	R3Mesh simpleMesh1;
	R3Mesh simpleMesh2;
	R3Mesh origMesh1;
	R3Mesh origMesh2;
	simpleMesh1.ReadFile(simplifiedMeshFrom.c_str());
	simpleMesh2.ReadFile(simplifiedMeshTo.c_str());	
	origMesh1.ReadFile(originalMeshFrom.c_str());
	origMesh2.ReadFile(originalMeshTo.c_str());
	
	// find 3d position of each sample on simplified meshes
	std::vector<R3Point> coarseFrom;
	std::vector<R3Point> coarseTo;	

	std::ifstream textStreamTri1(trianglesFrom.c_str());
	if(!textStreamTri1.is_open())
		std::cout<<"[ERROR] Could not open "<<trianglesFrom.c_str()<<std::endl;
	std::ifstream textStreamBary1(baryCoordsFrom.c_str());
	if(!textStreamBary1.is_open())
		std::cout<<"[ERROR] Could not open "<<baryCoordsFrom.c_str()<<std::endl;
	std::ifstream textStreamTri2(trianglesTo.c_str());
	if(!textStreamTri2.is_open())
		std::cout<<"[ERROR] Could not open "<<trianglesTo.c_str()<<std::endl;
	std::ifstream textStreamBary2(baryCoordsTo.c_str());
	if(!textStreamBary2.is_open())
		std::cout<<"[ERROR] Could not open "<<baryCoordsTo.c_str()<<std::endl;
	assert(textStreamTri1.is_open() && textStreamBary1.is_open()
		   && textStreamTri2.is_open() && textStreamBary2.is_open());
	
	while (!textStreamTri1.eof())
	{
		assert(!textStreamTri1.eof());
		assert(!textStreamBary1.eof());
		assert(!textStreamTri2.eof());
		assert(!textStreamBary2.eof());		
		int tri1ID;
		double b1[3];
		textStreamTri1>>tri1ID;
		textStreamBary1>>b1[0]>>b1[1];
		b1[2] = 1. - b1[0] - b1[1];

		int tri2ID;
		double b2[3];
		textStreamTri2>>tri2ID;
		textStreamBary2>>b2[0]>>b2[1];
		b2[2] = 1. - b2[0] - b2[1];
		
		coarseFrom.push_back(simpleMesh1.FacePoint(simpleMesh1.Face(tri1ID-1), b1));
		coarseTo.push_back(simpleMesh2.FacePoint(simpleMesh2.Face(tri2ID-1), b2));
	}
	
	std::vector<int> vFrom;
	std::vector<int> vTo;	
	// for each pair of 3d positions (in correspondence)
	for (int i=0; i<(int)coarseFrom.size(); i++)
	{
		R3Point p = coarseFrom[i];
		double minDist = FLT_MAX;
		int closestVertex=-1;
		// find nearest vertex (tom's cor file)		
		for (int j=0; j<origMesh1.NVertices(); j++)
		{
			double dist = (origMesh1.VertexPosition(origMesh1.Vertex(j))-p).Length();
			if (minDist > dist)
			{
				minDist = dist;
				closestVertex=j;
			}
		}
		
		assert(closestVertex>=0);
		vFrom.push_back(closestVertex);
	}
	
	// for each pair of 3d positions (in correspondence)
	for (int i=0; i<(int)coarseTo.size(); i++)
	{
		R3Point p = coarseTo[i];
		double minDist = FLT_MAX;
		int closestVertex=-1;
		// find nearest vertex (tom's cor file)		
		for (int j=0; j<origMesh2.NVertices(); j++)
		{
			double dist = (origMesh2.VertexPosition(origMesh2.Vertex(j))-p).Length();
			if (minDist > dist)
			{
				minDist = dist;
				closestVertex=j;
			}
		}
		
		assert(closestVertex>=0);
		vTo.push_back(closestVertex);
	}
	
	// write cor file
	std::ofstream textStream(outputCorFile.c_str());
	assert(textStream.is_open());
	assert(vFrom.size()==vTo.size());
	for (int i=0; i<(int)vFrom.size(); i++)
		textStream<<vFrom[i]<<"\n"<<vTo[i]<<"\n";
	textStream.close();
}

void ConvertDefDrivCorsZhang08(const VKString & smfFileFrom, const VKString & smfFileTo,
							   const VKString & corrZhangFile, const VKString & origFromFile,
							   const VKString & origToFile, const VKString & funkCoarseCorFile)
{
	double x, y, z;
	std::ifstream inputSmf1(smfFileFrom.c_str());
	assert(inputSmf1.is_open());
	std::ifstream inputSmf2(smfFileTo.c_str());
	assert(inputSmf2.is_open());
	std::ifstream inputZhangCorr(corrZhangFile.c_str());
	assert(inputZhangCorr.is_open());
	std::ifstream inputOrigOff1(origFromFile.c_str());
	assert(inputOrigOff1.is_open());
	std::ifstream inputOrigOff2(origToFile.c_str());	
	assert(inputOrigOff2.is_open());
	std::ofstream outputCoarseCorrs(funkCoarseCorFile.c_str());	
	assert(outputCoarseCorrs.is_open());
	
	std::vector<double*> coordsSmf1;
	while (!inputSmf1.eof())
	{
		char type;
		inputSmf1>>type;
		if (type=='v')
		{
			double * vals = new double[3];
			coordsSmf1.push_back(vals);
			inputSmf1>>vals[0]>>vals[1]>>vals[2];
		}
	}
	
	std::vector<double*> coordsSmf2;
	while (!inputSmf2.eof())
	{
		char type;
		inputSmf2>>type;
		if (type=='v')
		{
			double * vals = new double[3];
			coordsSmf2.push_back(vals);
			inputSmf2>>vals[0]>>vals[1]>>vals[2];
		}
	}	
	
	std::string tempStr;
	inputOrigOff1>>tempStr;			assert(strcmp(tempStr.c_str(), "OFF")==0);
	int nV1, nF1, nE1;
	inputOrigOff1>>nV1>>nF1>>nE1;
	RNArray<TriVertex*> vertexCoords1;
	for (int i=0; i<nV1; i++)
	{
		inputOrigOff1>>x>>y>>z;
		vertexCoords1.InsertTail(new TriVertex(vertexCoords1.NEntries(), x, y, z));
	}
	R3Kdtree<TriVertex*> kdTree1(vertexCoords1, &TriVertexPositionGlobalFn, NULL);
	
	inputOrigOff2>>tempStr;			assert(strcmp(tempStr.c_str(), "OFF")==0);
	int nV2, nF2, nE2;
	inputOrigOff2>>nV2>>nF2>>nE2;
	RNArray<TriVertex*> vertexCoords2;
	for (int i=0; i<nV2; i++)
	{
		inputOrigOff2>>x>>y>>z;
		vertexCoords2.InsertTail(new TriVertex(vertexCoords2.NEntries(), x, y, z));
	}
	R3Kdtree<TriVertex*> kdTree2(vertexCoords2, &TriVertexPositionGlobalFn, NULL);
	
	std::map<int, int> corrs;
	VKString lineInCorrFile;
	while (!inputZhangCorr.eof())
	{
		lineInCorrFile.readLine(inputZhangCorr);
		if (lineInCorrFile=="############### Deformation call ###############")
		{
			corrs.clear();	// only take the last one
			std::cout<<"RESTARTING!"<<std::endl;
		}
		else if (lineInCorrFile.startsWith("Feature pair"))
		{
			lineInCorrFile.replace("Feature pair ", "");
			VKStringList tokens = lineInCorrFile.split(" ");
			int fromID, toID;
			assert(tokens[0].startsWith("#"));
			assert(tokens[1]=="is");
			assert(tokens[2].startsWith("("));
			bool ok;
			fromID = tokens[2].replace("(", "").replace(", ", "").toInt(&ok);
			assert(ok);
			toID = tokens[3].replace(")", "").toInt(&ok);
			assert(ok);
			corrs[fromID] = toID;
			std::cout<<"\t"<<corrs[fromID]<<" -> "<<toID<<std::endl;
		}
	}
	
	for (std::map<int, int>::iterator iter = corrs.begin(); iter!=corrs.end(); iter++)
	{
		double * coords1 = coordsSmf1[iter->first];
		double * coords2 = coordsSmf2[iter->second];
		R3Point p1SMF(coords1[0], coords1[1], coords1[2]);
		R3Point p2SMF(coords2[0], coords2[1], coords2[2]);
		int vOrig1 = kdTree1.FindClosest(p1SMF)->id;
		int vOrig2 = kdTree2.FindClosest(p2SMF)->id;
		outputCoarseCorrs<<vOrig1<<"\n"<<vOrig2<<"\n";
	}
	
}


void Coarse2FineByProjection(const VKString & meshCoarse1, const VKString & meshCoarse2,
							 const VKString & fullCoarseMap,
							 const VKString & meshFine1, const VKString & meshFine2, 
							 const VKString & fullFineMap)
{
	std::cout<<"Coarse2FineByProjection"<<std::endl;
	R3Mesh coarse1;
	coarse1.ReadFile(meshCoarse1.c_str());
	R3Mesh coarse2;
	coarse2.ReadFile(meshCoarse2.c_str());

	R3Mesh fine1;
	fine1.ReadFile(meshFine1.c_str());
	R3Mesh fine2;
	fine2.ReadFile(meshFine2.c_str());
	
	std::ifstream fullcoarseMapText(fullCoarseMap.c_str());
	assert(fullcoarseMapText.is_open());
	std::ofstream fullFineMapText(fullFineMap.c_str());
	assert(fullFineMapText.is_open());

//	RNArray<TriVertex*> vertexCoarse1;	
	RNArray<TriVertex*> vertexFine2;
//	for (int i=0; i<coarse1.NVertices(); i++)
//	{
//		R3Point p = coarse1.VertexPosition(coarse1.Vertex(i));
//		vertexCoarse1.InsertTail(new TriVertex(i, p.X(), p.Y(), p.Z()));
//	}
	for (int i=0; i<fine2.NVertices(); i++)
	{
		R3Point p = fine2.VertexPosition(fine2.Vertex(i));
		vertexFine2.InsertTail(new TriVertex(i, p.X(), p.Y(), p.Z()));
	}
	
	R3MeshSearchTree searchTree1(&coarse1);	
//	R3Kdtree<TriVertex*> kdTree1(vertexCoarse1, &TriVertexPositionGlobalFn, NULL);
	R3Kdtree<TriVertex*> kdTree2(vertexFine2, &TriVertexPositionGlobalFn, NULL);	
	
	std::map<int, int> coarseMap;
	for (int i=0; i<coarse1.NVertices(); i++)
	{
		int v2;
		fullcoarseMapText>>v2;
		coarseMap[i] = v2;
	}
	
	for (int j=0; j<fine1.NVertices(); j++)
	{
		// find nearest triangle in coarse1
		R3MeshIntersection closest;
		searchTree1.FindClosest(fine1.VertexPosition(fine1.Vertex(j)), closest);
		
		assert(closest.type!=R3_MESH_NULL_TYPE);
		// map vertex of triangle by coarse map
		int mapv1 = coarseMap[coarse1.VertexID(coarse1.VertexOnFace(closest.face, 0))];
		int mapv2 = coarseMap[coarse1.VertexID(coarse1.VertexOnFace(closest.face, 1))];
		int mapv3 = coarseMap[coarse1.VertexID(coarse1.VertexOnFace(closest.face, 2))];
		
		R3Point mapP1 = coarse2.VertexPosition(coarse2.Vertex(mapv1));
		R3Point mapP2 = coarse2.VertexPosition(coarse2.Vertex(mapv2));
		R3Point mapP3 = coarse2.VertexPosition(coarse2.Vertex(mapv3));
		
		R3Point avePnt((mapP1.X()+mapP2.X()+mapP3.X())/3.,
					   (mapP1.Y()+mapP2.Y()+mapP3.Y())/3.,
					   (mapP1.Z()+mapP2.Z()+mapP3.Z())/3.);
		
		
		// find nearest vertex in fine 2
		int corrV = kdTree2.FindClosest(avePnt)->id;
		if (j<100)
			std::cout<<"Map["<<j<<"]->"<<corrV<<std::endl;
		fullFineMapText<<corrV<<"\n";
	}
	fullFineMapText.close();
}
