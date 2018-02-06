#include "MeshProcessor.h"
#include <map>
#include <fstream>

R3Point vertexPositionFunction(R3MeshVertex * v, void * mesh)
{
	return ((R3Mesh*)mesh)->VertexPosition(v);
}

MeshProcessor::MeshProcessor(R3Mesh * mesh)
{
	m_pMesh = mesh;
}

void MeshProcessor::CopyAndForget()
{
	R3Mesh * newMesh = new R3Mesh();
	for (int i=0; i<m_pMesh->NVertices(); i++)
		newMesh->CreateVertex(m_pMesh->VertexPosition(m_pMesh->Vertex(i)));
	
	for (int i=0; i<m_pMesh->NFaces(); i++)
	{
		int v1 = m_pMesh->VertexID(m_pMesh->VertexOnFace(m_pMesh->Face(i), 0));
		int v2 = m_pMesh->VertexID(m_pMesh->VertexOnFace(m_pMesh->Face(i), 1));
		int v3 = m_pMesh->VertexID(m_pMesh->VertexOnFace(m_pMesh->Face(i), 2));		
		newMesh->CreateFace(newMesh->Vertex(v1), newMesh->Vertex(v2), newMesh->Vertex(v3));
	}
	m_pMesh = newMesh;
}

void MeshProcessor::Smooth(int numIterations)
{
	double * xCoords = new double[m_pMesh->NVertices()];
	double * yCoords = new double[m_pMesh->NVertices()];
	double * zCoords = new double[m_pMesh->NVertices()];	
	double weight;
	for (int i=0; i<numIterations; i++)
	{
		for (int v=0; v<m_pMesh->NVertices(); v++)	// mean shift
		{
			double norm=0;			
			R3MeshVertex * vp = m_pMesh->Vertex(v);
			weight = m_pMesh->VertexValence(vp);			
			
			xCoords[v] = m_pMesh->VertexPosition(vp).X() * weight;
			yCoords[v] = m_pMesh->VertexPosition(vp).Y() * weight;
			zCoords[v] = m_pMesh->VertexPosition(vp).Z() * weight;			
			norm+=weight;
			for (int nb=0; nb<m_pMesh->VertexValence(vp); nb++)
			{
				weight = 1;
				R3Point pos = m_pMesh->VertexPosition(m_pMesh->VertexAcrossEdge(m_pMesh->EdgeOnVertex(vp, nb), vp));
				xCoords[v] += pos.X() * weight;
				yCoords[v] += pos.Y() * weight;
				zCoords[v] += pos.Z() * weight;				
				norm+=weight;				
			}
			xCoords[v] /= norm;
			yCoords[v] /= norm;
			zCoords[v] /= norm;			
		}
		
		for (int v=0; v<m_pMesh->NVertices(); v++)	// re-assign values
			m_pMesh->SetVertexPosition(m_pMesh->Vertex(v), R3Point(xCoords[v], yCoords[v], zCoords[v]));
	}
	
	delete [] xCoords;
	delete [] yCoords;
	delete [] zCoords;	
}

void MeshProcessor::SingleConnectedComponent()
{
	assert(false);
	std::vector<RNArray<R3MeshFace*> > connectedComps;
	NumConnectedComponents(connectedComps);
	
	std::cout<<"\tConnected Components Found: "<<connectedComps.size()<<std::endl;
	assert(connectedComps.size()==1);
	while (connectedComps.size()>1)
	{
		int largestClusterID=0;
		for (int i=1; i<(int)connectedComps.size(); i++)
			if (connectedComps[i].NEntries() > connectedComps[largestClusterID].NEntries())
				largestClusterID = i;
		
		RNArray<R3MeshVertex *> vts;
		std::map<R3MeshVertex *, int> addedVertex;
		for (int i=0; i<connectedComps[largestClusterID].NEntries(); i++)
		{
			R3MeshVertex * v1 = m_pMesh->VertexOnFace(connectedComps[largestClusterID][i], 0);
			R3MeshVertex * v2 = m_pMesh->VertexOnFace(connectedComps[largestClusterID][i], 1);
			R3MeshVertex * v3 = m_pMesh->VertexOnFace(connectedComps[largestClusterID][i], 2);			
			if (addedVertex.find(v1)==addedVertex.end())
				{		addedVertex[v1] = vts.NEntries();			vts.Insert(v1);			}
			if (addedVertex.find(v2)==addedVertex.end())
				{		addedVertex[v2] = vts.NEntries();			vts.Insert(v2);			}
			if (addedVertex.find(v3)==addedVertex.end())
				{		addedVertex[v3] = vts.NEntries();			vts.Insert(v3);			}
		}
		R3Kdtree<R3MeshVertex*> verticesTree(vts, &vertexPositionFunction, m_pMesh);
		
		addedVertex.clear();
		for (int i=0; i<(int)connectedComps.size(); i++)
		{
			if (i==largestClusterID)
				continue;
			
			std::cout<<"\tFaces in a small CC: "<<connectedComps[i].NEntries()<<std::flush;
			std::cout<<" Largest="<<connectedComps[largestClusterID].NEntries()<<std::endl;
			// TODO: REMOVE CC
		}
	}
	connectedComps.clear();
	assert(NumConnectedComponents(connectedComps)==1);
	int vertexCollapsed=0;
	for (int i=0; i<m_pMesh->NVertices(); i++)
	{
		if (m_pMesh->VertexValence(m_pMesh->Vertex(i))==0)
		{
			m_pMesh->DeleteVertex(m_pMesh->Vertex(i--));
			vertexCollapsed++;
		}
	}
	
	std::cout<<"\tVertices Removed: "<<vertexCollapsed<<std::endl;
}

int MeshProcessor::NumConnectedComponents(std::vector<RNArray<R3MeshFace*> > & connectedComps)
{
	std::map<R3MeshFace*, int> faceToCluster;
	for (int i=0; i<m_pMesh->NFaces(); i++)
	{
		if (faceToCluster.find(m_pMesh->Face(i))==faceToCluster.end())
		{
			int currCluster = (int)connectedComps.size();
			connectedComps.push_back(RNArray<R3MeshFace*>());
			m_pMesh->FindConnectedFaces(m_pMesh->Face(i), connectedComps[currCluster]);
			for (int j=0; j<connectedComps[currCluster].NEntries(); j++)
				faceToCluster[connectedComps[currCluster][j]] = currCluster;
			faceToCluster[m_pMesh->Face(i)] = (int)connectedComps.size()-1;
		}
	}
	return (int)connectedComps.size();
}

void MeshProcessor::ConsistentFaceOrientation(const char * filename)
{
	assert(false);
	// only off format is supported. 
	std::cout<<"Loading OFF file"<<std::endl;
	assert(strcmp(&(filename[strlen(filename)-4]), ".off")==0);
	std::ifstream textStream(filename);
	assert(textStream.is_open());
	
	std::string tempStr;
	textStream>>tempStr;		assert(strcmp(tempStr.c_str(), "OFF")==0);
	int vertexCount, faceCount, edgeCount;
	
	textStream>>vertexCount>>faceCount>>edgeCount;
	assert(vertexCount>=0 && faceCount>=0 && edgeCount>=0);
	
	double * vertexCoords = new double[3*vertexCount];
	int * faceIds = new int[3*faceCount];
	int * flip = new int[faceCount];	// 0 - unexplored, 1 - keep, 2 - flip
	std::vector<int> * vertexToFace = new std::vector<int>[vertexCount];
	
	for (int i=0; i<vertexCount * 3; i++)
		textStream>>vertexCoords[i];
	
	for (int i=0; i<faceCount; i++)
	{
		int numVert;
		textStream>>numVert;
		assert(numVert==3);
		textStream>>faceIds[3*i+0];
		textStream>>faceIds[3*i+1];
		textStream>>faceIds[3*i+2];		
		vertexToFace[faceIds[3*i+0]].push_back(i);
		vertexToFace[faceIds[3*i+1]].push_back(i);
		vertexToFace[faceIds[3*i+2]].push_back(i);		
		flip[i] = 0;	// unexplored
	}
	
	std::cout<<"Finding necessary flips"<<std::endl;
	std::map<int, bool> visited;
	int numFacesFlipped=0;
	for (int i=0; i<faceCount; i++)
	{
		std::vector<int> facesToVisit;
		if (visited.find(i)!=visited.end())
			continue;
		std::cout<<"New seed: "<<i<<std::endl;
		facesToVisit.push_back(i);
		flip[i] = 1;
		visited[i] = true;
		while (facesToVisit.size()>0)
		{
			int faceID = facesToVisit[0];
			assert(flip[faceID]!=0);
			R3Vector currNorm = GetFaceNorm(faceID, vertexCoords, faceIds, flip);
			std::map<int, bool> adjacentFaces;
			for (int vertID=0; vertID<3; vertID++)
			{
				int v = faceIds[faceID*3 + vertID];
				for (int adjFaceID = 0; adjFaceID < (int)vertexToFace[v].size(); adjFaceID++)
					adjacentFaces[vertexToFace[v][adjFaceID]] = true;
			}
			
			for (std::map<int, bool>::iterator iter=adjacentFaces.begin(); 
				 iter!=adjacentFaces.end(); iter++)
			{
				int adjFaceID = iter->first;
				R3Vector adjNorm = GetFaceNorm(adjFaceID, vertexCoords, faceIds, flip);
				if (adjNorm.Dot(currNorm) < 0 && flip[adjFaceID]==0)
				{
	//				assert(flip[adjFaceID]==0);	// verify that face is not set yet
					flip[adjFaceID] = 2;		// flip
					numFacesFlipped++;
				}
				else if (flip[adjFaceID]==0)
				{
					flip[adjFaceID] = 1;	// keep as is (unless already set to something)
					numFacesFlipped++;
				}
								
				if (visited.find(adjFaceID)==visited.end())
				{
					visited[adjFaceID] = true;
					facesToVisit.push_back(adjFaceID);
				}
			}			
			
			facesToVisit.erase(facesToVisit.begin());
		}
	}
	
	std::cout<<"Finding minority to flip"<<std::endl;
	int requiredFlips = 0;
	for (int i=0; i<faceCount; i++)
	{
		assert(flip[i]!=0);
		if (flip[i]==2)
			requiredFlips++;
	}
	if (requiredFlips > faceCount / 2)
	{
		for (int i=0; i<faceCount; i++)
		{
			if (flip[i]==2)
				flip[i]=1;
			else
				flip[i]=2;
		}
	}
	
	std::cout<<"Creating new R3Mesh need to flip: min("<<requiredFlips<<", "<<(faceCount-requiredFlips)<<")"<<std::endl;	
	m_pMesh = new R3Mesh();
	for (int i=0; i<vertexCount; i++)	// create vertices
		m_pMesh->CreateVertex(R3Point(vertexCoords[3*i+0], vertexCoords[3*i+1], vertexCoords[3*i+2]));
	
	for (int i=0; i<faceCount; i++)		// create faces
	{
		R3MeshVertex * v1 = m_pMesh->Vertex(faceIds[3*i+0]);
		R3MeshVertex * v2 = m_pMesh->Vertex(faceIds[3*i+1]);
		R3MeshVertex * v3 = m_pMesh->Vertex(faceIds[3*i+2]);		
		
		if (flip[i]==1)	// keep as is
			m_pMesh->CreateFace(v1, v2, v3);
		else if (flip[i]==2)
			m_pMesh->CreateFace(v3, v2, v1);
		else
			assert(false);
	}
	
	delete [] vertexCoords;
	delete [] faceIds;
	delete [] flip;
	delete [] vertexToFace;
}

R3Vector MeshProcessor::GetFaceNorm(int faceID, double * vertexCoords, int * faceIDs, int * flips)
{
	int v1 = faceIDs[faceID*3+0];
	int v2 = faceIDs[faceID*3+1];
	int v3 = faceIDs[faceID*3+2];	
	
	double x1 = vertexCoords[3*v1+0];
	double y1 = vertexCoords[3*v1+1];
	double z1 = vertexCoords[3*v1+2];	

	double x2 = vertexCoords[3*v2+0];
	double y2 = vertexCoords[3*v2+1];
	double z2 = vertexCoords[3*v2+2];	

	double x3 = vertexCoords[3*v3+0];
	double y3 = vertexCoords[3*v3+1];
	double z3 = vertexCoords[3*v3+2];	
	
	R3Vector v12(x2-x1, y2-y1, z2-z1);
	R3Vector v13(x1-x3, y1-y3, z1-z3);	
//	R3Vector v13(x3-x1, y3-y1, z3-z1);
	
	v12.Normalize();
	v13.Normalize();
	
//	std::cout<<"["<<v12.X()<<", "<<v12.Y()<<", "<<v12.Z()<<"] x ["<<v13.X()<<", "<<v13.Y()<<", "<<v13.Z()<<"]  ";
	v12.Cross(v13);
//	std::cout<<" = "<<"["<<v12.X()<<", "<<v12.Y()<<", "<<v12.Z()<<"]"<<std::endl;
	R3Vector norm = v12;
	norm.Normalize();
	
	if (flips[faceID]==0 || flips[faceID]==1)
		return norm;
	else
	{
		assert(flips[faceID]==2);
		return -norm;
	}
	
}

void MeshProcessor::DelaunayTriangByEdgeFlip()
{
	assert(m_pMesh!=NULL);
	// implementation of An Algorithm for the Construction of Intrinsic Delaunay Triangulations 
	// with Applications to Digital Geometry Processing. Matthew Fisher et al. 
	
	// stack of edges / marks of edges
	std::vector<int> edgeStack;
	bool * edgeMarks = new bool[m_pMesh->NEdges()];
	
	int nonDelaunay=0;
	for (int i=0; i<m_pMesh->NEdges(); i++)
	{
		edgeStack.push_back(i);
		edgeMarks[i] = true;
		if (!Delaunay(i))
			nonDelaunay++;
	}
	int numFlips=0;
	int MAX_ITERATIONS = 10000000;
	std::cout<<"\tInitial Non-Delaunay Edges = "<<nonDelaunay<<" / "<<m_pMesh->NEdges()<<std::endl;
	while (edgeStack.size()>0)
	{
		MAX_ITERATIONS--;
		if (MAX_ITERATIONS<0)
			break;
		int edgeID = edgeStack[0];
		edgeStack.erase(edgeStack.begin());
		edgeMarks[edgeID] = false;
		if (!Delaunay(edgeID) && CanSwapEdge(edgeID))
		{
			numFlips++;
			if (numFlips%10000==0)
			{
				std::cout<<"\r"<<std::flush;
				int MAX_BARS = 50;
				int CURR_BAR = (m_pMesh->NEdges() - edgeStack.size()) * MAX_BARS / m_pMesh->NEdges();
				std::cout<<"["<<std::flush;
				for (int i=0; i<MAX_BARS; i++)
					if (i<CURR_BAR)
						std::cout<<"."<<std::flush;
					else
						std::cout<<" "<<std::flush;
				std::cout<<"] "<<edgeStack.size()<<std::flush;
			}
			std::vector<int> adjcentEdges;			
			EdgesOnAdjacentFaces(edgeID,  adjcentEdges);
			
			m_pMesh->SwapEdge(m_pMesh->Edge(edgeID));

			for (int i=0; i<(int)adjcentEdges.size(); i++)
			{
				if (!edgeMarks[adjcentEdges[i]])
				{
					edgeMarks[adjcentEdges[i]] = true;
					edgeStack.push_back(adjcentEdges[i]);
				}
			}
		}
	}
	std::cout<<"\r"<<std::flush;	
	std::cout<<"\tDelaunay Triangulation Flipped Edges = "<<numFlips<<std::endl;

	nonDelaunay=0;
//	std::cout<<"\tNon-delaunay edges: "<<std::flush;
	for(int i=0; i<m_pMesh->NEdges(); i++)
	{
		if (!Delaunay(i))
		{
//			std::cout<<" "<<i<<std::flush;
			nonDelaunay++;
		}
	}
//	std::cout<<" | TOTAL="<<nonDelaunay<<std::endl;
	std::cout<<"\tRemaining Non-Delaunay Edges = "<<nonDelaunay<<std::endl;
}

bool MeshProcessor::Delaunay(int edgeID)
{
//	std::cout<<"\tDelaunay "<<edgeID<<" ? "<<std::endl;
	std::vector<int> adjcentEdges;
	EdgesOnAdjacentFaces(edgeID, adjcentEdges);
	
	if (adjcentEdges.size()<4)
		return true;
	
	double alpha = AngleBetweenEdges(adjcentEdges[0], adjcentEdges[1]);
	double beta = AngleBetweenEdges(adjcentEdges[2], adjcentEdges[3]);
	
//	std::cout<<"\t\tAdjacent Edges: ["<<adjcentEdges[0]<<", "<<adjcentEdges[1]<<", ";
//	std::cout<<adjcentEdges[2]<<", "<<adjcentEdges[3]<<"] "<<std::endl;
//	std::cout<<"\t\tAlpha="<<acos(e11.Dot(e12))<<" + "<<acos(e21.Dot(e22))<<" Del=";
//	std::cout<<((alpha + beta) <= 3.14159266)<<std::endl;
	
	return (alpha + beta) <= 3.14159266;
}

void MeshProcessor::EdgesOnAdjacentFaces(int edgeID, std::vector<int> & adjcentEdges)
{
	R3MeshEdge * edge = m_pMesh->Edge(edgeID);
	R3MeshFace * face1 = m_pMesh->FaceOnEdge(edge, 0);
	R3MeshFace * face2 = m_pMesh->FaceOnEdge(edge, 1);	
	if (face1==NULL && face2==NULL)
	{
		std::cout<<"[WARNING] Null edge in Delaunay Triangulation procedure"<<std::endl;
		return;
	}
	if (face1==NULL)
	{
		face1 = face2;
		face2 = NULL;
	}
	assert(face1!=NULL);
	
	R3MeshEdge * e1 = m_pMesh->EdgeOnFace(face1, 0);
	R3MeshEdge * e2 = m_pMesh->EdgeOnFace(face1, 1);
	R3MeshEdge * e3 = m_pMesh->EdgeOnFace(face1, 2);	
	if (e1==edge)
		e1 = e3;
	else if (e2==edge)
		e2 = e3;
	else
		assert(e3==edge);
	
	adjcentEdges.push_back(m_pMesh->EdgeID(e1));
	adjcentEdges.push_back(m_pMesh->EdgeID(e2));	
	
	if (face2!=NULL)
	{
		e1 = m_pMesh->EdgeOnFace(face2, 0);
		e2 = m_pMesh->EdgeOnFace(face2, 1);
		e3 = m_pMesh->EdgeOnFace(face2, 2);
		
		if (e1==edge)
			e1 = e3;
		else if (e2==edge)
			e2 = e3;
		else
			assert(e3==edge);
		
		adjcentEdges.push_back(m_pMesh->EdgeID(e1));
		adjcentEdges.push_back(m_pMesh->EdgeID(e2));	
	}
}

double MeshProcessor::AngleBetweenEdges(int e1ID, int e2ID)
{
	R3MeshEdge * e1 = m_pMesh->Edge(e1ID);
	R3MeshEdge * e2 = m_pMesh->Edge(e2ID);	
	R3Vector v1, v2;

	R3MeshVertex * v11 = m_pMesh->VertexOnEdge(e1, 0);
	R3MeshVertex * v12 = m_pMesh->VertexOnEdge(e1, 1);
	R3MeshVertex * v21 = m_pMesh->VertexOnEdge(e2, 0);
	R3MeshVertex * v22 = m_pMesh->VertexOnEdge(e2, 1);	
	v1 = m_pMesh->VertexPosition(v12)-m_pMesh->VertexPosition(v11);
	v2 = m_pMesh->VertexPosition(v22)-m_pMesh->VertexPosition(v21);		
	
	if (v11==v22)
		v2 = -1 * v2;
	else if (v12==v21)
		v1 = -1 * v1;
	else if (v12==v22)
	{
		v1 = -1 * v1;
		v2 = -1 * v2;
	}
	else
		assert(v11==v21);
	v1.Normalize();
	v2.Normalize();
	return acos(v1.Dot(v2));
}

bool MeshProcessor::CanSwapEdge(int edgeID)
{
	// verify if already in a degenrate case
	R3MeshEdge * edge = m_pMesh->Edge(edgeID);
	R3MeshVertex *v[2];
	for (int i=0; i<2; i++)
		v[i] = m_pMesh->VertexOnEdge(edge, i);
	
	// Get faces on edge
	R3MeshFace *f[2];
	for (int i=0; i<2; i++)
	{
		f[i] = m_pMesh->FaceOnEdge(edge, i);
		if (!f[i])
			return false;
	}

	// Get edges and vertices across faces
	R3MeshEdge *ef[2][2] = { { NULL, NULL }, { NULL, NULL } };
	R3MeshVertex *vf[2] = { NULL, NULL };
	ef[0][0] = m_pMesh->EdgeAcrossVertex(v[0], edge, f[0]);
	ef[0][1] = m_pMesh->EdgeAcrossVertex(v[1], edge, f[0]);
	vf[0] = m_pMesh->VertexAcrossEdge(ef[0][1], v[1]);
	assert(vf[0] == m_pMesh->VertexAcrossEdge(ef[0][0], v[0]));
	ef[1][0] = m_pMesh->EdgeAcrossVertex(v[0], edge, f[1]);
	ef[1][1] = m_pMesh->EdgeAcrossVertex(v[1], edge, f[1]);
	vf[1] = m_pMesh->VertexAcrossEdge(ef[1][1], v[1]);
	assert(vf[1] == m_pMesh->VertexAcrossEdge(ef[1][0], v[0]));
	if (m_pMesh->EdgeBetweenVertices(vf[0], vf[1])) 
		return false;
	
	// check that edge swap does not create non-manifold
	// Checks if there is another vertex other than v[1] adjacent to: 
	//		v[0] (old edge), vf[0] (new edge), vf[1] (new edge). 
	for (int vid=0; vid < 2; vid++)
	{
		for (int i=0; i<m_pMesh->VertexValence(v[vid]); i++)
		{
			R3MeshEdge * adjE = m_pMesh->EdgeOnVertex(v[vid], i);
			assert(m_pMesh->IsVertexOnEdge(v[vid], adjE));
			R3MeshVertex * adjV = m_pMesh->VertexAcrossEdge(adjE, v[vid]);
			assert(adjV!=NULL);
			
			if (adjV==v[1 - vid])
				continue;
			else if (m_pMesh->EdgeBetweenVertices(adjV, vf[0]) 
					 && m_pMesh->EdgeBetweenVertices(adjV, vf[1]))
				return false;
		}
		
	}
	
	return true;
}

int MeshProcessor::NumIsolatedVertices(R3Mesh * mesh)
{
	int numIsolated=0;
	for (int i=0; i<mesh->NVertices(); i++)
	{
		if (mesh->VertexValence(mesh->Vertex(i))==0)
			numIsolated++;
	}
	return numIsolated;
}

int MeshProcessor::NumConnectedComponents(R3Mesh *mesh)
{
	// Get fresh mark value
	R3mesh_mark++;
	
	// Iterate finding connected components
	int count = 0;
	int start = 0;
	for (;;) {
		// Find an unmarked face
		R3MeshFace *seed = NULL;
		for (; start < mesh->NFaces(); start++) {
			R3MeshFace *face = mesh->Face(start);
			if (mesh->FaceMark(face) != R3mesh_mark) {
				seed = face;
				break;
			}
		}
		
		// Check if found a new component
		if (!seed) break;
		else count++;
		
		// Mark connected component 
		RNArray<R3MeshFace *> stack;
		stack.InsertTail(seed);
		while (!stack.IsEmpty()) {
			R3MeshFace *face = stack.Tail();
			stack.RemoveTail();
			mesh->SetFaceMark(face, R3mesh_mark);
			for (int i = 0; i < 3; i++) {
				R3MeshFace *neighbor = mesh->FaceOnFace(face, i);
				if ((neighbor) && (mesh->FaceMark(neighbor) != R3mesh_mark)) {
					stack.InsertTail(neighbor);
				}
			}
		}
	}
	
	// Return number of connected components
	return count;
}


