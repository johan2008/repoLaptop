/***************************************************************************
 *   Copyright (C) 2008 by Vladimir Kim, Princeton University			   *
 *	 NOTE: PARTS OF THIS CODE COME FROM GAPS LIBRARY BY Tom Funkhouser	   *
 *   Author: Vladimir Kim (kim_vladimir@yahoo.com, vk@princeton.edu)       *
 *                                  ****                                   *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef __MESH_PROCESSOR_H
#define __MESH_PROCESSOR_H


#include <gaps.h>
#include <vector>
#include <iostream>
#include <assert.h>

class MeshProcessor
{
	public:
		MeshProcessor(R3Mesh * mesh);
		void Smooth(int numIterations);
		void SingleConnectedComponent();
		void CopyAndForget();
		void DelaunayTriangByEdgeFlip();
		
		void ConsistentFaceOrientation(const char * filename);
		
	
		R3Mesh * m_pMesh;
		
		int NumConnectedComponents(std::vector<RNArray<R3MeshFace*> > & connectedComps);
		static int NumConnectedComponents(R3Mesh * mesh);
		static int NumIsolatedVertices(R3Mesh * mesh);

	protected:
		R3Vector GetFaceNorm(int faceID, double * vertexCoords, int * faceIDs, int * flips);
		bool Delaunay(int edgeID);
		void EdgesOnAdjacentFaces(int edgeID, std::vector<int> & adjcentEdges);
		double AngleBetweenEdges(int e1, int e2);
		bool CanSwapEdge(int edgeID);
};


#endif

