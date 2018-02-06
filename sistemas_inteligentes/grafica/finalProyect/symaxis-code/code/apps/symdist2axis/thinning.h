#ifndef _MY_THINNING_ALG
#define _MY_THINNING_ALG
#include <vector>
#include <iostream>
#include "R3Shapes/R3Shapes.h"
extern int print_debug;
extern int print_verbose;
extern R3MeshProperty* extern_binary_property;
using namespace std;
static void ConnectedComponentsFromBinaryLabel(R3MeshProperty* property, vector< vector<int> >& _cc)
{
	R3Mesh* mesh = property->mesh;
	int _n = mesh->NVertices();
	_cc.clear();
	vector<int> visited;
	visited.resize(_n);
	for (int i=0; i<visited.size(); i++) visited[i] = 0;
	vector<int> _queue;
	while (1)
	{
		int _exist = -1;
		for (int i=0; i<visited.size(); i++)
		{
			if (visited[i] ==0)
			{
				_exist = i;
				break;
			}
		}
		if (_exist == -1) break;
		int label = property->VertexValue(_exist);
		_queue.clear();
		_queue.push_back(_exist);
		visited[_exist] = 1;
		int _p = 0;
		while (_p<_queue.size())
		{
			R3MeshVertex* vertex = mesh->Vertex(_queue[_p]);
			for (int e=0; e<mesh->VertexValence(vertex); e++)
			{
				R3MeshVertex* _v = mesh->VertexAcrossEdge(mesh->EdgeOnVertex(vertex, e), vertex);
				int _v_id = mesh->VertexID(_v);
				if (!visited[_v_id] && (property->VertexValue(_v_id) == label))
				{
					_queue.push_back(_v_id);
					visited[_v_id] = 1;
				}
			}
			_p++;
		}
		vector<int> component(_queue);
		_cc.push_back(component);
	}
}


// check whether a vertex in on the valley
static int CheckVertexOnValley(R3MeshProperty* property, int idx, double x1, double x2, double * edge=NULL)
{
	R3Mesh* mesh=property->mesh;
	R3MeshVertex* vertex = mesh->Vertex(idx);
	int label1 = 0, label2=0;
	for (int e=0; e<mesh->VertexValence(vertex); e++)
	{
		R3MeshVertex* _v = mesh->VertexAcrossEdge(mesh->EdgeOnVertex(vertex, e), vertex);
		double value = property->VertexValue(_v);
		if (value == x1) { if (edge) *edge = x1; label1 = 1;}
		if (value == x2) { if (edge) *edge = x2; label2 = 1;}
	}
	return (label1&&label2);
}


// Truncate a raw value
static R3MeshProperty* TruncateRawValue(R3MeshProperty* raw, double thres)
{
	R3MeshProperty* trunc = new R3MeshProperty(raw->mesh, "binary");
	int _n = raw->mesh->NVertices();
	int counter = 0;
	for (int i=0; i<_n; i++)
	{
		double value = ((raw->VertexValue(i)<thres)?1.0:0.0);
		if (value == 0.0) counter ++;
		trunc->SetVertexValue(i, value);
	}
//	cout<<counter<<" vertices are truncated!"<<endl;
	return trunc;
}


// Improve binary value
///////// 0-component and 1,2-component must be connected
// 0, 1, 2
static R3MeshProperty* ImproveBinaryValue(R3MeshProperty* binary)
{
	vector< vector<int> > connected_components;
	R3MeshProperty* property = new R3MeshProperty(*binary);
	property->SetName("ImprovedBinary");
	ConnectedComponentsFromBinaryLabel(property, connected_components);

//	cout<<"Before Improve: "<<connected_components.size()<<endl;
	// find two 0-components with largest size
	int b0, b1, s0=0, s1=0;
	for (int i=0; i<connected_components.size(); i++)
	{	
		if (property->VertexValue(connected_components[i][0]) == 1.0) continue;
		if (connected_components[i].size()>s0)
		{
			s1 = s0;
			b1 = b0;
			s0 = connected_components[i].size();
			b0 = i;
		}
		else if (connected_components[i].size()>s1)
		{
			s1 = connected_components[i].size();
			b1 = i;
		}
	}
	for (int i=0; i<connected_components.size(); i++)
	{
		if (property->VertexValue(connected_components[i][0]) == 1.0) continue;
		if (i == b0 || i == b1) continue;
		for (int j=0; j<connected_components[i].size(); j++)
		{
			property->SetVertexValue(connected_components[i][j], 1.0);
		}
	}
	ConnectedComponentsFromBinaryLabel(property, connected_components);

	// find the largest 1-components
	s0 = 0;
	for (int i=0; i<connected_components.size(); i++)
	{
		if (property->VertexValue(connected_components[i][0]) == 0.0) continue;
		if (connected_components[i].size()>s0)
		{
			s0 = connected_components[i].size();
			b0 = i;
		}
	}
	for (int i=0; i<connected_components.size(); i++)
	{
		if (property->VertexValue(connected_components[i][0]) == 0.0) continue;
		if (i == b0) continue;		
		for (int j=0; j<connected_components[i].size(); j++)
		{
			property->SetVertexValue(connected_components[i][j], 0.0);
		}
	}
	ConnectedComponentsFromBinaryLabel(property, connected_components);
//	cout<<"After Improve: "<<connected_components.size()<<endl;

	for (int i=0; i<connected_components.size(); i++)
	{
		if (property->VertexValue(connected_components[i][0]) == 1.0) continue;
		for (int j=0; j<connected_components[i].size(); j++)
		{
			property->SetVertexValue(connected_components[i][j], 2.0);
		}
		break;
	}
		
	return property;	
}

struct CurveData
{
	R3MeshVertex* Vertex;
	double Value;
	CurveData** heappointer;
};

static R3MeshProperty* CurveExtract(R3MeshProperty* raw, R3MeshProperty* binary, double curve = 1.0)
{
	R3Mesh* Mesh = raw->Mesh();
	R3MeshProperty* curve_prop = new R3MeshProperty(Mesh, "curve");
	for (int i=0; i<Mesh->NVertices(); i++)
		curve_prop->SetVertexValue(i, 0.0);
	CurveData* vertex_data = new CurveData[Mesh->NVertices()];
	for (int i=0; i<Mesh->NVertices(); i++)
	{
		vertex_data[i].Vertex = Mesh->Vertex(i);
		vertex_data[i].Value =  raw->VertexValue(i);
		vertex_data[i].heappointer = NULL;
	}
	CurveData tmp;
	RNHeap<CurveData*> heap(&tmp, &(tmp.Value), &(tmp.heappointer), FALSE);
    // Find all vertices on the edge
	for (int i=0; i<Mesh->NVertices(); i++)
	{
		R3MeshVertex* vertex = vertex_data[i].Vertex;
		double binary_value = binary->VertexValue(i);//vertex_data[i].BinaryValue;
		if (binary_value != curve) continue;
		for (int e=0; e<Mesh->VertexValence(vertex); e++)
		{
			R3MeshEdge* edge = Mesh->EdgeOnVertex(vertex, e);
			R3MeshVertex* _v = Mesh->VertexAcrossEdge(edge, vertex);
			double _b = binary->VertexValue(_v);
			if (_b != curve)
			{
				// on the edge
				heap.Push(&vertex_data[i]);
				break;
			}
		}
	}

	while (!heap.IsEmpty())
	{
		CurveData* top = heap.Pop();
		// to tell if top is on the valley
		int idx = Mesh->VertexID(top->Vertex);
		double side;
		int valley = CheckVertexOnValley(binary, idx, 0.0, 2.0, &side);
		if (valley)
			curve_prop->SetVertexValue(idx, 1.0);
		else
		{
			// erosion, add its neighbors
			binary->SetVertexValue(idx, side);
			for (int e=0; e<Mesh->VertexValence(top->Vertex); e++)
			{
				R3MeshEdge* edge = Mesh->EdgeOnVertex(top->Vertex, e);
				R3MeshVertex* _v = Mesh->VertexAcrossEdge(edge, top->Vertex);
				CurveData * _data = &vertex_data[Mesh->VertexID(_v)];
				if (binary->VertexValue(_v) == curve && curve_prop->VertexValue(_v) == 0.0 && _data->heappointer == NULL)
				{
					// in "Curve" Region && not decided && not in the heap
					heap.Push(_data);
				}
			}
		}
   

	}
	delete[] vertex_data;
	return curve_prop;
}


// trace a curve
static RNArray<R3MeshVertex*>* TraceCurve(R3MeshProperty* property, int val)
{
	R3Mesh* mesh = property->mesh;
	R3MeshVertex* seed = NULL;
	RNArray<R3MeshVertex*> candidates;
	for (int i=0; i<mesh->NVertices(); i++)
	{
		if (property->VertexValue(i)==val)
		{
			candidates.Insert(mesh->Vertex(i));
			if (seed == NULL)
				seed = mesh->Vertex(i);
		}
	}
	RNArray<R3MeshVertex*>* array = new RNArray<R3MeshVertex*>;
	if (seed == NULL) return array;
	array->Insert(seed);
	R3MeshVertex* prev = NULL;
	R3MeshVertex* cur = seed;
	while (1)
	{
		int suc = 0;
		for (int e=0; e<mesh->VertexValence(cur); e++)
		{
			R3MeshVertex * _v = mesh->VertexAcrossEdge(mesh->EdgeOnVertex(cur, e), cur);
			if (property->VertexValue(_v)==val && _v!=prev)
			{
				prev = cur;
				cur = _v;
				suc = 1;
				break;
			}
		}
		if (suc == 0)
		{
			// What if there is not adjacent point with val
			RNLength* dists = mesh->DijkstraDistances(cur);
			RNLength min_dist = FLT_MAX;
			R3MeshVertex* closest_vertex = NULL;
			for (int i=0; i<candidates.NEntries(); i++)
			{
				R3MeshVertex* _v = candidates.Kth(i);
				if (dists[mesh->VertexID(_v)] < min_dist)
				{
					min_dist = dists[mesh->VertexID(_v)];
					closest_vertex = _v;
				}
			}
			if (!closest_vertex)
			{
				cout<<"Unable to continue at vertex "<<mesh->VertexID(cur);
				exit(-1);
			}
			delete[] dists;
			RNArray<R3MeshEdge*> edges;
			mesh->DijkstraDistance(closest_vertex, cur, &edges);
			for (int i=0; i<edges.NEntries(); i++)
			{
				R3MeshEdge* edge = edges.Kth(i);
				cur = mesh->VertexAcrossEdge(edge, cur);
				array->Insert(cur);
				if (candidates.FindEntry(cur))
					candidates.Remove(cur);
			}
			continue;
		}
		if (cur == seed)
			break;
		array->Insert(cur);
		if (candidates.FindEntry(cur))
			candidates.Remove(cur);
//		if (print_debug)
//			cout<<"NCandiates : "<<candidates.NEntries()<<endl;
//		getchar();
	}
	cout<<"NCandidates remaining: "<<candidates.NEntries()<<endl;
	return array;
}

static RNArray<R3MeshVertex*>* Thinning(R3MeshProperty* _property, int flagifmin = 1)
{
	R3MeshProperty* property = new R3MeshProperty(*_property);
	R3Mesh* mesh = property->Mesh();

	if (!flagifmin)
	{
		double max = property->Maximum();
		for (int i=0; i<mesh->NVertices(); i++)
		{
			double val = property->VertexValue(i);
			property->SetVertexValue(i, max - val);
		}
		property->minimum = RN_UNKNOWN;
		double min = property->Minimum();
		for (int i=0; i<mesh->NVertices(); i++)
		{
			double val = property->VertexValue(i);
			property->SetVertexValue(i, val - min);
		}
		property->minimum = RN_UNKNOWN;
		property->maximum = RN_UNKNOWN;
	}
//	cout<<property->Maximum()<<" : "<<property->Minimum()<<endl;
	R3MeshProperty* BinaryValue = NULL;
	if (!extern_binary_property)
	{
		double BinaryThreshold = property->Mean();
		BinaryValue = TruncateRawValue(property, BinaryThreshold);
	}
	else
		BinaryValue = extern_binary_property;
	

	if (print_debug)
	{
		cout<<"Write binaryvalue to file "<<"binaryvalue.val"<<endl;
		BinaryValue->Write("binaryvalue.val");
	}
//	cout<<"Truncate Ok! ";
	R3MeshProperty* ImprovedBinaryValue = ImproveBinaryValue(BinaryValue);
	if (print_debug)
	{
		cout<<"Write improvedbinaryvalue to "<<"impbinaryvalue.val"<<endl;
		ImprovedBinaryValue->Write("impbinaryvalue.val");
	}
//	cout<<"Improve Ok! ";
	R3MeshProperty* CurveValue = CurveExtract(property, ImprovedBinaryValue, 1.0);
	if (print_debug)
	{
		cout<<"Write Curve value to "<<"curvevalue.val"<<endl;
		CurveValue->Write("curvevalue.val");
	}
//	cout<<"Extract Ok! ";
	delete property;
	delete BinaryValue;
	delete ImprovedBinaryValue;
	RNArray<R3MeshVertex*>* array = TraceCurve(CurveValue, 1.0);
	if (array->NEntries() == 0)
	{
		cout<<"Length is 0"<<endl;
	}
	delete CurveValue;
	return array;
}


#endif
