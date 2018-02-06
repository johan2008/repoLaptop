#ifndef __MY_BASICHEAP_H
#define __MY_BASICHEAP_H
#include <vector>
#include <map>
/**
 * Data structure of (priority) heap
 */
template <class T>
class BasicHeap
{
public:
	// initialize
	BasicHeap(void);
	BasicHeap(int flag_least_first);
	BasicHeap(const std::vector<T*>& array, int flag_least_first = 0);
	// mount a vector<T*> and make heap
	int Mount(const std::vector<T*>& array, int flag_least_first = 0);
	// maintain the property of heap
	int BubbleUp(int idx);
	int BubbleDown(int idx);
	// operations
	T* Peek() const;
	T* Pop();
	int Add(T* item);
	int Remove(T* item);
	int RemoveTail(void);
	int MakeHeap(void);
	// statistics of heap
	int NEntries(void) const;
	int IsEmpty(void) const;
	// for convenience
	T* Head(void) const;
	T* Tail(void) const;
	// get item
	T* Kth(int idx) const;
	// get index
	int GetIndex(T* item);
private:
	std::vector<T*> m_DataArray;
	std::map<T*, int> m_MapDataIndex;
	int m_FlagLeastFirst;
private:
	// Utility
	int Swap(int idx0, int idx1);
};

#include "BasicHeap.cpp"

#endif
