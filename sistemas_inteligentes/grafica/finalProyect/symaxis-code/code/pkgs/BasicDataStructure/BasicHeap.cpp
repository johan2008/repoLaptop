#ifndef __MY_BASICHEAP_C
#define __MY_BASICHEAP_C
template <class T>
BasicHeap<T>::BasicHeap(void)
	:m_FlagLeastFirst(0)
{}

template <class T>
BasicHeap<T>::
BasicHeap(int flag_least_first)
	:m_FlagLeastFirst(flag_least_first)
{}

template <class T>
BasicHeap<T>::
BasicHeap(const std::vector<T*>& array, int flag_least_first)
	:m_DataArray(array), m_FlagLeastFirst(flag_least_first)
{
	// establish a map from data (T*) to the index in the tree
	for (int i=0; i<NEntries(); i++)
	{
		m_MapDataIndex[m_DataArray[i]] = i;
	}
	// make heap
	MakeHeap();
}

template <class T>
int BasicHeap<T>::
Mount(const std::vector<T*>& array, int flag_least_first)
{
	m_DataArray.resize(array.size());
	for (unsigned i=0; i<array.size(); i++)
		m_DataArray[i] = array[i];
	m_FlagLeastFirst = flag_least_first;
	for (int i=0; i<NEntries(); i++)
	{
		m_MapDataIndex[m_DataArray[i]] = i;
	}
	MakeHeap();
	return 1;
}

template <class T>
int BasicHeap<T>::
BubbleUp(int idx)
{
	while (idx > 0)
	{
		int parent = (idx - 1) / 2;
		if (!m_FlagLeastFirst && (*m_DataArray[parent]) >= (*m_DataArray[idx])) break;
		if (m_FlagLeastFirst && (*m_DataArray[parent]) <= (*m_DataArray[idx])) break;
		Swap(parent, idx);	
		idx = parent;
	}
	// return where the bubble stops
	return idx;
}

template <class T>
int BasicHeap<T>::
BubbleDown(int idx)
{
	int children[2];
	children[0] = 2*idx + 1;
	children[1] = 2*idx + 2;
	while (children[0] < NEntries())
	{
		if (!m_FlagLeastFirst)
		{
			int largest = idx;
			if ((*m_DataArray[children[0]]) > (*m_DataArray[idx]))
			{
				largest = children[0];
			}
			if (children[1] < NEntries() && 
					(*m_DataArray[children[1]]) > (*m_DataArray[largest]))
			{
				largest = children[1];
			}
			if (largest == idx) break;
			Swap(largest, idx);
			idx = largest;
		}
		else
		{
			int least = idx;
			if ((*m_DataArray[children[0]]) < (*m_DataArray[idx]))
			{
				least = children[0];
			}
			if (children[1] < NEntries() &&
					(*m_DataArray[children[1]]) < (*m_DataArray[least]))
			{
				least = children[1];
			}
			if (least == idx) break;
			Swap(least, idx);
			idx = least;
		}
		children[0] = 2*idx + 1;
		children[1] = 2*idx + 2;
	}
	return idx;
}

template <class T>
T* BasicHeap<T>::
Peek(void) const
{
	assert(!IsEmpty());
	return m_DataArray[0];
}

template <class T>
T* BasicHeap<T>::
Pop(void)
{
	assert(!IsEmpty());
	T* ret = m_DataArray[0];
	Swap(0, NEntries()-1);
	RemoveTail();
	// maintain the heap
	if (!IsEmpty())
	{
		BubbleDown(0);
	}
	return ret;
}

template <class T>
int BasicHeap<T>::
Add(T* item)
{
	// make sure the item doesn't exist
	assert(m_MapDataIndex.find(item) == m_MapDataIndex.end());
	// add to the vector
	m_DataArray.push_back(item);
	m_MapDataIndex[item] = NEntries() - 1;
	// maintain the heap
	BubbleUp(NEntries()-1);
	return 1;
}

template <class T>
int BasicHeap<T>::
Remove(T* item)
{
	assert(m_MapDataIndex.find(item) != m_MapDataIndex.end());
	int index = m_MapDataIndex[item];
	// swap item and the last item
	Swap(index, NEntries()-1);
	// remove 
	RemoveTail();
	// maintain the heap
	if (!IsEmpty())
	{
		int _idx = BubbleUp(index);
		if (_idx == index)
		{
			BubbleDown(index);
		}
	}
	return 1;
}

template <class T>
int BasicHeap<T>::
RemoveTail(void)
{
	assert(!IsEmpty());
	// remove the last one
	m_MapDataIndex.erase(m_DataArray[NEntries()-1]);
	m_DataArray.pop_back();
	return 1;
}

template <class T>
int BasicHeap<T>::
MakeHeap(void)
{
	// repeatedly bubble-up
	for (int i=1; i<NEntries(); i++)
	{
		BubbleUp(i);
	}
	return 1;
}

template <class T>
int BasicHeap<T>::
NEntries(void) const
{
	return int(m_DataArray.size());
}

template <class T>
int BasicHeap<T>::
IsEmpty(void) const
{
	return NEntries()==0? 1:0;
}

template <class T>
T* BasicHeap<T>::
Head(void) const
{
	assert(!IsEmpty());
	return m_DataArray[0];
}

template <class T>
T* BasicHeap<T>::
Tail(void) const
{
	assert(!IsEmpty());
	return m_DataArray[NEntries()-1];
}

template <class T>
T* BasicHeap<T>::
Kth(int idx) const
{
	assert(idx>=0 && idx<NEntries());
	return m_DataArray[idx];
}

template <class T>
int BasicHeap<T>::
GetIndex(T* item)
{
	assert(m_MapDataIndex.find(item) != m_MapDataIndex.end());
	return m_MapDataIndex[item];
}


template <class T>
int BasicHeap<T>::
Swap(int idx0, int idx1)
{
	if (idx0 == idx1) return 1;
	T* temp = m_DataArray[idx0];
	m_DataArray[idx0] = m_DataArray[idx1];
	m_DataArray[idx1] = temp;
	m_MapDataIndex[m_DataArray[idx0]] = idx0;
	m_MapDataIndex[m_DataArray[idx1]] = idx1;
	return 1;
}
#endif

