#include <assert.h>

#ifndef __SORTABLE_ENTRY_H
#define __SORTABLE_ENTRY_H

class Sortable
{
	public:
		Sortable(double val, void * ptr, int id = -1)
		{
			this->value = val;
			this->id = id;
			this->ptr = ptr;
		}
		double value;
		void * ptr;
		int id;
		
		int operator<(const Sortable &B) const
		{
			if (value < B.value)
				return 1;
			else
				return 0;
		}
};

#endif
