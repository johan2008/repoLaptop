#ifndef __ANALYSIS_STATS_H
#define __ANALYSIS_STATS_H

#include <fstream>
#include <sstream>

#include "gaps.h"

#include "TimeProfiler.h"


using namespace std;

class AnalysisStats
{
	public:
		AnalysisStats();
		
		void PrintStats();
		void Increase(const VKString & counterName);
	
		TimeProfiler m_Timing;
		map<VKString, int> m_Counters;
	
		static AnalysisStats m_GlobalStats;
	
};

#endif




