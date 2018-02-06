#include "AnalysisStats.h"

AnalysisStats AnalysisStats::m_GlobalStats;

AnalysisStats::AnalysisStats()
{
}

void AnalysisStats::Increase(const VKString & counterName)
{
	map<VKString, int>::iterator iter = m_Counters.find(counterName);
	if (iter==m_Counters.end())
		m_Counters[counterName] = 1;
	else
		m_Counters[counterName]++;
}

void AnalysisStats::PrintStats()
{
	std::cout<<"AnalysisStats:"<<std::endl;
	std::cout<<AnalysisStats::m_GlobalStats.m_Timing.profile().c_str()<<std::endl;	
	for (map<VKString, int>::iterator iter=m_Counters.begin(); iter!=m_Counters.end(); iter++)
		std::cout<<iter->first.c_str()<<"\t"<<iter->second<<std::endl;
}

