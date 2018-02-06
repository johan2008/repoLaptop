/***************************************************************************
 *   Copyright (C) 2008 by Vladimir Kim, Princeton University			   *
 *   Author: Vladimir Kim (kim_vladimir@yahoo.com, vk@princeton.edu)       *
 *   Research Supervisor: Tom Funkhouser (funk@cs.princeton.edu)           *
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

#include "TimeProfiler.h"
#ifdef CALCULATE_TIMING
	TimeProfiler * TimeProfiler::m_pGlobalSuperProfiler = new TimeProfiler();;
#endif

#include "VKStringList.h"
TimeProfiler::Process::Process()
{
	m_iCurrStartTime = -1;
	m_iTotalRunningTime = 0;
}

TimeProfiler::TimeProfiler()
{
}

void TimeProfiler::startedProcess(const VKString & processID)
{
	if (m_mIDToProcessMap[processID].m_iCurrStartTime!=-1)
		finishedProcess(processID);			//if process is currently running - finish it first
	
	m_mIDToProcessMap[processID].m_iCurrStartTime = clock();
	
	if (m_mIDToProcessMap[processID].m_sProcessName.isEmpty())
		m_mIDToProcessMap[processID].m_sProcessName = processID;
}

void TimeProfiler::finishedProcess(const VKString& processID, bool printTime)
{
	Process & process = m_mIDToProcessMap[processID];
	long clockStamp = clock();
	if (process.m_iCurrStartTime==-1)
		process.m_iCurrStartTime = clockStamp;
	long diff = clockStamp - process.m_iCurrStartTime;
	process.m_iTotalRunningTime += (diff-1);
	process.m_iCurrStartTime = -1;

	if (printTime)
		std::cout<<">>>>>>> Finished: "<<process.m_sProcessName.toStdString()<<" - "<<getTotalProcessLengthInSeconds(process.m_sProcessName).toStdString()<<std::endl;
}

long TimeProfiler::getTotalProcessLengthInCycles(const VKString& processID)
{
	return m_mIDToProcessMap[processID].m_iTotalRunningTime;
}

VKString TimeProfiler::getTotalProcessLengthInSeconds(const VKString& processID)
{
	long remainder = m_mIDToProcessMap[processID].m_iTotalRunningTime % CLOCKS_PER_SEC;
	remainder = (remainder * 1000) / CLOCKS_PER_SEC;
	VKString retVal= VKString::number(m_mIDToProcessMap[processID].m_iTotalRunningTime / CLOCKS_PER_SEC) + "." + VKString::number(remainder) + "s";
	return retVal;

}

int TimeProfiler::getTotalProcessLengthInSecondsInt(const VKString& processID)
{
	return m_mIDToProcessMap[processID].m_iTotalRunningTime / CLOCKS_PER_SEC;
}

VKString TimeProfiler::profile()
{
	if (m_mIDToProcessMap.size()==0)
		return "Empty Time Profile";
	VKString profile = "Time Profile:\n";
	for (std::map<VKString, Process>::iterator iter = m_mIDToProcessMap.begin(); 
		 iter!=m_mIDToProcessMap.end(); iter++)
	{
		VKString procName = iter->first;
		procName.replace(" ", "_");
		profile += VKString("    ") + procName;
		profile += VKString(" \t ") + getTotalProcessLengthInSeconds(iter->first) + "\n";
	}
	return profile;
}

void TimeProfiler::WriteProgress(const VKString & process, int currProgress, int maxProgress)
{
	finishedProcess(process);	
	int TICK_SLOTS = 40;
	int currTicks = TICK_SLOTS * currProgress / maxProgress;
	std::cout<<"\r[";
	for (int i=0; i<TICK_SLOTS; i++)
	{
		if (i < currTicks)
			std::cout<<"+"<<std::flush;
		else
			std::cout<<" "<<std::flush;
	}
	
	int currTime = getTotalProcessLengthInSecondsInt(process);
	int remaining = (maxProgress-currProgress);
	int remainingTime = ((currProgress==0) ? (-1) : (currTime * remaining / currProgress));
	std::cout<<"] N="<<remaining<<" T={"<<currTime<<"s | "<<remainingTime<<"s}    "<<std::flush;
	startedProcess(process);	
}

void TimeProfiler::WriteIntermediateTimes(VKStringList & processes, int N)
{
	std::cout<<"\r[";	
	for (int i=0; i<processes.count(); i++)
	{
		finishedProcess(processes[i]);
		std::cout<<(getTotalProcessLengthInSecondsInt(processes[i]) / N)<<"s\t"<<std::flush;
		startedProcess(processes[i]);	
	}
		
	std::cout<<"] / "<<N<<std::flush;
}







