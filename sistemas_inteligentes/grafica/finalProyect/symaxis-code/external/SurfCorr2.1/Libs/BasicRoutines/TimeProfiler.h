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

#ifndef __TIME_PROFILER_H
#define __TIME_PROFILER_H

#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include <assert.h>
#include "VKString.h"


//#define CALCULATE_TIMING

class TimeProfiler
{
	public:
		TimeProfiler();
		void startedProcess(const VKString & processID);
		void finishedProcess(const VKString & processID, bool printTime = false);
		long getTotalProcessLengthInCycles(const VKString & processID);
		VKString getTotalProcessLengthInSeconds(const VKString & processID);
		int getTotalProcessLengthInSecondsInt(const VKString & processID);
		VKString profile();
	
		void WriteProgress(const VKString & process, int currProgress, int maxProgress);
		void WriteIntermediateTimes(VKStringList & processes, int N);

#ifdef CALCULATE_TIMING
		static TimeProfiler * m_pGlobalSuperProfiler;
#endif

	protected:
		struct Process
		{
			Process();
			long m_iCurrStartTime;
//			std::vector<long> m_iRunningTimes;
			long m_iTotalRunningTime;
			VKString m_sProcessName;
		};
		std::map<VKString, Process> m_mIDToProcessMap;
};

#endif

