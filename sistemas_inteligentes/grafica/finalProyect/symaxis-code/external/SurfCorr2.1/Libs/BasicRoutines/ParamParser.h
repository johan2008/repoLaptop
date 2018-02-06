#include <vector>
#include <map>
#include <assert.h>
#include <iostream>
#include <fstream>

#include "VKString.h"
#include "VKStringList.h"



#ifndef __PARAM_PARESER_H
#define __PARAM_PARESER_H

class ParamParser
{
	public:
		ParamParser(const VKString & inFilename, const char ** argv = NULL, int argc = 0);
		ParamParser(){}
		void PrintParams();
		void WriteParams(const VKString & outFile);

		std::map<VKString, std::map<VKString, VKStringList> > m_ParamModules;

		const VKStringList & GetStrValues(const VKString & module, 
								  const VKString &  valName, bool & valid);
	
		VKString GetStrValue(const VKString & module, 
							 const VKString &valName, bool & valid);
	
		double GetDoubleValue(const VKString & module, 
							  const VKString & valName, bool & valid);
	
		int GetIntValue(const VKString & module, 
						const VKString & valName, bool & valid);
	
		std::vector<double> GetDoubleValues(const VKString & module, 
											const VKString & valName, bool & valid);
	
		std::vector<int> GetIntValues(const VKString & module, 
									  const VKString & valName, bool & valid);
		VKStringList m_EmptyList;
};

#endif

