#include "ParamParser.h"
#include <assert.h>

ParamParser::ParamParser(const VKString & inFilename, const char ** argv, int argc)
{
	std::ifstream file(inFilename.c_str());

	if (!file.is_open())
	{
		if (!argc>0)	//all params must be here
		{
			std::cout<<"[ERROR] Param File does not exist and no command line arguments provided: "<<inFilename.c_str()<<std::endl;
			assert(false);
		}
	}

	VKString currParam = "General";
		//read file, fill values
	VKStringList commands;
	VKString tempString(file);

	commands = tempString.replace("\n", " ").replace("\t", " ").split(" ");
	for (int i=0; i<argc; i++)	//read arguments, add or refill arguments
		commands.push_back(argv[i]);

	for (int i=0; i<commands.count(); i++)
	{		
//		std::cout<<"command["<<i<<"] = "<<commands[i].c_str()<<std::endl;
		if (commands[i].startsWith("-") || commands[i].endsWith(":"))
		{
			currParam = commands[i].replace("-", "").replace(":","");
//			std::cout<<"\tModule="<<currParam.c_str()<<std::endl;
		}
		else if (currParam!="Invalid")
		{
			VKString field = commands[i++];
//			std::cout<<"\tField="<<field.c_str()<<std::endl;
			m_ParamModules[currParam][field].clear();
			if (commands.count()<=i)
			{
				std::cout<<"Failed to read file: "<<inFilename.c_str()<<std::endl;
				assert(commands.count()>i);
			}
			if (commands[i].startsWith("[") && !commands[i].endsWith("]"))				//loop over multiple param values
			{
				assert(commands.count()>i);
				if (commands[i]=="[")
					i++;
				while(commands[i]!="]")
				{
					bool terminate = commands[i].endsWith("]");
					m_ParamModules[currParam][field].push_back(commands[i].replace("[", "").replace("]", ""));
//					std::cout<<"\tVal: "<<m_ParamModules[currParam][field][m_ParamModules[currParam][field].count()-1].c_str()<<std::flush;
					if (terminate)	break;
					i++;
				}

			}
			else		//set single param value
			{
				m_ParamModules[currParam][field].push_back(commands[i].replace("[", "").replace("]", ""));
//				std::cout<<"\tSingleton="<<m_ParamModules[currParam][field][m_ParamModules[currParam][field].count()-1].c_str()<<std::endl;				
			}
		}
	}
	//std::cout<<"done"<<std::endl;
}

void ParamParser::WriteParams(const VKString & outFile)
{
	std::ofstream textStream(outFile.c_str());
	assert(textStream.is_open());
	for (std::map<VKString, std::map<VKString, VKStringList> >::iterator iter = m_ParamModules.begin();
		 iter!=m_ParamModules.end(); iter++)
	{
		textStream<<iter->first.c_str()<<":\n";
		for (std::map<VKString, VKStringList>::iterator iter2 = iter->second.begin();
			 iter2!=iter->second.end(); iter2++)
		{
			textStream<<"\t"<<iter2->first.c_str()<<"\t["<<iter2->second.join(", ").c_str()<<"]\n";
		}
	}	
	textStream.close();
}

void ParamParser::PrintParams()
{
	for (std::map<VKString, std::map<VKString, VKStringList> >::iterator iter = m_ParamModules.begin();
		 iter!=m_ParamModules.end(); iter++)
	{
		std::cout<<"Module = "<<iter->first.c_str()<<std::endl;
		for (std::map<VKString, VKStringList>::iterator iter2 = iter->second.begin();
			 iter2!=iter->second.end(); iter2++)
		{
			std::cout<<"\t'"<<iter2->first.c_str()<<"' = [";
			std::cout<<iter2->second.join(", ").c_str()<<"]"<<std::endl;
		}
	}
}

const VKStringList & ParamParser::GetStrValues(const VKString & module, 
									   const VKString & valName, bool & valid)
{
	valid = m_ParamModules[module].find(valName)!=m_ParamModules[module].end();

	if (valid)
		return m_ParamModules[module][valName];
	else
		return m_EmptyList;
}

VKString ParamParser::GetStrValue(const VKString & module, 
								  const VKString & valName, bool & valid)
{	
	valid = m_ParamModules[module].find(valName)!=m_ParamModules[module].end();

	if (valid)
		return m_ParamModules[module][valName][0];
	else
		return VKString();
}

double ParamParser::GetDoubleValue(const VKString & module, 
								   const VKString & valName, bool & valid)
{
	valid = m_ParamModules[module].find(valName)!=m_ParamModules[module].end();
	if (valid)
		return m_ParamModules[module][valName][0].toDouble(&valid);
	else
		return -1.;
}

int ParamParser::GetIntValue(const VKString & module, 
							 const VKString & valName, bool & valid)
{
	valid = m_ParamModules[module].find(valName)!=m_ParamModules[module].end();
	if (valid)
		return m_ParamModules[module][valName][0].toInt(&valid);
	else
		return -1;
}

std::vector<double> ParamParser::GetDoubleValues(const VKString & module, 
												 const VKString & valName, bool & valid)
{
	valid = m_ParamModules[module].find(valName)!=m_ParamModules[module].end();
	std::vector<double> newVals;
	if (valid)
	{
		VKStringList list = m_ParamModules[module][valName];

		bool currValid;
		for (int i=0; i<list.count(); i++)
		{
			newVals.push_back(list[i].toDouble(&currValid));
			valid = currValid && valid;
		}
	}
	return newVals;
}

std::vector<int> ParamParser::GetIntValues(const VKString & module, 
										   const VKString & valName, bool & valid)
{
	valid = m_ParamModules[module].find(valName)!=m_ParamModules[module].end();
	std::vector<int> newVals;

	if (valid)
	{
		VKStringList list = m_ParamModules[module][valName];

		bool currValid;
		for (int i=0; i<list.count(); i++)
		{
			newVals.push_back(list[i].toInt(&currValid));
			valid = currValid && valid;
		}
	}
	return newVals;
}


