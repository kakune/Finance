#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "params.h"

void Parameters::readParameters(const std::string& filename) {
    std::ifstream param_file(filename);
    if ( !param_file )
    {
        std::cerr << "Could not open parameter file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::string section = "DEFAULT";
    while ( std::getline(param_file, line) )
    {
        std::istringstream iss(line);
        std::string name, delim;
        std::string svalue;
        if(line[0] == '[' && line[line.size() - 1] == ']')
        {
            section = line.substr(1, line.size() - 2);
        }
        if(iss >> name >> delim >> svalue)
        {
            try
            {
                data[section][name] = std::stod(svalue);
            }
            catch(const std::exception& e)
            {
                sdata[section][name] = svalue;
            }
        }
    }
    param_file.close();
}
void Parameters::setVariable(const std::string& variablesection_, const std::string& variableindex, const std::string& workingsection_)
{
    variablesection = variablesection_;
    variablename = sdata[variablesection][variableindex];
    workingsection = workingsection_;
    index_variable = 0;
    start = data[variablesection]["start"];
    end = data[variablesection]["end"];
    num_variable = data[variablesection]["num_variable"];

    data[workingsection][variablename] = start;
    step = (end - start) / double(num_variable - 1);
}
bool Parameters::updateVariable()
{
    if ( index_variable >= num_variable ) {return false;}
    data[workingsection][variablename] = start + step * index_variable;
    index_variable++;
    return true;
}