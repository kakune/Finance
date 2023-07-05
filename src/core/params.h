#ifndef PARAMS_MODULE_H
#define PARAMS_MODULE_H
#include <map>
#include <string>

class Parameters {
private:
    std::string variablesection, variablename, workingsection;
    double start, end, num_variable, step;
    int index_variable;
public:
    std::map< std::string, std::map<std::string, double> > data;
    std::map< std::string, std::map<std::string, std::string> > sdata;
    void readParameters(const std::string& filename);
    void setVariable(const std::string&, const std::string&, const std::string&);
    bool updateVariable();
    double getVariable(){return data[workingsection][variablename];}
};

#endif