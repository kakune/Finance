#ifndef EUROPEAN_MODULE_H
#define EUROPEAN_MODULE_H

#include <vector>
#include <cstddef>
#include "option_core.h"

class EuropeanOption : public Option
{
public:
    void set_spot_paths(std::vector< std::vector< double > > spot_paths_) override
    {
        spot_paths = spot_paths_;
        num_path = spot_paths.size();
        present_spots.resize(num_path);
        for(int i = 0; i < num_path; i++)
        {
            present_spots.at(i) = spot_paths.at(i).back();
        }
    }
    void calc_payoff_mean() override
    {
        double sum = 0;
        for(int i = 0; i < num_path; i++)
        {
            sum += payoff(present_spots.at(i));
        }
        payoff_mean = sum / double(num_path);
    }
};

class EuropeanCallOption : public EuropeanOption, public Call
{
public:
    double payoff(double x) override { return Call::payoff(x); }
};

class EuropeanPutOption : public EuropeanOption, public Put
{
public:
    double payoff(double x) override { return Put::payoff(x); }
};



#endif