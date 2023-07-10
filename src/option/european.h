#ifndef OPTION_EUROPEAN_H
#define OPTION_EUROPEAN_H

#include <vector>
#include <cstddef>
#include "option_core.h"

class EuropeanOption : public Option
{
public:
    void calcPayoffMean() override
    {
        double sum = 0;
        for ( auto spot : getPresentSpots() )
        {
            sum += payoff( spot );
        }
        setPayoffMean( sum / double( getNumPath() ) );
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