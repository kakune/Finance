#ifndef ANALYTICAL_BS_H
#define ANALYTICAL_BS_H

#include <cmath>
#include <iostream>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <optional>
// #include "option_core.h"

class AnalyticalBlackScholes
{
private:
    double spot,strike,rate,vol,maturity;

public:
    void setSpot(double spot_){spot = spot_;}
    void setStrike(double strike_){strike = strike_;}
    void setRate(double rate_){rate = rate_;}
    void setVol(double vol_){vol = vol_;}
    void setMaturity(double maturity_){maturity = maturity_;}
    double getSpot(){return spot;}
    double getStrike(){return strike;}
    double getRate(){return rate;}
    double getVol(){return vol;}
    double getMaturity(){return maturity;}

    double dPlus()
    {
        return (std::log(spot/strike) + (rate + 0.5*vol*vol)*maturity) / (vol*sqrt(maturity));
    }

    double dMinus()
    {
        return dPlus() - vol*sqrt(maturity);
    }
    double payoffCall()
    {
        boost::math::normal_distribution<> dist(0, 1);
        double cdf_d1 = boost::math::cdf(dist, dPlus());
        double cdf_d2 = boost::math::cdf(dist, dMinus());
        return spot*cdf_d1 - strike*std::exp(-rate*maturity)*cdf_d2;
    }
    double payoffPut()
    {
        boost::math::normal_distribution<> dist(0, 1);
        double cdf_d1 = boost::math::cdf(dist, -dPlus());
        double cdf_d2 = boost::math::cdf(dist, -dMinus());
        return - spot*cdf_d1 + strike*std::exp(-rate*maturity)*cdf_d2;
    }
};

#endif