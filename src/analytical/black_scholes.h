#ifndef ANALYTICAL_BS_H
#define ANALYTICAL_BS_H

#include <cmath>
#include <iostream>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <optional>

class AnalyticalBlackScholes
{
private:
    double spot, strike, vol, rate, maturity;

public:
    AnalyticalBlackScholes ( double spot_, double strike_, double vol_, double rate_, double maturity_ )
        : spot( spot_ ), strike( strike_ ), vol( vol_ ), rate( rate_ ), maturity( maturity_ ) {}
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
    double payoffCallForVol( double vol_ )
    {
        double temp = vol;
        vol = vol_;
        double result = payoffCall();
        vol = temp;
        return result;
    }
};

class AnalyticalBlackScholesBuilder
{
private:
    double spot, strike, vol, maturity;
    double rate = 0.0;
public:
    AnalyticalBlackScholesBuilder& setSpot( double spot_ )
    {
        spot = spot_;
        return *this;
    }
    AnalyticalBlackScholesBuilder& setStrike( double strike_ )
    {
        strike = strike_;
        return *this;
    }
    AnalyticalBlackScholesBuilder& setVol( double vol_ )
    {
        vol = vol_;
        return *this;
    }
    AnalyticalBlackScholesBuilder& setRate( double rate_ )
    {
        rate = rate_;
        return *this;
    }
    AnalyticalBlackScholesBuilder& setMaturity( double maturity_ )
    {
        maturity = maturity_;
        return *this;
    }
    AnalyticalBlackScholes build()
    {
        return AnalyticalBlackScholes( spot, strike, vol, rate, maturity );
    }
};


#endif