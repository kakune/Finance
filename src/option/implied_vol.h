#ifndef IMPLIED_VOL_MODULE_H
#define IMPLIED_VOL_MODULE_H

#include <cmath>
#include <iostream>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <optional>
#include "option_core.h"


class AnalyticalBlackScholesCall
{
private:
    double spot,strike,rate,vol,maturity;

public:
    void set_spot(double spot_){spot = spot_;}
    void set_strike(double strike_){strike = strike_;}
    void set_rate(double rate_){rate = rate_;}
    void set_vol(double vol_){vol = vol_;}
    void set_maturity(double maturity_){maturity = maturity_;}

    double dPlus()
    {
        return (std::log(spot/strike) + (rate + 0.5*vol*vol)*maturity) / (vol*sqrt(maturity));
    }

    double dMinus()
    {
        return dPlus() - vol*sqrt(maturity);
    }
    double payoff()
    {
        boost::math::normal_distribution<> dist(0, 1);
        double cdf_d1 = boost::math::cdf(dist, dPlus());
        double cdf_d2 = boost::math::cdf(dist, dMinus());
        return spot*cdf_d1 - strike*std::exp(-rate*maturity)*cdf_d2;
    }
};

class ImpliedVolatility
{
private:
    int maxIter;
    double tol, rate;
    double vol_min, vol_max;
public:
    void set_tol(double tol_ = 1e-5){tol = tol_;}
    void set_reg(double vol_min_, double vol_max_){vol_min = vol_min_; vol_max = vol_max_;}
    void set_maxIter(int maxIter_ = 100){maxIter = maxIter_;}
};

#endif