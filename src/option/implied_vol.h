#ifndef OPTION_IMPLIED_VOL_H
#define OPTION_IMPLIED_VOL_H

#include <analytical/black_scholes.h>
#include <numerical/find_root.h>

class ImpliedVolatility
{
private:
    std::size_t maxIter;
    double tol, rate, maturity;
    double spot, strike;
    double volMin, volMax;
    AnalyticalBlackScholes BSObj;
    BrentSolverBuilder SolverBuilder;

public:
    ImpliedVolatility( double spot_, double strike_, double rate_, double maturity_, double volMin_, double volMax_, double tol_, std::size_t maxIter_ );

    double calcImpliedVol( double callValue );
};

class ImpliedVolatilityBuilder
{
private:
    double rate = 0.0; 
    double maturity;
    double spot, strike;
    double volMin, volMax;
    double tol = 1e-5;
    std::size_t maxIter = 100;
public:
    ImpliedVolatilityBuilder& setSpot( double spot_ )
    {
        spot = spot_;
        return *this;
    }
    ImpliedVolatilityBuilder& setStrike( double strike_ )
    {
        strike = strike_;
        return *this;
    }
    ImpliedVolatilityBuilder& setRate( double rate_ )
    {
        rate = rate_;
        return *this;
    }
    ImpliedVolatilityBuilder& setMaturity( double maturity_ )
    {
        maturity = maturity_;
        return *this;
    }
    ImpliedVolatilityBuilder& setVolReg( double volMin_, double volMax_ )
    {
        volMin = volMin_;
        volMax = volMax_;
        return *this;
    }
    ImpliedVolatilityBuilder& setTol( double tol_ )
    {
        tol = tol_;
        return *this;
    }
    ImpliedVolatilityBuilder& setMaxIter( std::size_t maxIter_ )
    {
        maxIter = maxIter_;
        return *this;
    }
    ImpliedVolatility build()
    {
        return ImpliedVolatility( spot, strike, rate, maturity, volMin, volMax, tol, maxIter );
    }
};

#endif