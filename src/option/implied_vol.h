#ifndef OPTION_IMPLIED_VOL_H
#define OPTION_IMPLIED_VOL_H

#include <analytical/black_scholes.h>
#include <numerical/find_root.h>

class ImpliedVolatility
{
private:
    int maxIter;
    double tol, rate, maturity;
    double callValue, spot, strike;
    double volMin, volMax;
    AnalyticalBlackScholes BSObj;
    BrentSolver SolverObj;

public:
    void setTol(double tol_ = 1e-5){tol = tol_; SolverObj.setTol(tol);}
    void setReg(double volMin_, double volMax_){volMin = volMin_; volMax = volMax_; SolverObj.setReg(volMin, volMax);}
    void setMaxIter(int maxIter_ = 100){maxIter = maxIter_; SolverObj.setMaxIter(maxIter);}
    void setRate(double rate_){rate = rate_; BSObj.setRate(rate);}
    void setMaturity(double maturity_){maturity = maturity_; BSObj.setMaturity(maturity);}
    void setCallValue(double callValue_){callValue = callValue_;}
    void setSpot(double spot_){spot = spot_; BSObj.setSpot(spot);}
    void setStrike(double strike_){strike = strike_; BSObj.setStrike(strike);}
    double getRate(){return rate;}
    double getMaturity(){return maturity;}
    double getCallValue(){return callValue;}
    double getStrike(){return strike;}

    double calcImpliedVol();


};

#endif