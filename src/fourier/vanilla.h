// This library is from Andersen & Piterbarg [2012] 

#ifndef FOURIER_VANILLA_H
#define FOURIER_VANILLA_H

#include <numerical/fft.h>
#include <analytical/black_scholes.h>
#include <complex>

typedef std::complex<double> comp;

class FourierPricing
{
protected:
    AnalyticalBlackScholes BSObj;
};

class FourierPricingSV : public FourierPricing
{
protected:
    double rate, maturity, lambda, b, theta, eta, L, initial_vol;
    double rho, rhom;
    double spot, strike;
    double wMax;
    double tol = 1e-4;
    std::size_t nOmega = 100;
public:
    void setSpot(double spot_){spot = spot_;}
    void setStrike(double strike_){strike = strike_;}
    void setRate(double rate_){rate = rate_;}
    void setRho(double rho_){rho = rho_; rhom = std::sqrt(1.0 - rho*rho);}
    void setTol(double tol_){tol = tol_;}
    void setNOmega(std::size_t nOmega_){nOmega = nOmega_;}
    void setParamsSpot(double lambda_, double b_, double L_)
    { lambda = lambda_; b = b_; L = L_; }
    void setParamsVol(double theta_, double eta_){theta = theta_; eta = eta_;}
    void setMaturity(double maturity_){maturity = maturity_;}
    double calcCallValue();
    comp q(comp);
    comp phiZ(comp, comp);
    comp integralR(comp, double, double);
};

#endif
