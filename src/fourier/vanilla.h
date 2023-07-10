// This library is from Andersen & Piterbarg [2012] 

#ifndef FOURIER_VANILLA_H
#define FOURIER_VANILLA_H

#include <numerical/fft.h>
#include <analytical/black_scholes.h>
#include <complex>

typedef std::complex<double> comp;

class FourierPricingSV
{
private:
    double rate, maturity, lambda, b, theta, eta, L;
    double rho, rhoc;
    double spot, strike;
    double wMax;
    double tol;
    std::size_t nOmega;
public:
    FourierPricingSV( double spot_, double strike_, double rate_, double maturity_, double rho_, double lambda_, double L_, double b_, double theta_, double eta_, double tol_, std::size_t nOmega)
        : spot( spot_ ), strike( strike_ ), rate( rate_ ), maturity( maturity_ ), rho( rho_ ), rhoc( std::sqrt(1.0 - rho*rho) ), lambda( lambda_ ), L( L_ ), b( b_ ), theta( theta_ ), eta( eta_ ), tol( tol_ ), nOmega( nOmega ) {}

    double calcCallValue();
    comp q( const comp& u) const;
    comp phiZ( const comp& v, const comp& u ) const;
    comp integralR( const comp& z, const double& a, const double& c ) const;
};

class FourierPricingSVBuilder
{
private:
    double spot, strike;
    double rate = 0.0;
    double maturity, rho;
    double lambda, b, L, theta, eta;
    double tol = 1e-5;
    std::size_t nOmega = 50;
public:
    FourierPricingSVBuilder& setSpot( double spot_ )
    {
        spot = spot_;
        return *this;
    }
    FourierPricingSVBuilder& setStrike( double strike_ )
    {
        strike = strike_;
        return *this;
    }
    FourierPricingSVBuilder& setRate( double rate_ )
    {
        rate = rate_;
        return *this;
    }
    FourierPricingSVBuilder& setMaturity( double maturity_ )
    {
        maturity = maturity_;
        return *this;
    }
    FourierPricingSVBuilder& setRho( double rho_ )
    {
        rho = rho_;
        return *this;
    }
    FourierPricingSVBuilder& setLambda( double lambda_ )
    {
        lambda = lambda_;
        return *this;
    }
    FourierPricingSVBuilder& setL( double L_ )
    {
        L = L_;
        return *this;
    }
    FourierPricingSVBuilder& setB( double b_ )
    {
        b = b_;
        return *this;
    }
    FourierPricingSVBuilder& setTheta( double theta_ )
    {
        theta = theta_;
        return *this;
    }
    FourierPricingSVBuilder& setEta( double eta_ )
    {
        eta = eta_;
        return *this;
    }
    FourierPricingSVBuilder& setTol( double tol_ )
    {
        tol = tol_;
        return *this;
    }
    FourierPricingSVBuilder& setNOmega( std::size_t nOmega_ )
    {
        nOmega = nOmega_;
        return *this;
    }
    FourierPricingSV build()
    {
        return FourierPricingSV( spot, strike, rate, maturity, rho, lambda, L, b, theta, eta, tol, nOmega );
    }
};

#endif
