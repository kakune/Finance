#ifndef MONTE_CARLO_PROCESS_H
#define MONTE_CARLO_PROCESS_H

#include "generator.h"
#include <cstddef>
#include <vector>
#include <cmath>
#include <random>

class SpotProcess
{
private:
    AntitheticGaussGenerator Generator;
    std::size_t numProcess, numTime;
    double dt, maturity;
    double initialSpot;
    std::vector< std::vector<double> > spots;
public:
    SpotProcess( double initialSpot_, double maturity_, std::size_t numProcess_, std::size_t numTime_, std::size_t numGeneratorPeriod_ );
    // void setSpots()
    // {
    //     spots.resize(numProcess);
    //     for(int i = 0; i < numProcess; i++)
    //     {
    //         spots.at(i).resize(numTime+1, initialSpot);
    //     }
    // }
    // void setTime( const std::size_t& numTime_, const double& maturity_)
    // {
    //     numTime = numTime_;
    //     maturity = maturity_;
    //     dt = maturity / double(numTime);
    //     Generator.setNumTime(numTime);
    // }
    void setSpotsElement( const std::size_t& i, const std::size_t& j, const double& x ) { spots.at(i).at(j) = x; }
    // void setInitialSpot( const double& initialSpot_){initialSpot = initialSpot_;}
    // void setNumProcess( const std::size_t& numProcess_){numProcess = numProcess_;}
    // void setGeneratorNumTime( const std::size_t& numTime_ ){ Generator.setNumTime( numTime_ ); }

    std::vector< std::vector<double> >& getSpots() { return spots; }
    std::size_t getNumProcess() const { return numProcess; }
    std::size_t getNumTime() const { return numTime; }
    double getSpotsElement( const std::size_t& i, const std::size_t& j ) const { return spots.at(i).at(j); }
    double getDt() const { return dt; }
    double getRandom() { return Generator(); }

    virtual std::vector< std::vector<double> >& generateSpots() = 0;
};

class BlackScholes : public SpotProcess
{
private:
    double vol, rate;
public:
    BlackScholes( double initialSpot_, double maturity_, std::size_t numProcess_, std::size_t numTime_, std::size_t numGeneratorPeriod_, double vol_, double rate_ );
    std::vector< std::vector<double> >& generateSpots();
};

class StochasticVolatilityProcess : public SpotProcess
{
private:
    std::vector< std::vector<double> > vols;
    double rho, rhoc, initialVol;
public:
    StochasticVolatilityProcess( double initialSpot_, double maturity_, std::size_t numProcess_, std::size_t numTime_, std::size_t numGeneratorPeriod_, double initialVol_, double rho_ );
    // void setVols()
    // {
    //     vols.resize( getNumProcess() );
    //     for(int i = 0; i < getNumProcess(); i++)
    //     {
    //         vols.at(i).resize( getNumTime() + 1, initialVol);
    //     }
    // }
    // void setRho( const double& rho_)
    // {
    //     rho = rho_;
    //     rhoc = std::sqrt(1.0 - rho*rho);
    // }
    // void setInitialVol( const double& initialVol_ ) { initialVol = initialVol_; }
    void setVolsElement( const std::size_t& i, const std::size_t& j, const double& x ) { vols.at(i).at(j) = x; }

    double getRho() const { return rho; }
    double getRhoc() const { return rhoc; }
    double getVolsElement( const std::size_t& i, const std::size_t& j ) const { return vols.at(i).at(j); }
};

class SABR : public StochasticVolatilityProcess
{
private:
    double rate, beta, volvol;
public:
    SABR( double initialSpot_, double maturity_, std::size_t numProcess_, std::size_t numTime_, std::size_t numGeneratorPeriod_, double initialVol_, double rho_, double rate_, double beta_, double volvol_ );
    std::vector< std::vector<double> >& generateSpots();
};

class SV : public StochasticVolatilityProcess
{
private:
    double rate, lambda, b, L, theta, eta;
public:
    SV( double initialSpot_, double maturity_, std::size_t numProcess_, std::size_t numTime_, std::size_t numGeneratorPeriod_, double rho_, double rate_, double lambda_, double b_, double L_, double theta_, double eta_ );
    // void setRate( const double& rate_ ) { rate = rate_; }
    // void setLambda( const double& lambda_ ) { lambda = lambda_; }
    // void setB( const double& b_ ) { b = b_; }
    // void setL( const double& L_ ) { L = L_; }
    // void setTheta( const double& theta_ ) { theta = theta_; }
    // void setEta( const double& eta_ ) { eta = eta_; }
    std::vector< std::vector<double> >& generateSpots();
};

class SpotProcessBuilder
{
protected:
    double initialSpot, maturity;
    std::size_t numProcess, numTime, numGeneratorPeriod;
    bool isSetPeriod = false;
public:
    SpotProcessBuilder& setInitialSpot( double initialSpot_ )
    {
        initialSpot = initialSpot_;
        return *this;
    }
    SpotProcessBuilder& setMaturity( double maturity_ )
    {
        maturity = maturity_;
        return *this;
    }
    SpotProcessBuilder& setNumProcess( std::size_t numProcess_ )
    {
        numProcess = numProcess_;
        return *this;
    }
    SpotProcessBuilder& setNumTime( std::size_t numTime_ )
    {
        numTime = numTime_;
        return *this;
    }
    SpotProcessBuilder& setNumGeneratorPeriod( std::size_t numGeneratorPeriod_ )
    {
        isSetPeriod = true;
        numGeneratorPeriod = numGeneratorPeriod_;
        return *this;
    }
};

class BlackScholesBuilder : public SpotProcessBuilder
{
protected:
    double vol, rate;
public:
    BlackScholesBuilder& setVol( double vol_ )
    {
        vol = vol_;
        return *this;
    }
    BlackScholesBuilder& setRate( double rate_ )
    {
        rate = rate_;
        return *this;
    }
    BlackScholes build()
    {
        if( !isSetPeriod ){ numGeneratorPeriod = numTime; }
        return BlackScholes( initialSpot, maturity, numProcess, numTime, numGeneratorPeriod, vol, rate );
    }
};

class StochasticVolatilityProcessBuilder : public SpotProcessBuilder
{
protected:
    double rho, initialVol;
public:
    StochasticVolatilityProcessBuilder& setInitialVol( double initialVol_ )
    {
        initialVol = initialVol_;
        return *this;
    }
    StochasticVolatilityProcessBuilder& setRho( double rho_ )
    {
        rho = rho_;
        return *this;
    }
};

class SABRBuilder : public StochasticVolatilityProcessBuilder
{
protected:
    double rate, beta, volvol;
public:
    SABRBuilder& setRate( double rate_ )
    {
        rate = rate_;
        return *this;
    }
    SABRBuilder& setBeta( double beta_ )
    {
        beta = beta_;
        return *this;
    }
    SABRBuilder& setVolvol( double volvol_ )
    {
        volvol = volvol_;
        return *this;
    }
    SABR build()
    {
        if( !isSetPeriod ){ numGeneratorPeriod = numTime * 2; }
        return SABR( initialSpot, maturity, numProcess, numTime, numGeneratorPeriod, initialVol, rho, rate, beta, volvol );
    }
};

class SVBuilder : public StochasticVolatilityProcessBuilder
{
protected:
    double rate, lambda, b, L, theta, eta;
public:
    SVBuilder& setRate( double rate_ )
    {
        rate = rate_;
        return *this;
    }
    SVBuilder& setLambda( double lambda_ )
    {
        lambda = lambda_;
        return *this;
    }
    SVBuilder& setB( double b_ )
    {
        b = b_;
        return *this;
    }
    SVBuilder& setL( double L_ )
    {
        L = L_;
        return *this;
    }
    SVBuilder& setTheta( double theta_ )
    {
        theta = theta_;
        return *this;
    }
    SVBuilder& setEta( double eta_ )
    {
        eta = eta_;
        return *this;
    }
    SV build()
    {
        if( !isSetPeriod ){ numGeneratorPeriod = numTime * 2; }
        return SV( initialSpot, maturity, numProcess, numTime, numGeneratorPeriod, rho, rate, lambda, b, L, theta, eta );
    }
};

#endif