#ifndef NUMERICAL_SPECIAL_FUNCTIONS_H
#define NUMERICAL_SPECIAL_FUNCTIONS_H

#include <complex>

typedef std::complex<double> comp;

// See https://www.sci.utah.edu/~vpegorar/research/2011_JGT.pdf
comp exponentialIntegral(comp z, double tol = 1e-6, std::size_t maxIter = 10000);
comp exponentialIntegralOne(comp z, double tol = 1e-6, std::size_t maxIter = 10000);
double exponentialIntegral(double);
comp eiAsymptoticSeries(comp);
comp eiContinuedFraction(comp);
comp eiPowerSeries(comp);

class ExponentialIntegral
{
private:
    comp z;
    double rez, imz, r;
    double tol;
    std::size_t maxIter;
    bool isConverge = false;

public:
    ExponentialIntegral( comp z_, double tol_ = 1e-6, std::size_t maxIter_ = 10000 )
        : z( z_ ), tol( tol_ ), maxIter( maxIter_ ) {}
    comp eiAsymptoticSeries();
    comp eiContinuedFraction();
    comp eiPowerSeries();
    comp calcExponentialIntegral();
    bool isConverged(comp, comp);
};

#endif