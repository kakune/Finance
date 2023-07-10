#ifndef NUMERICAL_INTEGRAL_H
#define NUMERICAL_INTEGRAL_H
#include <vector>
#include <functional>
#include <gsl/gsl_integration.h>

double oneDimIntegralCquadGSL(const std::function<double(double)>&, double, double);
double oneDimIntegralQagiGSL(const std::function<double(double)>&);
double oneDimIntegralQagiuGSL(const std::function<double(double)>&, double);
double oneDimIntegralQagilGSL(const std::function<double(double)>&, double);

class OneDimIntegrationGSL
{
private:
    double regMin, regMax;
    std::function< double( double ) > integrand;
    double result, error;
    size_t nevals;
    gsl_function F;
    static double gslFunction(double x, void* params)
    {
        OneDimIntegrationGSL* integrator = static_cast<OneDimIntegrationGSL*>(params);
        return integrator->integrand(x);
    }
public:
    OneDimIntegrationGSL( const std::function<double(double)>& integrand_, double regMin_, double regMax_ )
        : integrand( integrand_ ), regMin( regMin_ ), regMax( regMax_ )
    {
        F.function = &OneDimIntegrationGSL::gslFunction;
        F.params = this;
    }

    double calcIntegrand( double x ) { return integrand(x); }
    double cquadIntegral();
    double qagiIntegral();
    double qagiuIntegral();
    double qagilIntegral();
};


#include "integral.tpp"

#endif