#ifndef NUMERICAL_INTEGRAL_H
#define NUMERICAL_INTEGRAL_H
#include <vector>
#include <functional>
#include <gsl/gsl_integration.h>

double oneDimIntegralCquadGSL(std::function<double(double)>, double, double);
double oneDimIntegralQagiGSL(std::function<double(double)>);
double oneDimIntegralQagiuGSL(std::function<double(double)>, double);
double oneDimIntegralQagilGSL(std::function<double(double)>, double);

class OneDimIntegrationGSL
{
protected:
    double regMin, regMax;
    std::function<double(double)> integrand;
    double result, error;
    size_t nevals;
    gsl_function F;
    static double gslFunction(double x, void* params)
    {
        OneDimIntegrationGSL* integrator = static_cast<OneDimIntegrationGSL*>(params);
        return integrator->integrand(x);
    }
public:
    void setReg(double regMin_, double regMax_){regMin = regMin_; regMax = regMax_;}
    void setIntegrand(const std::function<double(double)>& integrand_){integrand = integrand_; F.function = &OneDimIntegrationGSL::gslFunction;
    F.params = this;}
    double cquadIntegral();
    double qagiIntegral();
    double qagiuIntegral();
    double qagilIntegral();
};


#include "integral.tpp"

#endif