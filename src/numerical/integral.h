#ifndef INTEGRAL_MODULE_H
#define INTEGRAL_MODULE_H
#include <vector>
#include <functional>

template <typename T>
class OneDimIntegration
{
protected:
    T regMin, regMax;
    std::size_t nStep;
    std::function<T(T)> integrand;
public:
    void setReg(T regMin_, T regMax_){regMin = regMin_; regMax = regMax_;}
    void setNStep(std::size_t nStep_){nStep = nStep_;}
    void setIntegrand(std::function<T(T)> integrand_){integrand = integrand_;}
    T calcTrapzIntegral();
};

#endif