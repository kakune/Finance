#include "integral.h"
#include <complex>
#include <iostream>

template <typename T>
T OneDimIntegration<T>::calcTrapzIntegral()
{
    T result;
    result = 0.5 * integrand(regMin);
    T width = (regMax - regMin) / double(nStep);
    T variable = regMin;
    for(int i_ = 0; i_ < nStep - 1; i_++)
    {
        variable += width;
        result += integrand(variable);
    }
    result += 0.5 * integrand(regMax);
    return result / double(nStep);
}

// int main()
// {
//     std::function<std::complex<double>(std::complex<double>)> f1 = [](std::complex<double> x)
//     {
//         std::complex<double> I(0.0, 1.0);
//         return std::exp(I*x);
//     };
//     OneDimIntegration<std::complex<double>> IntObj;
//     std::complex<double> regMin(0.0, 0.0), regMax(1.0, 0.0);
//     IntObj.setReg(regMin, regMax);
//     IntObj.setNStep(1000);
//     IntObj.setIntegrand(f1);
//     std::cout << IntObj.calcTrapzIntegral() << std::endl;
//     return 0;
// }