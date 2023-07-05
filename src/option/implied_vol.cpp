#include "implied_vol.h"

double ImpliedVolatility::calcImpliedVol()
{
    auto objectiveFunction = [&](double sigma)
    {
        BSObj.setVol(sigma);
        return BSObj.payoffCall() - callValue;
    };
    SolverObj.setFunction(objectiveFunction);
    return SolverObj.solve();
}