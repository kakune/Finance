#include "vanilla.h"
#include <numerical/special_functions.h>
#include <numerical/integral.h>
#include <functional>
#define PI 3.1415926535897932384626433832795028841971

comp FourierPricingSV::phiZ ( comp v, comp u )
{
    comp thetap = theta - rho*eta*lambda*b*u;
    comp gamma = std::sqrt(thetap * thetap - 2.0 * eta * eta * v);
    comp egt = std::exp(-gamma * maturity);
    comp A = (theta / (eta*eta) ) * (2.0 * std::log(2.0*gamma / (thetap + gamma - egt * (thetap - gamma))) + (thetap - gamma) * maturity);
    comp B = 2.0 * v * ( 1.0 - egt ) / ((thetap + gamma) * (1.0 - egt) + 2.0 * gamma * egt);
    return std::exp(A + B);
}

comp FourierPricingSV::q ( comp u )
{
    comp var = 0.5 * lambda * lambda * b * b * u * (u - 1.0);
    return phiZ(var, u) - std::exp(var*maturity);
}

comp FourierPricingSV::integralR ( comp z, double a, double c )
{
    comp ia(0.0, a);
    comp eiaz = std::exp(ia*z);
    return (exponentialIntegralOne(z*(c - ia))/eiaz - exponentialIntegralOne(z*(c + ia))*eiaz) * 0.5 / ia;
}

double FourierPricingSV::calcCallValue()
{
    comp I(0.0, 1.0);
    double spotp = b*spot + (1.0 - b) * L;
    double strikep = b*strike + (1.0 - b) * L;
    double logsk = std::log(spotp / strikep);
    wMax = 5.0 / (eta * lambda * b * rhom * maturity) + 1.0;
    double thetaHat = theta - rho*eta*lambda*b*0.5;
    comp qZero = (rhom + I * rho)*thetaHat*(1.0 + theta*maturity) / (rhom * eta * eta) + (2.0*theta/(eta*eta)) * (std::log(2.0*rhom) + I * std::atan(rho / rhom));
    comp qInf = (rhom + I * rho)*(1.0 + theta*maturity) * lambda*b / eta;
    comp iThree = std::exp(qZero+logsk*0.5) * integralR(qInf - I*logsk, 0.5, wMax);

    std::size_t nOmega = 500;
    double deltaOmega = 0.1 / (lambda * b * std::sqrt(maturity));
    std::function<comp(comp)> integrandTwo = [&](comp w){return std::exp((0.5+I*w)*logsk)*(q(0.5 + I*w) - std::exp(-qInf*w + qZero)) / (w*w +0.25); };
    comp iTwo = oneDimIntegral(integrandTwo, comp(wMax, 0.0), comp(wMax, 0.0) + double(nOmega)*deltaOmega, nOmega);

    std::function<comp(comp)> integrandOne = [&](comp w){return std::exp((0.5+I*w)*logsk)*(q(0.5 + I*w)) / (w*w +0.25); };
    comp iOne = oneDimIntegral(integrandOne, comp(0.0, 0.0), comp(wMax, 0.0), nOmega);

    // #include <iostream>
    // std::cout << "real : " << std::real(iOne + iTwo + iThree) << std::endl;
    // std::cout << "imag : " << std::imag(iOne + iTwo + iThree) << std::endl;
    
    // std::function<comp(comp)> integrandPoor = [&](comp w){return std::exp((0.5 + I*w) * logsk)*(q(0.5 + I*w)) / (w*w + 0.25); };
    // comp iPoor =  oneDimIntegral(integrandPoor, comp(-100, 0.0), comp(100, 0.0), 100000);

    BSObj.setRate(rate);
    BSObj.setSpot(spotp);
    BSObj.setStrike(strikep);
    BSObj.setMaturity(maturity);
    BSObj.setVol(lambda*b);
    // return BSObj.payoffCall() / b - (strikep / (2.0 * PI * b)) * std::real(iPoor);
    return BSObj.payoffCall() / b - (strikep / (2.0 * PI * b)) * std::real(iOne + iTwo + iThree);
}


