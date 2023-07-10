#include "vanilla.h"
#include <numerical/special_functions.h>
#include <numerical/integral.h>
#include <functional>

comp FourierPricingSV::phiZ( const comp& v, const comp& u ) const
{
    comp thetap = theta - rho*eta*lambda*b*u;
    comp gamma = std::sqrt(thetap * thetap - 2.0 * eta * eta * v);
    comp egt = std::exp(-gamma * maturity);
    comp A = (theta / (eta*eta) ) * (2.0 * std::log(2.0*gamma / (thetap + gamma - egt * (thetap - gamma))) + (thetap - gamma) * maturity);
    comp B = 2.0 * v * ( 1.0 - egt ) / ((thetap + gamma) * (1.0 - egt) + 2.0 * gamma * egt);
    return std::exp(A + B);
}

comp FourierPricingSV::q( const comp& u ) const
{
    comp var = 0.5 * lambda * lambda * b * b * u * (u - 1.0);
    return phiZ(var, u) - std::exp(var*maturity);
}

comp FourierPricingSV::integralR( const comp& z, const double& a, const double& c ) const
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
    double thetaHat = theta - rho*eta*lambda*b*0.5;
    double mu = thetaHat*thetaHat*0.5 / (eta*eta*lambda*lambda*b*b*rhoc*rhoc*rhoc) + 0.125 / rhoc;

    comp qZero = (rhoc + I * rho)*thetaHat*(1.0 + theta*maturity) / (rhoc * eta * eta) + (2.0*theta/(eta*eta)) * (std::log(2.0*rhoc) + I * std::atan(rho / rhoc));
    comp qInf = (rhoc + I * rho)*(1.0 + theta*maturity) * lambda*b / eta;
    comp qm = (theta/(eta*eta))*(maturity*mu*eta*lambda*b + 2.0*thetaHat*(rhoc * I*rho) / (rhoc*rhoc*eta*lambda*b)) + mu*lambda*b/eta;

    wMax = 5.0 / (eta * lambda * b * rhoc * maturity);
    while(std::abs(qm) > std::abs(-qInf*wMax + qZero)*tol*wMax)
    {
        wMax *= 2.0;
    }

    double iThree = std::real(std::exp(qZero+logsk*0.5) * integralR(qInf - I*logsk, 0.5, wMax));
    double deltaOmega = 0.1 / (lambda * b * std::sqrt(maturity));

    std::function<double(double)> integrandTwo = [&](double w){return std::real(std::exp((0.5+I*w)*logsk)*(q(0.5 + I*w) - std::exp(-qInf*w + qZero))) / (w*w +0.25); };
    double iTwo = oneDimIntegralCquadGSL(integrandTwo, wMax, wMax + double(nOmega)*deltaOmega);

    std::function<double(double)> integrandOne = [&](double w){return std::real(std::exp((0.5 + I*w) * logsk)*(q(0.5 + I*w)) / (w*w + 0.25)); };
    double iOne = oneDimIntegralCquadGSL(integrandOne, 0.0, wMax);


    AnalyticalBlackScholesBuilder BSBuilder;
    BSBuilder.setRate( rate );
    BSBuilder.setSpot( spotp );
    BSBuilder.setStrike( strikep );
    BSBuilder.setMaturity( maturity );
    BSBuilder.setVol( lambda*b );
    
    return BSBuilder.build().payoffCall() / b - (strikep / (M_PI * b)) * (iOne + iTwo + iThree);
}


