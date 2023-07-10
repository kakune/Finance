#include "special_functions.h"
#define PI 3.14159265358979323846264338327950288
#define GAMMA 0.57721566490153286060651209008240243


double exponentialIntegralOne( double x )
{
    return -std::expint( -x );
}

comp exponentialIntegral( comp z, double tol, std::size_t maxIter )
{
    ExponentialIntegral EiObj( z, tol, maxIter );
    return EiObj.calcExponentialIntegral();
}

comp exponentialIntegralOne(comp z, double tol, std::size_t maxIter)
{
    comp dif(0.0, ((std::arg(z) > 0) - ((std::arg(z) < 0)))*PI);
    return - exponentialIntegral(-z, tol, maxIter) - dif;
}

comp ExponentialIntegral::calcExponentialIntegral()
{
    rez = std::real(z);
    imz = std::imag(z);
    r = std::abs(z);
    if ( imz == 0.0 )
    {
        return comp(std::expint(rez), 0.0);
    }
    if ( r > 2.0 - 1.035 * std::log(tol) )
    {
        comp ei = eiAsymptoticSeries();
        if ( isConverge ) { return ei; }
    }
    if ( r > 1.0 && ( rez < 0 || std::abs(imz) > 1 ))
    {
        return eiContinuedFraction();
    }
    return eiPowerSeries();
}

bool ExponentialIntegral::isConverged(comp x, comp y)
{
    return ((std::abs(std::real(x - y))) <= tol * std::abs(std::real(x))) && ((std::abs(std::imag(x - y))) <= tol * std::abs(std::imag(x)));
}

comp ExponentialIntegral::eiAsymptoticSeries()
{
    comp ei(0.0, ((imz > 0)-(imz < 0)) * PI);
    comp tmp = std::exp(z) / z;
    comp old;
    for( int k = 1; k <= floor(r)+1; k++)
    {
        old = ei;
        ei += tmp;
        if ( isConverged(ei, old) )
        {
            isConverge = true;
            return ei;
        }
        tmp *= double(k) / z;
    }
    return ei;
}

comp ExponentialIntegral::eiContinuedFraction()
{
    comp ei(0.0, ((imz > 0)-(imz < 0)) * PI);
    comp c,d,old;
    if ( std::abs(ei) != 0.0 )
    {
        c = 1.0 / (1.0 - z - std::exp(z) / ei);
        d = 1.0 / (1.0 - z);
    }
    else
    {
        c = 0.0;
        d = 1.0 / (1.0 - z);
        ei = d * (- std::exp(z));
    }
    for ( int k = 1; k < maxIter; k++ )
    {
        c = 1.0 / (2.0 * k + 1.0 - z - double(k*k) * c);
        d = 1.0 / (2.0 * k + 1.0 - z - double(k*k) * d);
        old = ei;
        ei *= d / c;
        if (isConverged(ei, old)) {break;}
    }
    return ei;
}

comp ExponentialIntegral::eiPowerSeries()
{
    comp ei( GAMMA + std::log(r), ((imz > 0)-(imz < 0)) * std::abs(std::arg(z)) );
    comp tmp(1.0, 0.0);
    comp old;
    for ( int k = 1; k < maxIter; k++ )
    {
        tmp *= z / double(k);
        old = ei;
        ei += tmp / double(k);
        if (isConverged(ei, old)) {break;}
    }
    return ei;
}


// int main()
// {
//     comp var(-0.5, 0.04);
//     std::cout << exponentialIntegralOne(var) << std::endl;
//     return 0;
// }