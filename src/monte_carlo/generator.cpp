#include "generator.h"

double AntitheticGaussGenerator::operator()() 
{
    if ( indexNext == getNumTime() )
    {
       indexNext = 0;
       isOdd = !isOdd;
       return (*this)();
    }
    if ( isOdd )
    {
        lastValues.at( indexNext ) = GaussGenerator::operator()();
        return lastValues.at( indexNext );
    }
    else
    {
        return -lastValues.at( indexNext );
    }
}