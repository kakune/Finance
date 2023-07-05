#include "generator.h"

double AntitheticGaussGenerator::operator()() 
{
    if (index_next == num_time)
    {
       index_next = 0;
       is_odd_call = !is_odd_call;
       return (*this)();
    }
    if (is_odd_call)
    {
        last_values.at(index_next) = dist(engine);
        return last_values.at(index_next);
    }
    else
    {
        return -last_values.at(index_next);
    }
}