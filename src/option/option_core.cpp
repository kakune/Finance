#include "option_core.h"
#include <algorithm>

double VirtualOption::payoff(double x)
{
    if ( isSum ) { return Option1->payoff(x) + Option2->payoff(x); }
    else { return Option1->payoff(x) - Option2->payoff(x); }
}

void VirtualOption::calcPayoffMean()
{
    Option1->calcPayoffMean();
    Option2->calcPayoffMean();
    if ( isSum )
    {
        setPayoffMean( Option1->getPayoffMean() + Option2->getPayoffMean() );
    }
    else
    {
        setPayoffMean( Option1->getPayoffMean() - Option2->getPayoffMean() );
    }
}


// std::shared_ptr<Option> operator-(Option& lhs, Option& rhs) {
//     return std::make_shared<VirtualOption>(&lhs, &rhs, false);
// }