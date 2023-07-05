#include "option_core.h"
#include <algorithm>

double VirtualOption::payoff(double x)
{
    if(isSum) { return Option1->payoff(x) + Option2->payoff(x); }
    else { return Option1->payoff(x) - Option2->payoff(x); }
}

void VirtualOption::calc_payoff_mean()
{
    Option1->calc_payoff_mean();
    Option2->calc_payoff_mean();
    if(isSum)
    {
        payoff_mean = Option1->get_payoff_mean() + Option2->get_payoff_mean();
    }
    else
    {
        payoff_mean = Option1->get_payoff_mean() - Option2->get_payoff_mean();
    }
}


// std::shared_ptr<Option> operator-(Option& lhs, Option& rhs) {
//     return std::make_shared<VirtualOption>(&lhs, &rhs, false);
// }