#include "find_root.h"
#include <cmath>
#include <stdexcept>
#include <utility>
#include <iostream>
#include <limits>

double BrentSolver::solve()
{
    double regMinTemp = getRegMin();
    double regMaxTemp = getRegMax();
    double fa = calcObjectiveFunction( regMinTemp );
    double fb = calcObjectiveFunction( regMaxTemp );

    if (fa * fb > 0)
    {
        // throw std::runtime_error("Root must be bracketed in BrentSolver.");
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (std::fabs(fa) < std::fabs(fb))
    {
        std::swap(regMinTemp, regMaxTemp);
        std::swap(fa, fb);
    }

    double c = regMinTemp, fc = fa;
    bool mflag = true;
    double s = 0, d = 0;

    for (int iter = 0; iter < getMaxIter(); iter++) {
        double fs = calcObjectiveFunction(s);

        if (fa != fc && fb != fc)
        {
            s = regMinTemp * fb * fc / ((fa - fb) * (fa - fc)) + regMaxTemp * fa * fc / ((fb - fa) * (fb - fc)) + c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            s = regMaxTemp - fb * (regMaxTemp - regMinTemp) / (fb - fa);
        }

        double tol1 = ( 3 * getTol() ) / 2;
        double tol2 = 2 * tol1;

        if ( !(3 * (regMinTemp + tol1) < s && s < regMaxTemp - tol2)
            || (mflag && std::fabs(s - regMaxTemp) >= std::fabs(regMaxTemp - c) / 2)
            || (!mflag && std::fabs(s - regMaxTemp) >= std::fabs(c - d) / 2)
            || (mflag && std::fabs(regMaxTemp - c) < tol1)
            || (!mflag && std::fabs(c - d) < tol1) ) 
        {
            s = (regMinTemp + regMaxTemp) / 2;
            mflag = true;
        }
        else
        {
            mflag = false;
        }

        fs = calcObjectiveFunction(s);
        d = c;
        c = regMaxTemp;

        if (fa * fs < 0)
        {
            regMaxTemp = s;
            fb = fs;
        }
        else
        {
            regMinTemp = s;
            fa = fs;
        }

        if ( std::fabs(fa) < std::fabs(fb) )
        {
            std::swap(regMinTemp, regMaxTemp);
            std::swap(fa, fb);
        }

        if ( std::fabs(regMaxTemp - regMinTemp) < getTol() )
        {
            return regMaxTemp;
        }
    }

    throw std::runtime_error("Maximum number of iterations exceeded in BrentSolver.");
}