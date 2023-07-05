#include "find_root.h"
#include <cmath>
#include <stdexcept>
#include <utility>
#include <iostream>
#include <limits>

double BrentSolver::solve(std::function<double(double)> func)
{
    double reg_min_temp = reg_min;
    double reg_max_temp = reg_max;
    double fa = func(reg_min_temp);
    double fb = func(reg_max_temp);

    if (fa * fb > 0) {
        // throw std::runtime_error("Root must be bracketed in BrentSolver.");
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (std::fabs(fa) < std::fabs(fb)) {
        std::swap(reg_min_temp, reg_max_temp);
        std::swap(fa, fb);
    }

    double c = reg_min_temp, fc = fa;
    bool mflag = true;
    double s = 0, d = 0;

    for (int iter = 0; iter < maxIter; iter++) {
        double fs = func(s);

        if (fa != fc && fb != fc) {
            s = reg_min_temp * fb * fc / ((fa - fb) * (fa - fc)) + reg_max_temp * fa * fc / ((fb - fa) * (fb - fc)) + c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            s = reg_max_temp - fb * (reg_max_temp - reg_min_temp) / (fb - fa);
        }

        double tol1 = (3 * tol) / 2;
        double tol2 = 2 * tol1;

        if (!(3 * (reg_min_temp + tol1) < s && s < reg_max_temp - tol2) || (mflag && std::fabs(s - reg_max_temp) >= std::fabs(reg_max_temp - c) / 2) || (!mflag && std::fabs(s - reg_max_temp) >= std::fabs(c - d) / 2) || (mflag && std::fabs(reg_max_temp - c) < tol1) || (!mflag && std::fabs(c - d) < tol1)) {
            s = (reg_min_temp + reg_max_temp) / 2;
            mflag = true;
        } else {
            mflag = false;
        }

        fs = func(s);
        d = c;
        c = reg_max_temp;

        if (fa * fs < 0) {
            reg_max_temp = s;
            fb = fs;
        } else {
            reg_min_temp = s;
            fa = fs;
        }

        if (std::fabs(fa) < std::fabs(fb)) {
            std::swap(reg_min_temp, reg_max_temp);
            std::swap(fa, fb);
        }

        if (std::fabs(reg_max_temp - reg_min_temp) < tol) {
            return reg_max_temp;
        }
    }

    throw std::runtime_error("Maximum number of iterations exceeded in BrentSolver.");
}