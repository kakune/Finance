#ifndef FIND_ROOT_MODULE_H
#define FIND_ROOT_MODULE_H
#include <functional>

class NonLinearSolver
{
protected:
    double reg_min,reg_max;
    double tol;
    int maxIter;
public:
    NonLinearSolver(double reg_min_, double reg_max_, double tol_ = 1e-5, int maxIter_ = 100): reg_min(reg_min_), reg_max(reg_max_), tol(tol_), maxIter(maxIter_) {}
    void set_reg(double reg_min_, double reg_max_){reg_min = reg_min_; reg_max = reg_max_;}
    void set_tol(double tol_){tol = tol_;}
    void set_maxIter(int maxIter_){maxIter = maxIter_;}
    virtual double solve(std::function<double(double)> func) = 0;
};

class BrentSolver : public NonLinearSolver
{
public:
    BrentSolver(double reg_min_, double reg_max_, double tol_ = 1e-5, int maxIter_ = 100) : NonLinearSolver(reg_min_, reg_max_, tol_, maxIter_) {}

    double solve(std::function<double(double)> func) override;
};

#endif