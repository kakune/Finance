#ifndef NUMERICAL_FIND_ROOT_H
#define NUMERICAL_FIND_ROOT_H
#include <functional>

class NonLinearSolver
{
protected:
    double reg_min,reg_max;
    double tol;
    int maxIter;
    std::function<double(double)> objectiveFunction;
public:
    NonLinearSolver(double tol_ = 1e-5, int maxIter_ = 100): tol(tol_), maxIter(maxIter_) {}
    void setReg(double reg_min_, double reg_max_){reg_min = reg_min_; reg_max = reg_max_;}
    void setTol(double tol_){tol = tol_;}
    void setMaxIter(int maxIter_){maxIter = maxIter_;}
    void setFunction(std::function<double(double)> func_){objectiveFunction = func_;}
    virtual double solve() = 0;
};

class BrentSolver : public NonLinearSolver
{
public:
    BrentSolver(double tol_ = 1e-5, int maxIter_ = 100) : NonLinearSolver(tol_, maxIter_) {}

    double solve() override;
};

#endif