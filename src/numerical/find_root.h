#ifndef NUMERICAL_FIND_ROOT_H
#define NUMERICAL_FIND_ROOT_H
#include <functional>

class NonLinearSolver
{
private:
    std::function< double( double ) > objectiveFunction;
    double regMin,regMax;
    double tol;
    std::size_t maxIter;

public:
    NonLinearSolver( const std::function< double( double ) >& objectiveFunction_,
                    double regMin_,
                    double regMax_,
                    double tol_,
                    std::size_t maxIter_ )
                    : objectiveFunction( objectiveFunction_ ),
                    regMin( regMin_ ),
                    regMax( regMax_ ),
                    tol( tol_ ),
                    maxIter( maxIter_ ) {}

    double getRegMin() const { return regMin; }
    double getRegMax() const { return regMax; }
    double getTol() const { return tol; }
    std::size_t getMaxIter() const { return maxIter; }
    double calcObjectiveFunction( double x ) const { return objectiveFunction( x ); }
    virtual double solve() = 0;
};

class BrentSolver : public NonLinearSolver
{
public:
    BrentSolver( const std::function< double( double ) >& objectiveFunction_,
                    double regMin_,
                    double regMax_,
                    double tol_,
                    std::size_t maxIter_ )
                    : NonLinearSolver( objectiveFunction_, regMin_, regMax_, tol_, maxIter_ ) {}

    double solve() override;
};



class NonLinearSolverBuilder
{
protected:
    double regMin, regMax;
    double tol = 1e-5;
    std::size_t maxIter = 100;
    std::function< double( double ) > objectiveFunction;
public:
    NonLinearSolverBuilder& setReg( double regMin_, double regMax_ )
    {
        regMin = regMin_;
        regMax = regMax_;
        return *this;
    }
    NonLinearSolverBuilder& setObjectiveFunction( const std::function< double( double ) >& objectiveFunction_ )
    {
        objectiveFunction = objectiveFunction_;
        return *this;
    }
    NonLinearSolverBuilder& setTol( double tol_ )
    {
        tol = tol_;
        return *this;
    }
    NonLinearSolverBuilder& setMaxIter( std::size_t maxIter_ )
    {
        maxIter = maxIter_;
        return *this;
    }
};

class BrentSolverBuilder : public NonLinearSolverBuilder
{
public:
    BrentSolver build()
    {
        return BrentSolver( objectiveFunction, regMin, regMax, tol, maxIter );
    }
};



#endif