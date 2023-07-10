#include "implied_vol.h"

ImpliedVolatility::ImpliedVolatility( double spot_, double strike_, double rate_, double maturity_, double volMin_, double volMax_, double tol_, std::size_t maxIter_ )
    : spot( spot_ ), strike( strike_ ), rate( rate_ ), maturity( maturity_ ), volMin( volMin_ ), volMax( volMax_ ), tol( tol_ ), maxIter( maxIter_ ),
    BSObj( AnalyticalBlackScholesBuilder().setSpot( spot ).setStrike( strike ).setRate( rate ).setMaturity( maturity ).setVol( 0.0 ).build())
{
    SolverBuilder.setReg( volMin, volMax );
    SolverBuilder.setTol( tol );
    SolverBuilder.setMaxIter( maxIter );
}

double ImpliedVolatility::calcImpliedVol( double callValue )
{
    auto objectiveFunction = [&]( double vol )
    {
        return BSObj.payoffCallForVol( vol ) - callValue;
    };
    SolverBuilder.setObjectiveFunction( objectiveFunction );
    BrentSolver SolverObj = SolverBuilder.build();
    return SolverObj.solve();
}