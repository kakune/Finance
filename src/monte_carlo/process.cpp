#include "process.h"
#include <vector>
#include <iostream>
#include <random>
#include <cmath>
#include <limits>

SpotProcess::SpotProcess( double initialSpot_, double maturity_, std::size_t numProcess_, std::size_t numTime_, std::size_t numGeneratorPeriod_ )
  : initialSpot( initialSpot_ ), maturity( maturity_ ), numProcess( numProcess_ ), numTime( numTime_ ), Generator( numGeneratorPeriod_ )
{
    spots.resize(numProcess);
    for(int i = 0; i < numProcess; i++)
    {
        spots.at(i).resize(numTime+1, initialSpot);
    }

    dt = maturity / double(numTime);
}

BlackScholes::BlackScholes( double initialSpot_, double maturity_, std::size_t numProcess_, std::size_t numTime_, std::size_t numGeneratorPeriod_, double vol_, double rate_ )
  : SpotProcess( initialSpot_, maturity_, numProcess_, numTime_, numGeneratorPeriod_ ), vol( vol_ ), rate( rate_ ) {}

StochasticVolatilityProcess::StochasticVolatilityProcess( double initialSpot_, double maturity_, std::size_t numProcess_, std::size_t numTime_, std::size_t numGeneratorPeriod_, double initialVol_, double rho_ )
  : SpotProcess( initialSpot_, maturity_, numProcess_, numTime_, numGeneratorPeriod_ ), initialVol( initialVol_ ), rho( rho_ ), rhoc( std::sqrt(1.0 - rho*rho) )
{
  vols.resize( getNumProcess() );
  for(int i = 0; i < getNumProcess(); i++)
  {
      vols.at(i).resize( getNumTime() + 1, initialVol);
  }
}

SABR::SABR( double initialSpot_, double maturity_, std::size_t numProcess_, std::size_t numTime_, std::size_t numGeneratorPeriod_, double initialVol_, double rho_, double rate_, double beta_, double volvol_ )
  : StochasticVolatilityProcess( initialSpot_, maturity_, numProcess_, numTime_, numGeneratorPeriod_, initialVol_, rho_ ), rate( rate_ ), beta( beta_ ), volvol( volvol_ ) {}

SV::SV( double initialSpot_, double maturity_, std::size_t numProcess_, std::size_t numTime_, std::size_t numGeneratorPeriod_, double rho_, double rate_, double lambda_, double b_, double L_, double theta_, double eta_ )
  : StochasticVolatilityProcess( initialSpot_, maturity_, numProcess_, numTime_, numGeneratorPeriod_, 1.0, rho_ ), rate( rate_ ), lambda( lambda_ ), b( b_ ), L( L_ ), theta( theta_ ), eta( eta_ ) {}

std::vector< std::vector<double> >& BlackScholes::generateSpots()
{
  // setSpots();
  double dt = getDt();
  double growth = (rate - 0.5 * vol * vol) * dt;
  double factor = vol * std::sqrt(dt);

  for(int i = 0; i < getNumProcess(); i++)
  {
    for(int j = 0; j < getNumTime(); j++)
    {
      setSpotsElement(i, j+1, 
      getSpotsElement(i, j) * std::exp(growth + factor * getRandom()));
      // spots.at(i).at(j+1) = spots.at(i).at(j) * (growth + distribution(generator) * factor);
    }
  }
  return getSpots();
}

std::vector< std::vector<double> >& SABR::generateSpots()
{
  // setSpots();
  // setVols();
  // setGeneratorNumTime( getNumTime() * 2 );

  double dt = getDt();
  double sqDt = std::sqrt(dt);

  double rho = getRho();
  double rhoc = getRhoc();

  double growth = 1.0 + rate * dt;
  double volGrowth = - 0.5 * volvol * volvol * dt;

  for(int i = 0; i < getNumProcess(); i++)
  {
    double logSpot = std::log(getSpotsElement( i, 0 ));

    for(int j = 0; j < getNumTime(); j++)
    {
      double brown1 = sqDt * getRandom();
      double brown2 = rho * brown1 + rhoc * sqDt * getRandom();
      double factor = getVolsElement( i, j ) * std::pow(getSpotsElement(i, j), beta - 1);

      logSpot += factor * brown1 + (rate - 0.5 * factor * factor) * dt;

      setVolsElement( i, j+1, getVolsElement( i, j ) * std::exp(volvol * brown2 + volGrowth) );
      setSpotsElement( i, j+1, std::exp(logSpot) );

    }
  }

  return getSpots();
}

std::vector< std::vector<double> >& SV::generateSpots()
{
  // setInitialVol(1.0);
  // setSpots();
  // setVols();
  // setGeneratorNumTime( getNumTime() * 2 );

  double dt = getDt();
  double sqDt = std::sqrt(dt);

  double rho = getRho();
  double rhoc = getRhoc();

  double discrep = (1.0 - b) * L;
  double volGrowth = theta - 0.5 * eta * eta;


  for(int i = 0; i < getNumProcess(); i++)
  {
    double logSpot = std::log(getSpotsElement( i, 0 ));
    double logVol = std::log(getVolsElement( i, 0 ) );

    for(int j = 0; j < getNumTime(); j++)
    {
      double brown1 = sqDt * getRandom();
      double brown2 = rho * brown1 + rhoc * sqDt * getRandom();
      double sqVol = std::sqrt(getVolsElement( i, j ));
      double factor = lambda * (b + discrep / getSpotsElement( i, j )) * sqVol;

      logSpot += factor * brown1  - 0.5 * factor * factor * dt + rate * dt;
      logVol += (eta / sqVol) * brown2 + (volGrowth / getVolsElement( i, j ) - theta) * dt;

      setVolsElement( i, j+1, std::max(std::exp(logVol), std::numeric_limits<double>::epsilon()) );
      setSpotsElement( i, j+1, std::exp(logSpot) );
      // spots.at(i).at(j+1) = spots.at(i).at(j) + lambda * (b*spots.at(i).at(j) + discrep) * sqVol * brown1;
      // vols.at(i).at(j+1) = std::max(vols.at(i).at(j) + theta * (1.0 - vols.at(i).at(j)) * dt + eta * sqVol * brown2, 0.0);
    }
  }
  return getSpots();
}