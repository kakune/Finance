#include "process.h"
#include <vector>
#include <iostream>
#include <random>
#include <cmath>
#include <limits>

void BlackScholes::generate_spots()
{
  set_spots();
  double growth = (rate - 0.5 * vol * vol) * dt;
  double factor = vol * std::sqrt(dt);

  for(int i = 0; i < num_process; i++)
  {
    for(int j = 0; j < num_time; j++)
    {
      spots.at(i).at(j+1) = spots.at(i).at(j) * std::exp(growth + factor * Generator());
      // spots.at(i).at(j+1) = spots.at(i).at(j) * (growth + distribution(generator) * factor);
    }
  }
}

void SABR::generate_spots()
{
  set_spots();
  set_vols();
  Generator.set_num_time(num_time*2);
  double growth = 1.0 + rate * dt;
  double sq_dt = std::sqrt(dt);
  double vol_growth = - 0.5 * volvol * volvol * dt;

  for(int i = 0; i < num_process; i++)
  {
    double log_spot = std::log(spots.at(i).at(0));

    for(int j = 0; j < num_time; j++)
    {
      double brown1 = sq_dt * Generator();
      double brown2 = rho * brown1 + rhom * sq_dt * Generator();
      double factor = vols.at(i).at(j) * std::pow(spots.at(i).at(j), beta - 1);

      log_spot += factor * brown1 + (rate - 0.5 * factor * factor) * dt;

      vols.at(i).at(j+1) = vols.at(i).at(j) * std::exp(volvol * brown2 + vol_growth);
      spots.at(i).at(j+1) = std::exp(log_spot);

    }
  }
}

void SV::generate_spots()
{
  set_spots();
  set_vols();
  Generator.set_num_time(num_time*2);

  double discrep = (1.0 - b) * L;
  double sq_dt = std::sqrt(dt);
  double vol_growth = theta - 0.5 * eta * eta;


  for(int i = 0; i < num_process; i++)
  {
    double log_spot = std::log(spots.at(i).at(0));
    double log_vol = std::log(vols.at(i).at(0));

    for(int j = 0; j < num_time; j++)
    {
      double brown1 = sq_dt * Generator();
      double brown2 = rho * brown1 + rhom * sq_dt * Generator();
      double sq_vol = std::sqrt(vols.at(i).at(j));
      double factor = lambda * (b + discrep / spots.at(i).at(j)) * sq_vol;

      log_spot += factor * brown1  - 0.5 * factor * factor * dt + rate * dt;
      log_vol += (eta / sq_vol) * brown2 + (vol_growth / vols.at(i).at(j) - theta) * dt;

      vols.at(i).at(j+1) = std::max(std::exp(log_vol), std::numeric_limits<double>::epsilon());
      spots.at(i).at(j+1) = std::exp(log_spot);
      // spots.at(i).at(j+1) = spots.at(i).at(j) + lambda * (b*spots.at(i).at(j) + discrep) * sq_vol * brown1;
      // vols.at(i).at(j+1) = std::max(vols.at(i).at(j) + theta * (1.0 - vols.at(i).at(j)) * dt + eta * sq_vol * brown2, 0.0);
    }
  }
}