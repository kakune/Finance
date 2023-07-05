#ifndef MONTE_CARLO_PROCESS_H
#define MONTE_CARLO_PROCESS_H

#include "generator.h"
#include <cstddef>
#include <vector>
#include <cmath>
#include <random>

class SpotProcess
{
protected:
    AntitheticGaussGenerator Generator;
    std::size_t num_process, num_time;
    double dt, maturity;
    double initial_spot;
    std::vector< std::vector<double> > spots;
    void set_spots()
    {
        spots.resize(num_process);
        for(int i = 0; i < num_process; i++)
        {
            spots.at(i).resize(num_time+1, initial_spot);
        }
    }
public:
    void set_time(std::size_t num_time_, double maturity_){num_time = num_time_; maturity = maturity_; dt = maturity / double(num_time);  Generator.set_num_time(num_time);}
    void set_initial_spot(double initial_spot_){initial_spot = initial_spot_;}
    void set_num_process(std::size_t num_process_){num_process = num_process_;}
    virtual void generate_spots() = 0;
    std::vector< std::vector<double> > get_spots(){return spots;}
};

class BlackScholes : public SpotProcess
{
private:
    double vol, rate;
public:
    void set_vol_rate(double vol_, double rate_){vol = vol_; rate = rate_;}
    void generate_spots();
};

class SABR : public SpotProcess
{
private:
    std::vector< std::vector<double> > vols;
    double volvol, rate, beta, rho, rhom, initial_vol;
public:
    void set_vols()
    {
        vols.resize(num_process);
        for(int i = 0; i < num_process; i++)
        {
            vols.at(i).resize(num_time+1, initial_vol);
        }
    }
    void set_rho_rate(double rho_, double rate_){rho = rho_; rhom = std::sqrt(1.0 - rho*rho); rate = rate_;}
    void set_volvol_beta(double volvol_, double beta_){volvol = volvol_; beta = beta_;}
    void set_initial_vol(double initial_vol_){initial_vol = initial_vol_;}
    void generate_spots();
};

class SV : public SpotProcess
{
private:
    std::vector< std::vector<double> > vols;
    double rate, lambda, b, theta, eta, L, initial_vol;
    double rho, rhom;
public:
    void set_vols()
    {
        vols.resize(num_process);
        for(int i = 0; i < num_process; i++)
        {
            vols.at(i).resize(num_time+1, 1.0);
        }
    }
    void set_rho_rate(double rho_, double rate_){rho = rho_; rhom = std::sqrt(1.0 - rho*rho); rate = rate_;}
    void set_params_spot(double lambda_, double b_, double L_)
    {lambda = lambda_; b = b_; L = L_;}
    void set_params_vol(double theta_, double eta_)
    { theta = theta_; eta = eta_;}
    void generate_spots();
};

#endif