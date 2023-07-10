#ifndef MONTE_CARLO_GENERATOR_H
#define MONTE_CARLO_GENERATOR_H

#include <random>
#include <vector>

class GaussGenerator
{
private:
    std::mt19937 engine;
    std::normal_distribution<> dist;
    std::size_t numTime;

public:
    GaussGenerator( std::size_t numTime_ = 1 )
        : engine( std::random_device()() ), dist( 0.0, 1.0 ), numTime( numTime_ ) {}
    virtual double operator()() { return dist( engine ); }
    std::size_t getNumTime() { return numTime; }
};

class AntitheticGaussGenerator : public GaussGenerator
{
private:
    std::vector <double> lastValues;
    bool isOdd;
    int indexNext;

public:
    AntitheticGaussGenerator( std::size_t numTime_ )
        : GaussGenerator( numTime_ ), isOdd(true)
    {
        lastValues.resize( getNumTime() );
    }
    double operator()() override;
};

#endif