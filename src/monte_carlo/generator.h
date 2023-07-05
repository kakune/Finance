#ifndef GENERATOR_MODULE_H
#define GENERATOR_MODULE_H

#include <random>
#include <vector>

class GaussGenerator
{
protected:
    std::mt19937 engine;
    std::normal_distribution<> dist;
    std::size_t num_time;

public:
    GaussGenerator()
        : engine(std::random_device()()), dist(0.0, 1.0), num_time(1) {}
    virtual void set_num_time(size_t num_time_) { num_time = num_time_; }
    virtual double operator()() { return dist(engine); }
};

class AntitheticGaussGenerator : public GaussGenerator
{
protected:
    std::vector <double> last_values;
    bool is_odd_call;
    int index_next;

public:
    AntitheticGaussGenerator()
        : GaussGenerator(), is_odd_call(true){}
    void set_num_time(size_t num_time_) override { num_time = num_time_; last_values.resize(num_time); }
    double operator()() override;
};

#endif