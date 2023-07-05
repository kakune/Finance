#ifndef OPTION_OPTION_CORE_H
#define OPTION_OPTION_CORE_H

#include <cstddef>
#include <vector>
#include <memory>

class Option : public std::enable_shared_from_this<Option>
{
protected:
    std::size_t num_path;
    std::vector< double > present_spots;
    std::vector< std::vector< double > > spot_paths;
    double payoff_mean;
public:
    virtual double payoff(double) = 0;
    virtual void calc_payoff_mean() = 0;
    double get_payoff_mean() const { return payoff_mean; }
    virtual void set_present_spots(std::vector< double > present_spots_) {present_spots = present_spots_; num_path = present_spots.size();}
    virtual void set_spot_paths(std::vector< std::vector< double > > spot_paths_)
    {
        spot_paths = spot_paths_;
        num_path = spot_paths.size();
        present_spots.resize(num_path);
        for(int i = 0; i < num_path; i++)
        {
            present_spots.at(i) = spot_paths.at(i).back();
        }
    }
    virtual std::vector< std::shared_ptr<Option> > get_all_Options()
    {
        std::vector< std::shared_ptr<Option> > result = {shared_from_this()};
        return result;
    }
};

class VirtualOption : public Option 
{
private:
    std::shared_ptr<Option> Option1, Option2;
    bool isSum;
public:
    void set_Options(std::shared_ptr<Option> Option1_, std::shared_ptr<Option> Option2_){Option1 = Option1_; Option2 = Option2_;}
    void set_isSum(bool isSum_){isSum = isSum_;}
    double payoff(double) override;
    void calc_payoff_mean() override;
    void set_spot_paths(std::vector< std::vector< double > > spot_paths_) override 
    {
        Option1->set_spot_paths(spot_paths_);
        Option2->set_spot_paths(spot_paths_);
    }
    void set_present_spots(std::vector< double > present_spots_) override
    {
        Option1->set_present_spots(present_spots_);
        Option2->set_present_spots(present_spots_);
    }
    std::vector< std::shared_ptr<Option> > get_all_Options() override
    {
        std::vector< std::shared_ptr<Option> > result;
        std::vector< std::shared_ptr<Option> > Option1s = Option1->get_all_Options();
        std::vector< std::shared_ptr<Option> > Option2s = Option2->get_all_Options();
        result.reserve( Option1s.size() + Option2s.size() ); // preallocate memory
        result.insert( result.end(), Option1s.begin(), Option1s.end() );
        result.insert( result.end(), Option2s.begin(), Option2s.end() );
        return result;
    }
    std::shared_ptr<Option> get_Option1(){return Option1;}
    std::shared_ptr<Option> get_Option2(){return Option2;}
};

class Payoff
{
public:
    virtual double payoff(double) = 0;
};

class Call : public Payoff
{
protected:
    double strike;
public:
    double payoff(double x) override
    {
        return std::max(x - strike, 0.0);
    }
    void set_strike(double strike_){strike = strike_;}
};

class Put : public Payoff
{
protected:
    double strike;
public:
    double payoff(double x) override
    {
        return std::max(strike - x, 0.0);
    }
    void set_strike(double strike_){strike = strike_;}
};

template <class T, class U>
std::shared_ptr<VirtualOption> operator+(std::shared_ptr<T>& obj1, std::shared_ptr<U>& obj2)
{
    static_assert(std::is_base_of<Option, T>::value, "T must be a descendant of A");
    static_assert(std::is_base_of<Option, U>::value, "U must be a descendant of A");
    std::shared_ptr<VirtualOption> resultObj = std::make_shared<VirtualOption>();
    resultObj->set_Options(obj1, obj2);
    resultObj->set_isSum(true);
    return resultObj;
}

template <class T, class U>
std::shared_ptr<VirtualOption> operator-(std::shared_ptr<T>& obj1, std::shared_ptr<U>& obj2)
{
    static_assert(std::is_base_of<Option, T>::value, "T must be a descendant of A");
    static_assert(std::is_base_of<Option, U>::value, "U must be a descendant of A");
    std::shared_ptr<VirtualOption> resultObj = std::make_shared<VirtualOption>();
    resultObj->set_Options(obj1, obj2);
    resultObj->set_isSum(false);
    return resultObj;
}

#endif