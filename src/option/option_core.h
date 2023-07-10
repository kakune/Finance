#ifndef OPTION_OPTION_CORE_H
#define OPTION_OPTION_CORE_H

#include <cstddef>
#include <vector>
#include <memory>

class Option : public std::enable_shared_from_this<Option>
{
private:
    std::size_t numPath;
    std::vector< double > presentSpots;
    std::vector< std::vector< double > > spotPaths;
    double payoffMean;
public:
    virtual double payoff(double) = 0;
    virtual void calcPayoffMean() = 0;
    void setPayoffMean( double payoffMean_ ) { payoffMean = payoffMean_; }

    void setNumPath( std::size_t numPath_ ) { numPath = numPath_; }
    virtual void setPresentSpots( const std::vector< double >& presentSpots_)
    {
        presentSpots = presentSpots_;
        setNumPath( presentSpots.size() );
    }
    virtual void setSpotPaths( const std::vector< std::vector< double > >& spotPaths_)
    {
        spotPaths = spotPaths_;
        setNumPath( spotPaths.size() );
        presentSpots.resize(numPath);
        for(int i = 0; i < numPath; i++)
        {
            presentSpots.at(i) = spotPaths.at(i).back();
        }
    }

    const std::vector< std::vector< double > >& getSpotPaths() { return spotPaths; }
    const std::vector< double >& getPresentSpots() { return presentSpots; }

    double getPayoffMean() const { return payoffMean; }
    std::size_t getNumPath() const { return numPath; }
    virtual std::vector< std::shared_ptr<Option> > getAllOptions()
    {
        std::vector< std::shared_ptr<Option> > result = { shared_from_this() };
        return result;
    }
};

class VirtualOption : public Option 
{
private:
    std::shared_ptr<Option> Option1, Option2;
    bool isSum;
public:
    void setOptions(std::shared_ptr<Option> Option1_, std::shared_ptr<Option> Option2_){Option1 = Option1_; Option2 = Option2_;}
    void setIsSum(bool isSum_){isSum = isSum_;}
    double payoff(double) override;
    void calcPayoffMean() override;
    void setSpotPaths( const std::vector< std::vector< double > >& spotPaths_) override 
    {
        Option1->setSpotPaths(spotPaths_);
        Option2->setSpotPaths(spotPaths_);
    }
    void setPresentSpots( const std::vector< double >& presentSpots_) override
    {
        Option1->setPresentSpots(presentSpots_);
        Option2->setPresentSpots(presentSpots_);
    }
    std::vector< std::shared_ptr<Option> > getAllOptions() override
    {
        std::vector< std::shared_ptr<Option> > result;
        std::vector< std::shared_ptr<Option> > Option1s = Option1->getAllOptions();
        std::vector< std::shared_ptr<Option> > Option2s = Option2->getAllOptions();
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
private:
    double strike;
public:
    double payoff(double x) override
    {
        return std::max( x - strike , 0.0);
    }
    void setStrike( double strike_ ) { strike = strike_; }
    double getStrike() { return strike; }
};

class Put : public Payoff
{
protected:
    double strike;
public:
    double payoff(double x) override
    {
        return std::max( strike - x , 0.0);
    }
    void setStrike( double strike_ ) { strike = strike_; }
    double getStrike() { return strike; }
};

template <class T, class U>
const std::shared_ptr<VirtualOption> operator+(const std::shared_ptr<T>& obj1, const std::shared_ptr<U>& obj2)
{
    static_assert(std::is_base_of<Option, T>::value, "T must be a descendant of A");
    static_assert(std::is_base_of<Option, U>::value, "U must be a descendant of A");
    std::shared_ptr<VirtualOption> resultObj = std::make_shared<VirtualOption>();
    resultObj->setOptions(obj1, obj2);
    resultObj->setIsSum(true);
    return resultObj;
}

template <class T, class U>
const std::shared_ptr<VirtualOption> operator-(const std::shared_ptr<T>& obj1, const std::shared_ptr<U>& obj2)
{
    static_assert(std::is_base_of<Option, T>::value, "T must be a descendant of A");
    static_assert(std::is_base_of<Option, U>::value, "U must be a descendant of A");
    std::shared_ptr<VirtualOption> resultObj = std::make_shared<VirtualOption>();
    resultObj->setOptions(obj1, obj2);
    resultObj->setIsSum(false);
    return resultObj;
}

#endif