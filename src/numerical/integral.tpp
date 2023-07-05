template <typename T>
class OneDimIntegration
{
protected:
    T regMin, regMax;
    std::size_t nStep;
    std::function<T(T)> integrand;
public:
    void setReg(T regMin_, T regMax_){regMin = regMin_; regMax = regMax_;}
    void setNStep(std::size_t nStep_){nStep = nStep_;}
    void setIntegrand(std::function<T(T)> integrand_){integrand = integrand_;}
    T calcTrapzIntegral();
};

template <typename T>
T OneDimIntegration<T>::calcTrapzIntegral()
{
    T result;
    result = 0.5 * integrand(regMin);
    T width = (regMax - regMin) / double(nStep);
    T variable = regMin;
    for(int i_ = 0; i_ < nStep - 1; i_++)
    {
        variable += width;
        result += integrand(variable);
    }
    result += 0.5 * integrand(regMax);
    return result / double(nStep);
}

template <typename T>
T oneDimIntegral(std::function <T(T)> func, T regMin, T regMax, std::size_t nStep = 50)
{
    OneDimIntegration<T> IntObj;
    IntObj.setReg(regMin, regMax);
    IntObj.setNStep(nStep);
    IntObj.setIntegrand(func);
    return IntObj.calcTrapzIntegral();
}