template <typename T>
class OneDimIntegration
{
private:
    T regMin, regMax;
    std::size_t nStep;
    std::function< T( T ) > integrand;
public:
    OneDimIntegration( const std::function< T ( T ) >& integrand_, double regMin_, double regMax_, std::size_t nStep_ = 50 )
        : integrand( integrand_ ), regMin( regMin_ ), regMax( regMax_ ), nStep( nStep_ ) {}
    T calcTrapzIntegral();
    T calcSimpsonIntegral();
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
    return result * (regMax - regMin) / double(nStep);
}

template <typename T>
T OneDimIntegration<T>::calcSimpsonIntegral()
{
    T result;
    result = 0.5 * integrand(regMin);
    T width = (regMax - regMin) / double(nStep);
    T variable = regMin;
    bool isOdd = true;
    for(int i_ = 0; i_ < nStep - 1; i_++)
    {
        variable += width;
        result += (1.0 + isOdd) * integrand(variable);
    }
    result += 0.5 * integrand(regMax);
    return result * (2.0 / 3.0) * ((regMax - regMin) / double(nStep));
}

template <typename T>
T oneDimTrapzIntegral(const std::function < T ( T ) >& func, T regMin, T regMax, std::size_t nStep = 50)
{
    OneDimIntegration<T> IntObj( func, regMin, regMax, nStep );
    return IntObj.calcTrapzIntegral();
}

template <typename T>
T oneDimSimpsonIntegral(const std::function < T ( T ) >& func, T regMin, T regMax, std::size_t nStep = 50)
{
    OneDimIntegration<T> IntObj( func, regMin, regMax, nStep );
    return IntObj.calcSimpsonIntegral();
}