template <typename T>
std::vector< std::complex<double> > calcDFT( const std::vector<T>& vec )
{
    std::size_t nPoint = vec.size();
    std::vector< std::complex<double> > result(nPoint, 0.0);
    std::complex factor(0.0,-2.0*M_PI/double(nPoint));
    for(int i = 0; i < nPoint; i++)
    {
        for(int j = 0; j < nPoint; j++)
        {
            result.at(i) += std::exp(factor * double(i) * double(j)) * vec.at(j);
        }
    }
    return result;
}

template <typename T>
std::vector< std::complex<double> > calcInverseDFT( const std::vector<T>& vec )
{
    std::size_t nPoint = vec.size();
    std::vector<T> vecConj(nPoint);
    for(int i = 0; i < nPoint; i++)
    {
        vecConj.at(i) = std::conj(vec.at(i));
    }
    std::vector< std::complex<double> > result = calcDFT(vecConj);
    for(int i = 0; i < nPoint; i++)
    {
        result.at(i) = std::conj(result.at(i)) / double(nPoint);
    }
    return result;
}

template <typename T>
void bitReverse( std::vector<T>& vec )
{
    for (int i = 0; i < vec.size(); i++)
    {
        int j = 0;
        int bL = 1, bR = vec.size() >> 1;
        for (; bL < vec.size(); bL <<= 1, bR >>= 1) {	// bit reverse
            if ((i & bL) != 0) {
                j |= bR;
            }
        }
        if (i < j) {
            std::swap(vec.at(i), vec.at(j));
        }
    }
}

template <typename T>
std::vector< std::complex<double> > calcFFT(const std::vector<T>& vec)
{
    std::size_t nPoint = vec.size();
    if(!isPow2( nPoint )){throw std::invalid_argument("Length must be pow of 2.");}
    std::vector< std::complex<double> > result( nPoint );
    for(int i = 0; i < nPoint; i++)
    {
        result.at(i) = std::complex<double>(vec.at(i));
    }
    calcFFTStep(nPoint, 0, result);
    bitReverse(result);
    return result;
}

template <typename T>
std::vector< std::complex<double> > calcInverseFFT(const std::vector<T>& vec)
{
    std::size_t nPoint = vec.size();
    std::vector<T> vecConj( nPoint );
    for(int i = 0; i < nPoint; i++)
    {
        vecConj.at(i) = std::conj(vec.at(i));
    }
    std::vector< std::complex<double> > result = calcFFT(vecConj);
    for(int i = 0; i < nPoint; i++)
    {
        result.at(i) = std::conj(result.at(i)) / double( nPoint );
    }
    return result;
}