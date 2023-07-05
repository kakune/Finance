#ifndef FFT_MODULE_H
#define FFT_MODULE_H
#include <complex>
#include <vector>

class DFT
{
protected:
    std::size_t nPoint;
public:
    template <typename T> std::vector< std::complex<double> > calcDFT(const std::vector<T>&);
    template <typename T> std::vector< std::complex<double> > calcInverseDFT(const std::vector<T>&);
};

class FFT : public DFT
{
public:
    bool isPow2(unsigned int x)
    {
        if (x == 0) {return false;}
        return (x & (x - 1)) == 0;
    }
    template <typename T> void bitReverse(std::vector<T>&);
    template <typename T> std::vector< std::complex<double> > calcFFT(const std::vector<T>&);
    void calcFFTStep(const int, const int, std::vector< std::complex<double> >&);
    template <typename T> std::vector< std::complex<double> > calcInverseFFT(const std::vector<T>&);

};

#endif