#include "fft.h"
#define PI 3.1415926535897932384626433832795028841971
#include <iostream>
#include <exception> 

template <typename T>
std::vector< std::complex<double> > DFT::calcDFT(const std::vector<T>& vec)
{
    nPoint = vec.size();
    std::vector< std::complex<double> > result(nPoint, (0.0, 0.0));
    std::complex factor(0.0,-2.0*PI/nPoint);
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
std::vector< std::complex<double> > DFT::calcInverseDFT(const std::vector<T>& vec)
{
    nPoint = vec.size();
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
void FFT::bitReverse(std::vector<T>& vec)
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
std::vector< std::complex<double> > FFT::calcFFT(const std::vector<T>& vec)
{
    nPoint = vec.size();
    if(!isPow2(nPoint)){throw std::invalid_argument("Length must be pow of 2.");}
    std::vector< std::complex<double> > result(nPoint);
    for(int i = 0; i < nPoint; i++)
    {
        result.at(i) = std::complex<double>(vec.at(i));
    }
    calcFFTStep(nPoint, 0, result);
    bitReverse(result);
    return result;
}
void FFT::calcFFTStep(const int length, const int start, std::vector< std::complex<double> >& result)
{
    const int half = length/2;
    const std::complex<double> factor(0.0, -2.0*PI/double(length));
    if(length > 1)
    {
        for(int p = 0; p < half; p++)
        {
            const std::complex<double> wp = std::exp(factor * double(p));
            const std::complex<double> tempA = result.at(start + p);
            const std::complex<double> tempB = result.at(start + p + half);
            result.at(start + p) = tempA + tempB;
            result.at(start + p + half) = (tempA - tempB) * wp;
        }
        calcFFTStep(half, start, result);
        calcFFTStep(half, start + half, result);
    }   
}

template <typename T>
std::vector< std::complex<double> > FFT::calcInverseFFT(const std::vector<T>& vec)
{
    nPoint = vec.size();
    std::vector<T> vecConj(nPoint);
    for(int i = 0; i < nPoint; i++)
    {
        vecConj.at(i) = std::conj(vec.at(i));
    }
    std::vector< std::complex<double> > result = calcFFT(vecConj);
    for(int i = 0; i < nPoint; i++)
    {
        result.at(i) = std::conj(result.at(i)) / double(nPoint);
    }
    return result;
}


// int main()
// {
//     FFT DFTObj1;
//     std::vector<double> xs{1.0,0.0,1.0,1.0};
//     std::cout << DFTObj1.calcInverseFFT(DFTObj1.calcFFT(xs)).at(1) << std::endl;
//     return 0;
// }