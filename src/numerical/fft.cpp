#include "fft.h"


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

// int main()
// {
//     FFT DFTObj1;
//     std::vector<double> xs{1.0,0.0,1.0,1.0};
//     std::cout << DFTObj1.calcInverseFFT(DFTObj1.calcFFT(xs)).at(1) << std::endl;
//     return 0;
// }