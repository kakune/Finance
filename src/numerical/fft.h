#ifndef NUMERICAL_FFT_H
#define NUMERICAL_FFT_H
#include <complex>
#include <vector>

template <typename T> std::vector< std::complex<double> > calcDFT( const std::vector<T>& );
template <typename T> std::vector< std::complex<double> > calcInverseDFT( const std::vector<T>& );


bool isPow2(unsigned int x);
template <typename T> void bitReverse(std::vector<T>&);
template <typename T> std::vector< std::complex<double> > calcFFT(const std::vector<T>&);
void calcFFTStep(const int, const int, std::vector< std::complex<double> >&);
template <typename T> std::vector< std::complex<double> > calcInverseFFT(const std::vector<T>&);


#include "fft.tpp"

#endif