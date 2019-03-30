#pragma once
#include <Windows.h>

void print_signal(double* signal, size_t len);
void print_complex(_Dcomplex num);
print_complex_fftw(fftw_complex num);
plot_signal(_Dcomplex* signal, size_t size);
void uSleep(int waitTime_us);