#include "stdafx.h"
#include "utils.h"

void print_signal(double* signal, size_t len) {
	for (size_t i = 0; i < len; i++) {
		printf("%e ", 10*log10(signal[i]));
	}
	printf("\n");
}

void print_complex(_Dcomplex num) {
	printf("%f%+fi ", num._Val[0], num._Val[1]);
}

plot_signal(_Dcomplex* signal, size_t size) {
	for (size_t i = 0; i < size; i++) {
		print_complex(signal[i]);
	}
	printf("\n");
}

print_complex_fftw(fftw_complex num) {
	printf("%f%+fi ", num[0], num[1]);
}

void uSleep(int waitTime_us) {
	__int64 time1 = 0, time2 = 0, freq = 0;

	QueryPerformanceCounter((LARGE_INTEGER *)&time1);
	QueryPerformanceFrequency((LARGE_INTEGER *)&freq);

	do {
		QueryPerformanceCounter((LARGE_INTEGER *)&time2);
	} while ((time2 - time1) < waitTime_us*freq / 1000000);
}