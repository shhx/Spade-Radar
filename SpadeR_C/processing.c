#include "stdafx.h"
#include "processing.h"

#define M_PI		    3.14159265358979323846   // pi
#define I_COMPLEX		 _DCOMPLEX_(0.0, 1.0)

static inline void _Cadd(_Dcomplex* X, const _Dcomplex Y);
static inline _Dcomplex _Cmult(const _Dcomplex X, const _Dcomplex Y);
static inline _Dcomplex _Cmult_constant(const _Dcomplex X, const double Y);
double* hamming(size_t length);

_Dcomplex zero = { 0.0,0.0 };


/*
	Hamming window generator
*/
double* hamming(size_t length) {
	double *out = (double*)malloc(length * sizeof(length));
	for (size_t i = 0; i < length; i++) {
		out[i] = 0.54 - 0.46*cos(2 * M_PI*i / (length - 1));
	}
	return out;
}

/*
	Generate a table of sines and cosines to downconvert a signal.
	fc -> central frequency of the signal
	fs ->	sampling frequency
	length -> length of the signal
*/
void generate_demod_table(double* demod_table, double fc, double fs, size_t length) {
	for (size_t i = 0; i < length; i++) {
		demod_table[i * 2] = 2 * cos(2 * M_PI*(fc / fs)*i) / (1 << 15);
		demod_table[i * 2 + 1] = -2 * sin(2 * M_PI*(fc / fs)*i) / (1 << 15);
	}
}

/*
	Matched filter generator.
	filter -> output filter. Already fliped, windowed and conjugated
	length -> length of the filter
	pulse_width -> pulse width of the signal
	bandwidth -> bandwidth of the signal
	FC -> central frecuency of the output filter
	FS -> sampling frecuency
	real_filter -> select betweeen real or complex filter
	decimate_factor -> factor of decimation of the filter

*/
void generate_matched_filter(_Dcomplex* filter, size_t length, double pulse_width, double bandwidth, double FC, double FS, int real_filter, int decimate_factor) {
	// Frequency slope
	double alpha = bandwidth / (pulse_width*FS*FS);
	_Dcomplex* reversed_filter = (_Dcomplex*)malloc(length * sizeof(_Dcomplex));
	// Signal generation
	for (size_t n = 0; n < length; n++) {
		reversed_filter[n]._Val[0] = cos(2 * M_PI*(alpha*n*n / 2 - (bandwidth / (2 * FS))*n));
		if (real_filter)
			reversed_filter[n]._Val[1] = 0;
		else
			reversed_filter[n]._Val[1] = sin(2 * M_PI*(alpha*n*n / 2 - (bandwidth / (2 * FS))*n));
	}
	// Move to frequency FC and conjugate
	for (size_t n = 0; n < length; n++) {
		_Dcomplex temp = cexp(_Cmulcr(I_COMPLEX, 2 * M_PI*(FC / FS)*n));
		reversed_filter[n] = conj(_Cmult(reversed_filter[n], temp));
	}
	size_t decimated_length = length / decimate_factor;
	double* hamm_win = hamming(decimated_length);
	double acc = 0;
	for (size_t i = 0; i < decimated_length; i++) {
		acc += hamm_win[i];
	}
	// Flip filter, decimate, normalize, apply window and copy to output 
	for (size_t i = 0; i < decimated_length; i++) {
		filter[i] = reversed_filter[length - 1 - i * decimate_factor];
		filter[i] = _Cmulcr(filter[i], hamm_win[i]);
		filter[i] = _Cmulcr(filter[i], decimated_length/acc);
	}
}

/*
	Downconversion using a pre-generated table. Outputs a complex value (IQ)
*/
void demodulate_IQ(const __int16* signal, size_t signal_len, const double* demod_table, _Dcomplex* result) {
	for (size_t i = 0; i < signal_len; i++) {
		result[i]._Val[0] = signal[i] * demod_table[i * 2];
		result[i]._Val[1] = signal[i] * demod_table[i * 2 + 1];
	}
}

double* decimate(const double* signal, size_t length_in, int factor) {
	double* out = (double*)malloc(length_in / factor * sizeof(double));
	for (size_t i = 0; i < length_in / factor; i++) {
		out[i] = signal[i*factor];
	}
	return out;
}

_Dcomplex* decimate_complex(const _Dcomplex* signal, size_t length_in, int factor) {
	_Dcomplex* out = (_Dcomplex*)malloc(length_in / factor * sizeof(_Dcomplex));
	for (size_t i = 0; i < length_in / factor; i++) {
		out[i] = signal[i*factor];
	}
	return out;
}

/*
	Simple complex convolution
*/
void convolve_complex(const _Dcomplex* signal, size_t signal_len, const _Dcomplex* h, size_t h_len, _Dcomplex* result) {
	size_t n;

	for (n = 0; n < signal_len + h_len - 1; n++) {
		size_t kmin, kmax, k;

		result[n] = zero;

		kmin = (n >= h_len - 1) ? n - (h_len - 1) : 0;
		kmax = (n < signal_len - 1) ? n : signal_len - 1;

		for (k = kmin; k <= kmax; k++) {
			_Dcomplex temp = _Cmult(signal[k], h[n - k]);
			_Cadd(result + n, temp);
		}
	}
}

/*
	Simple convolution
*/
void convolve(const double* signal, size_t signal_len, const double* h, size_t h_len, double* result) {
	size_t n;

	for (n = 0; n < signal_len + h_len - 1; ++n) {
		size_t kmin, kmax, k;

		result[n] = 0;

		kmin = (n >= h_len - 1) ? n - (h_len - 1) : 0;
		kmax = (n < signal_len - 1) ? n : signal_len - 1;

		for (k = kmin; k <= kmax; ++k) {
			result[n] += signal[k] * h[n - k];
		}
	}
}

/*
	Filter and decimate using a polyphasic filter. Decimate factor = 4
*/
void filt_polyphasic_4(data_mem_t* data_mem, const _Dcomplex* signal, size_t signal_len, const double* h, size_t h_len, _Dcomplex* result) {
	int D = 4;
	double h0[4], h1[4], h2[4], h3[4];
	for (size_t i = 0; i < h_len - D; i += D) {
		h0[i / 4] = h[i];
		h1[i / 4] = h[i + 1];
		h2[i / 4] = h[i + 2];
		h3[i / 4] = h[i + 3];
	}
	h0[3] = h[12];
	h1[3] = h[13];
	h2[3] = h[14];
	h3[3] = 0;
	_Dcomplex acc0, acc1, acc2, acc3;
	size_t length = (D*D + (D - 1) * 2 + signal_len);
	//_Dcomplex *temp_signal = (_Dcomplex*)malloc(length * sizeof(_Dcomplex));
	memset(data_mem->temp_signal, 0, (D - 1 + signal_len + h_len - 1) * sizeof(_Dcomplex));
	memcpy(data_mem->temp_signal + D - 1, signal, signal_len * sizeof(_Dcomplex));

	size_t kmin, kmax, k;
	size_t n;
	n = 0;

	for (n = 0; n < length / D; n++) {
		acc0 = zero;
		acc1 = zero;
		acc2 = zero;
		acc3 = zero;
		int t = (int)ceil((double)(h_len + 1) / D);
		kmin = (n >= t - 1) ? n - (t - 1) : 0;
		kmax = (n < signal_len - 1) ? n : signal_len - 1;
		_Dcomplex temp;
		for (k = kmin; k <= kmax; k++) {
			int a = k * D;
			temp = _Cmult_constant(data_mem->temp_signal[k*D], h3[n - k]);
			_Cadd(&acc0, temp);
			temp = _Cmult_constant(data_mem->temp_signal[k*D + 1], h2[n - k]);
			_Cadd(&acc1, temp);
			temp = _Cmult_constant(data_mem->temp_signal[k*D + 2], h1[n - k]);
			_Cadd(&acc2, temp);
			temp = _Cmult_constant(data_mem->temp_signal[k*D + 3], h0[n - k]);
			_Cadd(&acc3, temp);
		}
		_Cadd(result + n, acc0);
		_Cadd(result + n, acc1);
		_Cadd(result + n, acc2);
		_Cadd(result + n, acc3);
	}
	//free(temp_signal);
}


void filt_FIR_fast(const _Dcomplex* signal, size_t signal_len, const double* h, size_t h_len, _Dcomplex* result) {
	size_t n;

	for (n = 0; n < signal_len + h_len - 1; n++) {
		size_t kmin, kmax, k;

		result[n] = zero;

		kmin = (n >= h_len - 1) ? n - (h_len - 1) : 0;
		kmax = (n < signal_len - 1) ? n : signal_len - 1;

		for (k = kmin; k <= kmax; k++) {
			_Dcomplex temp = _Cmult_constant(signal[k], h[n - k]);
			_Cadd(result + n, temp);
		}
	}
}

/*
	Complex multiply by constant
*/
static inline _Dcomplex _Cmult_constant(const _Dcomplex X, const double Y) {
	_Dcomplex temp;
	temp._Val[0] = X._Val[0] * Y;
	temp._Val[1] = X._Val[1] * Y;
	return temp;
}

/*
	Complex multiply
*/
static inline _Dcomplex _Cmult(const _Dcomplex X, const _Dcomplex Y) {
	_Dcomplex temp;
	temp._Val[0] = X._Val[0] * Y._Val[0] - Y._Val[1] * X._Val[1];
	temp._Val[1] = Y._Val[0] * X._Val[1] + Y._Val[1] * X._Val[0];
	return temp;
}

/*
	Complex add
*/
static inline void _Cadd(_Dcomplex* X, const _Dcomplex Y) {
	X->_Val[0] += Y._Val[0];
	X->_Val[1] += Y._Val[1];
}

/*
	Covolution using the overlap add method. Parameters are fixed and obtained using Matlab.
*/
void convolve_overlap_add(const _Dcomplex* signal, size_t signal_len, const fftw_complex* h_FFT, size_t h_len, _Dcomplex* result, const fftw_plan p_forward, const fftw_plan p_backward, size_t nfft) {

	size_t P = nfft - h_len + 1;
	fftw_complex *in, *signal_FFT, *temp_result;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nfft);
	signal_FFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nfft);
	temp_result = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nfft);

	size_t i = 0;
	while (i < signal_len) {
		size_t il = (i + P) < signal_len ? (i + P) : signal_len;

		memcpy(in, signal + i, (il - i) * sizeof(fftw_complex));							// Copy chunk of the original signal
		fftw_execute_dft(p_forward, in, signal_FFT);										// FFT		
		for (size_t j = 0; j < nfft; j++) {													// Multiply both responses
			double temp = signal_FFT[j][0];
			signal_FFT[j][0] = signal_FFT[j][0] * h_FFT[j][0] - h_FFT[j][1] * signal_FFT[j][1];
			signal_FFT[j][1] = signal_FFT[j][1] * h_FFT[j][0] + h_FFT[j][1] * temp;
		}
		fftw_execute_dft(p_backward, signal_FFT, temp_result);								// inverse FFT
		size_t k = (i + nfft - 1) < (h_len + signal_len - 1) ? (i + nfft - 1) : (h_len + signal_len - 1); //Maximum length
		for (size_t ii = i; ii < k; ii++) {													// Accumulate results
			result[ii]._Val[0] = result[ii]._Val[0] + temp_result[ii - i][0] * 1 / nfft / h_len;
			result[ii]._Val[1] = result[ii]._Val[1] + temp_result[ii - i][1] * 1 / nfft / h_len;
		}
		i += P;
	}
	fftw_free(in);
	fftw_free(signal_FFT);
	fftw_free(temp_result);
}

void convolve_fft(data_mem_t* data_mem, size_t signal_len, const fftw_complex* h_FFT, size_t h_len, _Dcomplex* result, const fftw_plan p_forward, const fftw_plan p_backward, size_t nfft) {
	
	fftw_execute_dft(p_forward, data_mem->filtered, data_mem->signal_FFT);							// FFT		
	for (size_t j = 0; j < nfft; j++) {																// Multiply both responses
		double temp = data_mem->signal_FFT[j][0];
		data_mem->signal_FFT[j][0] = data_mem->signal_FFT[j][0] * h_FFT[j][0] - h_FFT[j][1] * data_mem->signal_FFT[j][1];
		data_mem->signal_FFT[j][1] = data_mem->signal_FFT[j][1] * h_FFT[j][0] + h_FFT[j][1] * temp;
	}
	fftw_execute_dft(p_backward, data_mem->signal_FFT, data_mem->temp_result);						// inverse FFT
}