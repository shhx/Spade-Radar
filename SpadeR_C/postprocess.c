#include "stdafx.h"
#include "postprocess.h"


matrix_t* matrix_data_init(size_t data_length, size_t integration_length) {
	matrix_t *m = (matrix_t*)malloc(sizeof(matrix_t));
	m->data_length = data_length;
	m->intgrt_length = integration_length;
	m->index = 0;
	m->data = (_Dcomplex*)malloc(data_length*integration_length*(sizeof(_Dcomplex)));
	memset(m->data, 0, data_length*integration_length*(sizeof(_Dcomplex)));
	return m;
}

size_t matrix_data_add(matrix_t *m, const _Dcomplex *data, size_t data_length) {
	memcpy(m->data + (m->index*m->data_length), data, data_length * sizeof(_Dcomplex));
	m->index++;
	return m->index;
}

double* matrix_integrate(matrix_t *m) {
	double *result = (double*)malloc(m->data_length * sizeof(double));
	memset(result, 0, m->data_length * sizeof(double));
	for (size_t i = 0; i < m->data_length; i++) {			// Integrate the matrix
		for (size_t j = 0; j < m->intgrt_length; j++) {
			double abs = cabs(m->data[j*m->data_length + i]);
			result[i] += abs * abs * 1 / m->intgrt_length;
		}
	}
	m->index = 0;
	return result;
}

void cfar(double* signal, size_t signal_len, double* h_cfar, size_t cfar_length, double* result) {

	size_t n;
	int l_cfar = 32;
	double Pfa = 1e-6;
	double alpha = l_cfar * (pow(Pfa, -1.0 / l_cfar) - 1);

	for (n = 0; n < signal_len + cfar_length - 1; n++) {
		size_t kmin, kmax, k;
		result[n] = 0;
		kmin = (n >= cfar_length - 1) ? n - (cfar_length - 1) : 0;
		kmax = (n < signal_len - 1) ? n : signal_len - 1;
		for (k = kmin; k <= kmax; k++) {
			double temp = fabs(signal[k]);
			result[n] += signal[k] * h_cfar[n - k];
		}
		result[n] *= alpha;
	}
}

size_t detect(double* signal, size_t signal_len, double* cfar, size_t* detected) {
	size_t n = 0;
	for (size_t i = 0; i < signal_len; i++) {
		double t1 = signal[i];
		double t2 = cfar[i];
		if (t1 > t2) {
			detected[n++] = i;
		}
	}
	return n;
}