#include "stdafx.h"
#include "processing.h"
#include "postprocess.h"
#include "utils.h"
#include "init.h"
#include "thpool.h"

#define N							40000
#define FS							250e6
#define FC							75.008e6F
#define PULSE_WIDTH					24e-6
#define BANDWIDTH					52e6//4*2.083e6 //En la version final seran ~50 MHz
#define DECIMATE_FACTOR				4
#define INTEGRATION_LENGTH			20
#define DATA_LEN					11490
#define NFFT						16384
#define OFFSET						0//50

#define FIR_LENGTH					15
#define CFAR_LENGTH					41

#define MAX_THREADS					3
#define MAX_QUEUE					50

tpool_t thread_pool;

//TODO read coefficients from file??
double LPF_filter[FIR_LENGTH] = { 0.0054601,0.010925,-0.0062984,-0.045283,-0.037309,0.091092,0.28855,0.38661,0.28855,0.091092,-0.037309,-0.045283,-0.0062984,0.010925,0.0054601 };
double h_cfar[CFAR_LENGTH] = { 0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250,0.031250 };

double timediff(clock_t t1, clock_t t2) {
	double elapsed;
	elapsed = ((double)t2 - t1) / CLOCKS_PER_SEC * 1000;
	return elapsed;
}

void process_data(__int16 *data, size_t read_size, data_mem_t *data_mem, config_t *config, struct common_data *common_data) {
	for (size_t num_data = 0; num_data < 2; num_data++)	{

		demodulate_IQ(data + OFFSET, read_size, common_data->demod_table, data_mem->demodulated);

		filt_polyphasic_4(data_mem, data_mem->demodulated, read_size, LPF_filter, FIR_LENGTH, data_mem->filtered);
		//plot_signal(data_mem->filtered, data_mem->filtered_size/ DECIMATE_FACTOR);

		//convolve_overlap_add(data_mem->filtered, data_mem->filtered_size / DECIMATE_FACTOR, common_data->matched_filter_fft, 
			//config->matched_filt_deci_length, data_mem->result, common_data->p_forward, common_data->p_backward, NFFT);
		convolve_fft(data_mem, data_mem->filtered_size / DECIMATE_FACTOR, common_data->matched_filter_fft, 
			config->matched_filt_deci_length, data_mem->result, common_data->p_forward, common_data->p_backward, NFFT);

		int result_length = data_mem->result_length - 2 * config->matched_filt_deci_length;
	}
}

int main() {
	FILE* handle;
	errno_t err = fopen_s(&handle, "H:/ch0_test.dat", "rb");
	//errno_t err = fopen_s(&handle, "H:/ch0_prueba16G.dat", "rb");
	//errno_t err = fopen_s(&handle, "C:/Users/Mario/Desktop/Programas ADC/Proyecto Spader/SpaderADC/test_data.bin", "rb");
	if (err == 0) {
		//printf("The file was opened\n");
	} else {
		printf("The file was not opened %d\n", err);
		return 0;
	}

	double* demod_table = (double*)malloc(N * 2 * sizeof(double)); // N values with real and imaginary parts
	generate_demod_table(demod_table, FC, FS, N);

	size_t matched_filter_length = (size_t)(PULSE_WIDTH * FS);
	_Dcomplex* matched_filter = (_Dcomplex*)malloc(matched_filter_length * sizeof(_Dcomplex));
	generate_matched_filter(matched_filter, matched_filter_length, PULSE_WIDTH, BANDWIDTH, 0, FS, 0, DECIMATE_FACTOR);
	size_t matched_filt_deci_length = matched_filter_length / DECIMATE_FACTOR;
	//for (size_t i = 0; i < matched_filter_length/DECIMATE_FACTOR; i++) {
	//	print_complex(matched_filter[i]);
	//}
	fftw_complex *in, *matched_filter_fft;
	fftw_plan p_forward, p_backward;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * NFFT);
	matched_filter_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * NFFT);
	p_forward = fftw_plan_dft_1d(NFFT, in, matched_filter_fft, FFTW_FORWARD, FFTW_MEASURE); // con la opcion FFTW_MEASURE tarda mas pero genera un plan mejor
	p_backward = fftw_plan_dft_1d(NFFT, in, matched_filter_fft, FFTW_BACKWARD, FFTW_MEASURE);
	memset(in, 0, NFFT * sizeof(fftw_complex));
	for (size_t i = 0; i < matched_filt_deci_length; i++) {
		in[i][0] = matched_filter[i]._Val[0];
		in[i][1] = matched_filter[i]._Val[1];
		//print_complex(matched_filter[i]);
	}
	fftw_execute_dft(p_forward, in, matched_filter_fft);
	//plot_signal((_Dcomplex*)matched_filter_fft, NFFT);

	config_t config;
	config.read_size = N;
	config.decimate_factor = DECIMATE_FACTOR;
	config.fir_length = FIR_LENGTH;
	config.matched_filt_deci_length = matched_filt_deci_length;
	config.nfft = NFFT;

	struct common_data common_data;
	common_data.demod_table = demod_table;
	common_data.matched_filter_fft = matched_filter_fft;
	common_data.p_forward = p_forward;
	common_data.p_backward = p_backward;

	tpool_init(&thread_pool, MAX_THREADS, MAX_QUEUE, 0, &config, &common_data);
	
	matrix_t *int_matrix = matrix_data_init(DATA_LEN, INTEGRATION_LENGTH);

	clock_t t1, t2;
	double dif;


	int i = 0;
	t1 = clock();
	while (i<20000 && !feof(handle)) {
		__int16 data_2channels[N * 2];
		size_t read_size = fread(data_2channels, sizeof(__int16), 2*N, handle);
		int TIMES = 100;
		for (size_t ii = 0; ii < TIMES; ii++) {
			i++;
			__int16 *data = (__int16*)malloc(N * sizeof(__int16));
			for (size_t i = 0; i < N; i++) {
				data[i] = data_2channels[2 * i];
			}
			//printf("%d\n", (int)read_size);
		
			read_size /= 2; // ignore channel 2
			read_size -= OFFSET;
			tpool_add_work(thread_pool, process_data, data, read_size);

			uSleep(160);
		}

		//int TIMES = 1;
		//for (size_t i = 0; i < TIMES; i++) {

		//	demodulate_IQ(data + OFFSET, read_size, demod_table, data_mem->demodulated);

		//	filt_polyphasic_4(data_mem, data_mem->demodulated, read_size, LPF_filter, FIR_LENGTH, data_mem->filtered);
		//	//plot_signal(data_mem->filtered, data_mem->filtered_size/ DECIMATE_FACTOR);

		//	//convolve_overlap_add(data_mem->filtered, data_mem->filtered_size / DECIMATE_FACTOR, matched_filter_fft, matched_filt_deci_length, data_mem->result, p_forward, p_backward, NFFT);
		//	convolve_fft(data_mem, data_mem->filtered_size / DECIMATE_FACTOR, matched_filter_fft, matched_filt_deci_length, data_mem->result, p_forward, p_backward, NFFT);
		//	

		//	int result_length = data_mem->result_length - 2 * matched_filt_deci_length;

		//	//plot_signal(data_mem->temp_result, data_mem->result_length);
		//	//if (matrix_data_add(int_matrix, data_mem->result + matched_filt_deci_length, result_length) == int_matrix->intgrt_length) {
		//	//	double *data_integrated = matrix_integrate(int_matrix);
		//	//	//print_signal(data_integrated, result_length);

		//	//	double* cfar_tresh = (double*)malloc((result_length + CFAR_LENGTH - 1) * sizeof(double));
		//	//	size_t* detected = (size_t*)malloc(result_length * sizeof(size_t));
		//	//	memset(detected, 0, result_length * sizeof(size_t));
		//	//	cfar(data_integrated, result_length, h_cfar, CFAR_LENGTH, cfar_tresh);
		//	//	//print_signal(cfar_tresh + CFAR_LENGTH / 2, result_length);

		//	//	size_t n_detected = detect(data_integrated, result_length, cfar_tresh + (CFAR_LENGTH - 1) / 2, detected);
		//	//	//printf("Detected: %lld\n", n_detected);
		//	//	for (size_t i = 0; i < n_detected; i++) {
		//	//		////printf("%d\n", detected[i]);
		//	//		//printf("%.3f m \n", detected[i] * 0.7*3e8 * 4 / 250e6);
		//	//	}
		//	//	free(cfar_tresh);
		//	//	free(data_integrated);
		//	//	free(detected);
		//	//}

		//}
		/*t2 = clock();
		dif = timediff(t1, t2);
		printf("Elapsed: %.6f ms\n", (float)dif / TIMES);*/
	}
	tpool_destroy(thread_pool, 1);
	t2 = clock();
	dif = timediff(t1, t2);
	printf("Elapsed: %f ms\n", dif);
	printf("Number of PRIs = %d\n", i-1);
	printf("%.2f us/frame\n", 1000*dif/(i - 1));
	fclose(handle);

	free(demod_table);
	free(matched_filter);
	free(int_matrix->data);
	free(int_matrix);
	fftw_destroy_plan(p_forward);
	fftw_destroy_plan(p_backward);
	//fftw_free(in); fftw_free(out);

	//while (1);
	return 0;
}
