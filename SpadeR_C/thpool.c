#include "thpool.h"


static void tpool_thread(tpool_t tpool);

void tpool_init(tpool_t *tpoolp, int num_worker_threads, int max_queue_size, int do_not_block_when_full, config_t* config, struct common_data *common_data) {
	int i;
	tpool_t tpool;

	//allocate a pool data structure 
	if ((tpool = (tpool_t)malloc(sizeof(struct tpool))) == NULL)
		perror("malloc"), exit(-1);

	//initilize the fields
	tpool->num_threads = num_worker_threads;
	tpool->max_queue_size = max_queue_size;
	tpool->do_not_block_when_full = do_not_block_when_full;
	if ((tpool->threads = (HANDLE *)malloc(sizeof(HANDLE)*num_worker_threads)) == NULL)
		perror("malloc"), exit(-1);
	tpool->cur_queue_size = 0;
	tpool->queue_head = NULL;
	tpool->queue_tail = NULL;
	tpool->queue_closed = 0;
	tpool->shutdown = 0;
	tpool->config = config;
	tpool->common_data = common_data;
	tpool->queue_lock = CreateMutex(NULL,FALSE,NULL);
	tpool->queue_not_full = CreateEvent(
		NULL,               // default security attributes
		TRUE,               // manual-reset event
		FALSE,              // initial state is nonsignaled
		NULL				// object name
	);
	if(tpool->queue_not_full == NULL)
		fprintf(stderr, "pthread_cond_init"), exit(-1);
	tpool->queue_not_empty = CreateEvent(NULL, TRUE, FALSE, NULL);
	if (tpool->queue_not_empty == NULL)
		fprintf(stderr, "pthread_cond_init"), exit(-1);
	tpool->queue_empty = CreateEvent(NULL, TRUE, FALSE, NULL);
	if (tpool->queue_empty == NULL)
		fprintf(stderr, "pthread_cond_init"), exit(-1);
	
	//Allocate memory for the threads and IDs
	tpool->data_mem_th = (data_mem_t*)malloc(num_worker_threads * sizeof(data_mem_t));
	tpool->thIDs = (DWORD*)malloc(sizeof(DWORD)*num_worker_threads);

	//create threads
	for (i = 0; i < num_worker_threads; i++) {
		tpool->threads[i] = CreateThread(
			NULL,                   // default security attributes
			0,                      // use default stack size  
			tpool_thread,			// thread function name
			(void *)(void*)tpool,   // argument to thread function 
			0,                      // use default creation flags 
			&tpool->thIDs[i]);   // returns the thread identifier 
		if(tpool->threads[i] == NULL)
			fprintf(stderr, "CreateThread %d", i), exit(-1);

		printf("ID[%d]: %d\n", i, tpool->thIDs[i]);
		//Init data_mem_t for each thread
		init_data_mem(&tpool->data_mem_th[i], config->read_size, config->decimate_factor, config->fir_length, config->matched_filt_deci_length, config->nfft);
	}

	*tpoolp = tpool;
}

static void tpool_thread(tpool_t tpool) {

	tpool_work_t *my_workp;
	DWORD dwthID = GetCurrentThreadId();
	int id;
	// Search index of the thread id
	for (int i = 0; i < tpool->num_threads; i++)
	{
		if (tpool->thIDs[i] == dwthID) {
			id = i;
			break;
		}
	}
	// Infinite loop of the thread
	for (;;) {
		WaitForSingleObject(tpool->queue_lock, INFINITE);  // no time-out 
		//pthread_mutex_lock(&tpool->queue_lock);
		while (tpool->cur_queue_size == 0 && !tpool->shutdown) {
			//pthread_cond_wait(&tpool->queue_not_empty, &tpool->queue_lock);
			ReleaseMutex(tpool->queue_lock);
			WaitForSingleObject(tpool->queue_not_empty, INFINITE);
			WaitForSingleObject(tpool->queue_lock, INFINITE);
		}
		if (tpool->shutdown) {
			ReleaseMutex(tpool->queue_lock);
			//pthread_mutex_unlock(&tpool->queue_lock);
			ExitThread(0);
			//pthread_exit(NULL);
		}

		my_workp = tpool->queue_head;
		tpool->cur_queue_size--;
		if (tpool->cur_queue_size == 0)
			tpool->queue_head = tpool->queue_tail = NULL;
		else
			tpool->queue_head = my_workp->next;
		if (!tpool->do_not_block_when_full && (tpool->cur_queue_size == tpool->max_queue_size - 1)) {
			SetEvent(tpool->queue_not_full);
			//pthread_cond_broadcast(&tpool->queue_not_full);
		}

		if (tpool->cur_queue_size == 0) {
			SetEvent(tpool->queue_empty);
			//pthread_cond_signal(&tpool->queue_empty);
		}
			
		ReleaseMutex(tpool->queue_lock);
		//pthread_mutex_unlock(&tpool->queue_lock);

		// Execute processing using the correspondig data_mem of this thread
		(*(my_workp->routine))(my_workp->arg1, my_workp->arg2, &tpool->data_mem_th[id], tpool->config, tpool->common_data);
		free(my_workp->arg1);
		free(my_workp);
	}

}

int tpool_add_work(tpool_t tpool, void(*routine)(), void *arg1, size_t arg2) {

	tpool_work_t *workp;
	WaitForSingleObject(tpool->queue_lock, INFINITE);
	//pthread_mutex_lock(&tpool->queue_lock);

	if ((tpool->cur_queue_size == tpool->max_queue_size) && (tpool->do_not_block_when_full)) {
		ReleaseMutex(tpool->queue_lock);
		//pthread_mutex_unlock(&tpool->queue_lock);
		return -1;
	}

	while ((tpool->cur_queue_size == tpool->max_queue_size) && (!(tpool->shutdown || tpool->queue_closed))) {
		ReleaseMutex(tpool->queue_lock);
		WaitForSingleObject(tpool->queue_not_full, INFINITE);
		WaitForSingleObject(tpool->queue_lock, INFINITE);
		//pthread_cond_wait(&tpool->queue_not_full, &tpool->queue_lock);
	}

	//allocate work structure
	workp = (tpool_work_t*)malloc(sizeof(tpool_work_t));
	workp->routine = routine;
	workp->arg1 = arg1;
	workp->arg2 = arg2;
	if (tpool->cur_queue_size == 0) {
		tpool->queue_tail = tpool->queue_head = workp;
		SetEvent(tpool->queue_not_empty);
		//pthread_cond_broadcast(&tpool->queue_not_empty);
	}
	else {
		tpool->queue_tail->next = workp;
		tpool->queue_tail = workp;
	}
	tpool->cur_queue_size++;
	ReleaseMutex(tpool->queue_lock);
	//pthread_mutex_unlock(&tpool->queue_lock);
	return 1;
}



int tpool_destroy(tpool_t tpoolp, int finish) {
	int i, rtn;
	tpool_work_t *cur_nodep;

	WaitForSingleObject(tpoolp->queue_lock, INFINITE);
	//pthread_mutex_lock(&tpoolp->queue_lock);

	if (tpoolp->queue_closed || tpoolp->shutdown) {
		ReleaseMutex(tpoolp->queue_lock);
		//pthread_mutex_unlock(&tpoolp->queue_lock);
		return 0;
	}

	tpoolp->queue_closed = 1;

	if (finish) {
		while (tpoolp->cur_queue_size != 0) {
			ReleaseMutex(tpoolp->queue_lock);
			WaitForSingleObject(tpoolp->queue_empty, INFINITE);
			WaitForSingleObject(tpoolp->queue_lock, INFINITE);
			//pthread_cond_wait(&tpoolp->queue_empty, &tpoolp->queue_lock);
		}
	}

	tpoolp->shutdown = 1;
	ReleaseMutex(tpoolp->queue_lock);
	//pthread_mutex_unlock(&tpoolp->queue_lock);
	SetEvent(tpoolp->queue_not_empty);
	SetEvent(tpoolp->queue_not_full);
	//pthread_cond_broadcast(&tpoolp->queue_not_empty);
	//pthread_cond_broadcast(&tpoolp->queue_not_full);

	WaitForMultipleObjects(tpoolp->num_threads, tpoolp->threads, TRUE, INFINITE);
	/*for (i = 0; i<tpoolp->num_threads; i++) {
		pthread_join(tpoolp->threads[i], NULL);
	}*/

	free(tpoolp->threads);
	while (tpoolp->queue_head != NULL) {
		cur_nodep = tpoolp->queue_head->next;
		tpoolp->queue_head = tpoolp->queue_head->next;
		free(cur_nodep);
	}
	free(tpoolp);
	return 0;
}
