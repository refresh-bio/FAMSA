//#pragma OPENCL EXTENSION cl_intel_printf : enable

#ifdef SEQS_AS_CHARS
	typedef char type_seq;	
#else
	typedef int type_seq;
#endif

#define UNKNOWN_SYMBOL 23

// divides a by b and gets ceil of the result
#define _ceildiv(a,b) (((a) + (b) - 1) / (b))

// perform ceil rounding to the closest multiplicity of b
#define _ceilround(a,b) (((a) + (b) - 1) / (b) * (b))

typedef struct {
	int x;
	int y;
} Task;


__kernel __attribute__((reqd_work_group_size(THREADS_PER_PAIR * PAIRS_PER_TASK, 1, 1)))
void calculateLCSs(
	global		const	type_seq*		sequences,
	global		const	int*			lengths,
	global		const	int*			offsets,
	global		const	unsigned int*	bitVectors,
	global		const	Task*			tasks,
	global				int*			lcs,
	local				unsigned int*	localBuffer)
{
	#define task_id		get_group_id(0)
	#define local_size  get_local_size(0)
	
	int pair_id = get_local_id(0) / THREADS_PER_PAIR;
	int pair_offset = get_local_id(0) % THREADS_PER_PAIR;
	int j = pair_offset;
	int x = tasks[task_id].x;
	int y = tasks[task_id].y + pair_id; // task gives only starting index

	int bit_length_x = _ceildiv(lengths[x], 32);
	int bv_offset = (bit_length_x % 2 == 0) ? bit_length_x : bit_length_x + 1 ; // round up to even number
	int length_y = (y < x) ? lengths[y] : 0;
	
	local unsigned int* X = localBuffer + pair_id * bv_offset;
	bitVectors += offsets[x];
	global const type_seq *seq_y = (y < x) ? (sequences + offsets[y]) : 0;
	
	#ifdef LOCAL_BIT_VECTORS
		local unsigned int* bv_x = localBuffer + bv_offset * PAIRS_PER_TASK;	
		for (int k = get_local_id(0); k < bv_offset * 24; k += local_size) {
			bv_x[k] = bitVectors[k];
		}
	#else
		global const unsigned int* bv_x = bitVectors;
	#endif

	unsigned int V, V2, tB, sB = (unsigned int)0;
	int end_i = max(THREADS_PER_PAIR, bit_length_x);
	
	local int reqIters;
	// int last_pair = min(PAIRS_PER_TASK - 1, (x - 1) - (y - pair_id));
	// if (pair_id == last_pair && j == 0) {

	// first pair requires most iterations
	if (get_local_id(0) == 0) {
		reqIters = (length_y / THREADS_PER_PAIR + 1) * end_i + (length_y % THREADS_PER_PAIR);
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	int itersCount = reqIters;
	
	// process all pairs
	int i = -j;
	for (int iter = 0; iter < itersCount; ++iter) {
		
		bool legal = (i >= 0 && i < bit_length_x && j < length_y);
		int c = legal ? seq_y[j] : UNKNOWN_SYMBOL;  
		
		if (c != UNKNOWN_SYMBOL) {
			V = (j == 0) ? ~(unsigned int)0 : X[i];		// initialisation condition
		} 

		barrier(CLK_LOCAL_MEM_FENCE);

		if (c != UNKNOWN_SYMBOL) { // unknown symbol or illegal index
			tB = V & bv_x[bv_offset * c + i];
			V2 = V + tB + sB;
			
	//		unsigned int sBprev = sB;
			sB = V2 < V;
			X[i] = V2 | (V - tB);

	//		printf("thread %d, GPU %d, bitv=%08x, V=%08x, tb=%08x, sBprev=%08x, V2=%08x, sB=%08x, X[%d]=%08x\n", 
	//			pair_offset, j, bv_x[bv_offset * c + i], V, tB, sBprev, V2, sB, i, X[i]);
			
		}

		barrier(CLK_LOCAL_MEM_FENCE);

		if (++i == end_i) {
			j += THREADS_PER_PAIR;
			i = 0;
			sB = (unsigned int)0;
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);
	

	int q = 0;
	for (int i = pair_offset; i < bit_length_x; i += THREADS_PER_PAIR) {
		#if __OPENCL_VERSION__  >= 120
			q += popcount(~X[i]);
		#else
			for(V = ~X[i]; V; V &= V-1) { ++q; }
		#endif
	}
	
	local int res[THREADS_PER_PAIR * PAIRS_PER_TASK];
	res[get_local_id(0)] = q;
	barrier(CLK_LOCAL_MEM_FENCE);

	if ((pair_offset == 0) && (y < x)) { // first thread for each pair
		q = 0;
		for(int i = 0; i < THREADS_PER_PAIR; ++i) {
			q += res[get_local_id(0) + i];		
		}

		int global_pair_id = (x * (x - 1) - tasks[0].x * (tasks[0].x - 1)) / 2 + y;
		lcs[global_pair_id] = q;
	}
}

