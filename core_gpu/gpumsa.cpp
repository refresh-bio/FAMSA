/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core_gpu/gpumsa.h"

#include "../opencl_utils/hardware/OpenCl.h"
#include "../opencl_utils/hardware/Buffer.h"
#include "../opencl_utils/kernel_repository/KernelRepository.h"
#include "../opencl_utils/kernel_repository/KernelFactory.h"
#include "../opencl_utils/common/mathex.h"
#include "../opencl_utils/common/Log.h"

#include <numeric>
#include <set>
#include <algorithm>
#include <cassert>
#include <atomic>

#undef min
#undef max

using namespace std;

#define SHOW_PROGRESS
//#define COMPARE_WITH_CPU

CGpuFAMSA::CGpuFAMSA(int platform, int device, int pairsPerWave, int pairsPerTask, int threadsPerPair)
	: pairsPerWave(pairsPerWave), pairsPerTask(pairsPerTask), threadsPerPair(threadsPerPair), localBitVectors(true)
{
	// opencl test stuff
	Log::getInstance(Log::LEVEL_DEBUG).enable();
	Log::getInstance(Log::LEVEL_NORMAL).enable();

	cout << "Getting OpenCL devices..." << endl;
	cout << clex::OpenCL::listDevices(clex::OpenCL::ANY_DEVICE) << endl;

	cl = std::shared_ptr<clex::OpenCL>(new clex::OpenCL(clex::OpenCL::ANY_DEVICE, platform, device, false));

	cout << "Following OpenCL device will be used:" << endl << *(cl->mainDevice->info) << endl;

	std::vector<std::string> defines;
	defines.push_back("THREADS_PER_PAIR=" + std::to_string(threadsPerPair));
	defines.push_back("PAIRS_PER_TASK=" + std::to_string(pairsPerTask));

	if (localBitVectors) {
		defines.push_back("LOCAL_BIT_VECTORS");
	}

	#ifdef SEQS_AS_CHARS
		defines.push_back("SEQS_AS_CHARS");
	#endif

	std::vector<int> files;
	files.push_back(KernelRepository::SingleLinkage);

	KernelFactory::instance(cl).setFastMath(true);
	kernelLCS = KernelFactory::instance(cl).create(files, "calculateLCSs", defines);
}


void CGpuFAMSA::UPGMA_CalculateDistances(UPGMA_dist_t * dist_matrix)
{
	cout << "UPGMA guide trees are not supported in a GPU mode.";
	exit(-1);

	::size_t n_seq = sequences.size();

	// calculate max sequence length to be calculated on the GPU
	::size_t maxGpuLength =
		(cl->mainDevice->info->localMemSize - (pairsPerTask * threadsPerPair * (sizeof(cl_int)))) /
		((24 + pairsPerTask) * sizeof(cl_uint));
	maxGpuLength *= 32;

	// process cpu portion
	auto it = std::lower_bound(sequences.begin(), sequences.end(), maxGpuLength, [](const CSequence& c1, size_t val)->bool {
		return c1.length > val;
	});

	::size_t cpuLo = 0;
	::size_t cpuHi = distance(sequences.begin(), it);
	::size_t gpuLo = cpuHi;
	::size_t gpuHi = n_seq;

	std::thread cpuTasksThread;

	if (cpuHi < 0) {
		cout << "CPU range: " << cpuLo << " to " << cpuHi << endl;
		cpuTasksThread = std::thread([this, cpuLo, cpuHi, dist_matrix]()->void {
			std::vector<int> lcs(cpuHi - cpuLo);
			calculateReferenceLCSs(cpuLo, cpuHi, lcs);
			for (int i = cpuLo, j = 0; i < cpuHi; ++i, ++j) {
			//	double indel = (*sequences)[row_id].length + (*sequences)[j * 4 + k].length - 2 * lcs_lens[k];
			//	dist_matrix[i] = 1.0 / (lcs[j] / pow(indel, indel_exp));
			}
		});
	}

	// process waves of given size
//	::size_t lo = gpuLo;
//	while (lo < gpuHi) {
//		::size_t hi = std::min(gpuHi, calculateRangeHi(lo, pairsPerWave));

//		this->calculateLCSs(lo, hi, lcs, postprocessingThread, events[0]);

}

void CGpuFAMSA::SingleLinkage() {
	
#ifdef COMPARE_WITH_CPU
	cout << "Calculating reference guide tree...";
	auto treeCopy = this->guide_tree;
	CFAMSA::SingleLinkage();
	auto refTree = this->guide_tree;
	this->guide_tree = treeCopy;
	cout << "done" << endl;
#endif

	int prefetch_offset = 64 * 2;
	int prefetch_offset_2nd = 128 * 2;

	// prepare OpenCL preprocessing
	preprocess();
	
	int i, j;
	int next;
	::size_t n_seq = sequences.size();

//	vector<int> pi(n_seq);
	vector<int> pi(n_seq + max(prefetch_offset, prefetch_offset_2nd), 0);
	vector<double> lambda(n_seq);
	vector<double> sim_vector(n_seq);
	std::vector<int> lcs(this->pairsPerWave);

	std::thread postprocessingThread;
	std::vector<cl::Event> events(1);

	// calculate max sequence length to be calculated on the GPU
	::size_t maxGpuLength =
		(cl->mainDevice->info->localMemSize - (pairsPerTask * threadsPerPair * (sizeof(cl_int)))) /
		((24 + pairsPerTask) * sizeof(cl_uint));
	maxGpuLength *= 32;

	// process cpu portion
	auto it = std::lower_bound(sequences.begin(), sequences.end(), maxGpuLength, [](const CSequence& c1, size_t val)->bool{
		return c1.length > val;
	});

	::size_t cpuLo = 0;
	::size_t cpuHi = distance(sequences.begin(), it);
	::size_t gpuLo = cpuHi;
	::size_t gpuHi = n_seq;
	
	std::thread cpuTasksThread;
	std::vector<int> cpuLcs(this->calculateRangeElements(cpuLo, cpuHi));
	if (cpuHi > 0) {
		std::string s = "Max sequence length for GPU variant is " + std::to_string(maxGpuLength) + "!";
		throw std::runtime_error(s);
		cout << "CPU range: " << cpuLo << " to " << cpuHi << endl;
		cpuTasksThread = std::thread([this, &cpuLcs, cpuLo, cpuHi, &pi, &lambda, &sim_vector]()->void {
			calculateReferenceLCSs(cpuLo, cpuHi, cpuLcs);
		});
	}

	// process waves of given size
	::size_t lo = gpuLo;
	while (lo < gpuHi) {
		::size_t hi = std::min(gpuHi, calculateRangeHi(lo, pairsPerWave));

		this->calculateLCSs(lo, hi, lcs, postprocessingThread, events[0]);

#ifdef COMPARE_WITH_CPU
		
		cl::WaitForEvents(events);
		
		std::vector<int> refLcs(calculateRangeElements(lo, hi));
		this->calculateReferenceLCSs(lo, hi, refLcs);
		std::vector<int> diffs(refLcs.size());
		
	
		
		cout << "CPU verification: " << endl;
		int t = 0;
		for (i = lo; i < hi; ++i) {
			for (j = 0; j < i; ++j, ++t) {
				int a = lcs[t];
				int b = refLcs[t];
				if (lcs[t] != refLcs[t]) {
					cout << t << ", " 
						<< i << "(" << sequences[i].id << ") and " 
						<< j << "(" << sequences[j].id << "):" 
						<< abs(a - b) << ", " << a << " vs " << b << endl;
				}
			}
		}
	//	cout << endl;
#endif

		postprocessingThread = std::thread([this, &lcs, lo, hi, &pi, &lambda, &sim_vector, &events]()->void {
			cl::WaitForEvents(events);
			this->partialLinkage(lcs, lo, hi, pi, lambda, sim_vector);
		});

		lo = hi;
	}

	if (postprocessingThread.joinable()) {
		postprocessingThread.join();
	}

	if (cpuTasksThread.joinable()) {
		cpuTasksThread.join();
		this->partialLinkage(cpuLcs, cpuLo, cpuHi, pi, lambda, sim_vector);
	}

	// final linkage
	vector<int> elements(n_seq - 1);
	for (i = 0; i < n_seq - 1; ++i)
		elements[i] = i;

	stable_sort(elements.begin(), elements.end(), [&](int x, int y){
		return lambda[x] > lambda[y];
	});

	vector<int> index(n_seq);
	for (i = 0; i < n_seq; ++i)
		index[i] = i;

	for (i = 0; i < n_seq - 1; ++i) {
		j = elements[i];
		next = pi[j];
		guide_tree.push_back(make_pair(index[j], index[next]));
		index[next] = n_seq + i;
	}

// fixme: guide tree comparison
#ifdef COMPARE_WITH_CPU
	bool ok = false;
	if (guide_tree.size() == refTree.size()) {
		std::vector<pair<int, int>> refSet(refTree.begin(), refTree.end());
		std::vector<pair<int, int>> outSet(guide_tree.begin(), guide_tree.end());
		
		// order elements in pair
		auto pairOrderer = [](const pair<int, int>& p)->pair<int, int> {
			return (p.first < p.second) ? p : make_pair(p.second, p.first); };
		std::transform(refSet.begin(), refSet.end(), refSet.begin(), pairOrderer);
		std::transform(outSet.begin(), outSet.end(), outSet.begin(), pairOrderer);
		
		// sort pairs
		auto pairComparer = [](const pair<int, int>& p1, const pair<int, int>& p2)->bool{
			return p1.first < p2.first;};
		std::stable_sort(refSet.begin(), refSet.end(), pairComparer);
		std::stable_sort(outSet.begin(), outSet.end(), pairComparer);

		// get differences
		std::set<pair<int,int>> restRef, restOut;
		std::set_difference(refSet.begin(), refSet.end(), outSet.begin(), outSet.end(), std::inserter(restRef, restRef.begin()));
		std::set_difference(outSet.begin(), outSet.end(), refSet.begin(), refSet.end(), std::inserter(restOut, restOut.begin()));
		ok = (restRef.size() == 0) && (restOut.size() == 0);
	}

	cout << "Comparison with reference guide tree: " << ok << endl;
#endif

}


void CGpuFAMSA::preprocess()
{
	int numSeqs = sequences.size();
	lengths.resize(numSeqs);
	offsets.resize(numSeqs);
	
	// fill lenghts and offsets arrays
	int curOffset = 0;
	for (int i = 0; i != sequences.size(); ++i) {
		auto &s = sequences[i];
		lengths[i] = s.length;
		offsets[i] = curOffset;
		curOffset += mathex::ceilround((int)s.length, 64); // round up for simplicity
	}
	
	// copy raw sequence data and bit vectors
	// bit vectors are of same size as raw data (32 vectors per sequence, each of length sequence / 32)
	rawSequences.resize(curOffset);
	bitVectors.resize(curOffset);

	int vid = 0;
	for (int s = 0; s < sequences.size(); ++s) {
		auto &seq = sequences[s];
		std::copy(seq.data.begin(), seq.data.end(), rawSequences.begin() + offsets[s]);

		for (int i = 0; i < seq.bit_masks.size(); ++i) {
			auto& mask = seq.bit_masks[i];
			for (int j = 0; j < mask.size(); ++j) {
				bitVectors[vid++] = (cl_uint)mask[j];
				bitVectors[vid++] = (cl_uint)(mask[j] >> 32);
				
			}
		}
	}

	assert(vid == bitVectors.size());
	assert(std::all_of(rawSequences.begin(), rawSequences.end(), [](char c)->bool{ return c < 24; }));

	// allocate all OpenCL buffers
	cl_mem_flags ronly_host = CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR;
	buffer_rawSequences = make_shared<clex::Buffer>(*cl, ronly_host, rawSequences.size() * sizeof(type_seq), rawSequences.data(), "buffer_rawSequences");
	buffer_lengths = make_shared<clex::Buffer>(*cl, ronly_host, lengths.size() * sizeof(cl_int), lengths.data(), "buffer_lengths");
	buffer_offsets = make_shared<clex::Buffer>(*cl, ronly_host, offsets.size() * sizeof(cl_int), offsets.data(), "buffer_offsets");
	buffer_bitVectors = make_shared<clex::Buffer>(*cl, ronly_host, bitVectors.size() * sizeof(cl_uint), bitVectors.data(), "buffer_bitVectors");
	 
	int maxTaskCount = (pairsPerWave / pairsPerTask) + numSeqs;
	buffer_tasks = make_shared<clex::Buffer>(*cl, CL_MEM_READ_WRITE, maxTaskCount * sizeof(Task), nullptr, "buffer_tasks");
	buffer_lcs = make_shared<clex::Buffer>(*cl, CL_MEM_WRITE_ONLY, pairsPerWave * sizeof(cl_int), nullptr, "buffer_lcs");

	clCall(kernelLCS->setArg(0, *buffer_rawSequences));
	clCall(kernelLCS->setArg(1, *buffer_lengths));
	clCall(kernelLCS->setArg(2, *buffer_offsets));
	clCall(kernelLCS->setArg(3, *buffer_bitVectors));
	clCall(kernelLCS->setArg(4, *buffer_tasks));
	clCall(kernelLCS->setArg(5, *buffer_lcs));
}


void CGpuFAMSA::calculateLCSs(
	::size_t verticalRangeLo, 
	::size_t verticalRangeHi, 
	std::vector<int>& lcs, 
	std::thread& postprocessingThread,
	cl::Event& finishedEvent)
{
	int numSeqs = sequences.size();
	
	// generate tasks
	::size_t prevCount = (verticalRangeLo * (verticalRangeLo - 1)) / 2;
	::size_t curCount = (verticalRangeHi * (verticalRangeHi - 1)) / 2;

	std::vector<Task> tasks; 
	tasks.reserve((curCount - prevCount) / pairsPerTask + verticalRangeHi - verticalRangeLo);
	
	int taskId = 0;
	for (int i = verticalRangeLo; i < (int) verticalRangeHi; ++i) {
		for (int j = 0; j < i; j += pairsPerTask) {
			Task t;
			t.x = i;
			t.y = j;
			tasks.push_back(t);
		}
	}
	size_t pairsCount = curCount - prevCount;

	//int maxLength = sequences[verticalRangeHi - 1].length;
	int maxLength = sequences[verticalRangeLo].length;

	int vectorLength = mathex::ceildiv(maxLength, 64) * 2;
	::size_t localMemBytes = vectorLength * pairsPerTask * sizeof(cl_uint);

	if (localBitVectors) {
		localMemBytes += vectorLength * 24 * sizeof(cl_uint);
	}

	//clCall(kernelLCS->setArg(6, cl::__local(localMemBytes)));
	clCall(kernelLCS->setArg(6, cl::__local(localMemBytes)));

	cl::NDRange offsetRange(0);
	cl::NDRange localRange(threadsPerPair * pairsPerTask);
	cl::NDRange globalRange(threadsPerPair * pairsPerTask * tasks.size());
	
	cl::Event finished;

	clCall(cl->mainDevice->mainQueue->enqueueWriteBuffer(*buffer_tasks, CL_FALSE, 0, tasks.size() * sizeof(Task), tasks.data()));
	clCall(cl->mainDevice->mainQueue->enqueueNDRangeKernel(*kernelLCS, offsetRange, globalRange, localRange));
	
	if (postprocessingThread.joinable()) {
		postprocessingThread.join();
	}
	
	clCall(cl->mainDevice->mainQueue->enqueueReadBuffer(*buffer_lcs, CL_FALSE, 0, pairsCount * sizeof(cl_int), lcs.data(), nullptr, &finishedEvent));
}

void CGpuFAMSA::partialLinkage(
	const std::vector<int> &lcs,
	::size_t verticalRangeLo,
	::size_t verticalRangeHi,
	vector<int>& pi,
	vector<double>& lambda,
	vector<double>& sim_vector) const
{
	int prefetch_offset = 64 * 2;
	int prefetch_offset_2nd = 128 * 2;

	int i, j;
	int next;
	::size_t n_seq = sequences.size();
	
	double indel_exp = params.indel_exp;

	auto lcs_itr = lcs.begin();
	for (i = (int) verticalRangeLo; i < (int) verticalRangeHi; ++i) {

#ifdef SHOW_PROGRESS
		if (n_seq > 1000 && i % 100 == 0)
		{
			cout << i << " / " << n_seq << "\r";
			fflush(stdout);
		}
#endif

		pi[i] = i;
//		lambda[i] = -infty;
		lambda[i] = -infty_double;

		for (j = 0; j < i; ++j) {
			double lcs_len = (double)(*lcs_itr++);
			double indel = sequences[i].length + sequences[j].length - 2 * lcs_len;
//			sim_vector[j] = lcs_len / (indel * indel);
			sim_vector[j] = lcs_len / pow(indel, indel_exp);
		}

/*		for (j = 0; j < i; ++j) {
			next = pi[j];
			if (lambda[j] > sim_vector[j])
				sim_vector[next] = max(sim_vector[next], sim_vector[j]);
			else {
				sim_vector[next] = max(lambda[j], sim_vector[next]);
				pi[j] = i;
				lambda[j] = sim_vector[j];
			}
		}

		for (j = 0; j < i; ++j) {
			next = pi[j];
			if (lambda[next] >= lambda[j])
				pi[j] = i;
		}*/
	
		auto p_lambda = lambda.begin();
		auto p_sim_vector = sim_vector.begin();
		auto p_pi = pi.begin();

		for (int j = 0; j < i; ++j)
		{
			next = pi[j];

#ifdef _MSC_VER					// Visual C++
			//			_mm_prefetch((const char*) &(*sim_vector)[pi[j + prefetch_offset]], 2);
			_mm_prefetch((const char*)&sim_vector[*(p_pi + prefetch_offset)], 2);
#endif
#ifdef __GNUC__
			//			__builtin_prefetch((&(*sim_vector)[pi[j + prefetch_offset]]), 1, 2);
			__builtin_prefetch((&sim_vector[*(p_pi + prefetch_offset)]), 1, 2);
#endif

			auto &x = sim_vector[next];

			//			if (lambda[j] > (*sim_vector)[j])
			//			if (*p_lambda > *p_sim_vector)
			if (isgreater(*p_lambda, *p_sim_vector))
			{
				//				(*sim_vector)[next] = max((*sim_vector)[next], (*sim_vector)[j]);
				x = max(x, *p_sim_vector);
			}
			else
			{
				//				(*sim_vector)[next] = max(lambda[j], (*sim_vector)[next]);
				x = max(x, *p_lambda);
				//				pi[j] = i;
				*p_pi = i;
				//				lambda[j] = (*sim_vector)[j];
				*p_lambda = *p_sim_vector;
			}

			++p_pi;
			++p_lambda;
			++p_sim_vector;
		}

		p_pi = pi.begin();
		p_lambda = lambda.begin();
		for (int j = 0; j < i; ++j)
		{
#ifdef _MSC_VER					// Visual C++
			//			_mm_prefetch((const char*) &(*sim_vector)[pi[j + prefetch_offset_2nd]], 3);
			_mm_prefetch((const char*)&lambda[*(p_pi + prefetch_offset_2nd)], 0);
#endif
#ifdef __GNUC__
			//			__builtin_prefetch((&(*sim_vector)[pi[j + prefetch_offset_2nd]]), 0, 3);
			__builtin_prefetch((&lambda[*(p_pi + prefetch_offset_2nd)]), 1, 0);
#endif

			//			next = pi[j];
			next = *p_pi;
			//			if (lambda[next] >= lambda[j])
			//			if (lambda[next] >= *p_lambda)
			if (isgreaterequal(lambda[next], *p_lambda))
				//				pi[j] = i;
				*p_pi = i;

			++p_pi;
			++p_lambda;
		}
	}


}


void CGpuFAMSA::calculateLCSs(
	::size_t verticalRangeLo,
	::size_t verticalRangeHi,
	std::vector<int>& lcs)
{
	struct Task {
		int id;
		int i;
		int j;
		Task(int id, int i, int j) : id(id), i(i), j(j) {}
	};

	// generate tasks
	std::vector<Task> tasks;
	int tid = 0;
	for (int i = verticalRangeLo; i < (int) verticalRangeHi; ++i) {
		for (int j = 0; j < i; ++j) {
			//int taskId = calculateRangeElements(verticalRangeLo, i) + j;
			tasks.emplace_back(tid++, i, j);
		}
	}
	
	std::vector<thread> threads(n_threads);
	std::atomic<int> currentTask(0);
	
	for (auto& t : threads) {
		t = thread([&tasks, &lcs, &currentTask, this]()->void {
			int taskId;
			while ((taskId = currentTask.fetch_add(1)) < tasks.size()) {
				const auto& task = tasks[taskId];
				uint32_t res;
				this->lcsbp_classic0.Calculate(&this->sequences[task.i], &this->sequences[task.j], res);
				lcs[task.id] = res;
			}
		});
	}

	for (auto& t : threads) {
		t.join();
	}

}


void CGpuFAMSA::calculateReferenceLCSs(
	::size_t verticalRangeLo,
	::size_t verticalRangeHi,
	std::vector<int>& lcs)
{
	for (int i = verticalRangeLo; i < (int)verticalRangeHi; ++i) {
		for (int j = 0; j < i; ++j) {
			int taskId = calculateRangeElements(verticalRangeLo, i) + j;
			uint32_t res;
			this->lcsbp_classic0.Calculate(&sequences[i], &sequences[j], res);
			lcs[taskId] = res;
		}
	}
}


::size_t CGpuFAMSA::calculateRangeHi(::size_t rangeLo, ::size_t elements) const
{
	float l = rangeLo;
	float n = elements;
	float h = sqrt(l*(l-1) + 2*n);
	
	return (::size_t)h;
}


::size_t CGpuFAMSA::calculateRangeElements(::size_t rangeLo, ::size_t rangeHi) const
{
	::size_t prevCount = (rangeLo * (rangeLo - 1)) / 2;
	::size_t curCount = (rangeHi * (rangeHi - 1)) / 2;

	return curCount - prevCount;
}

