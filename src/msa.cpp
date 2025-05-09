/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "msa.h"
#include "./tree/GuideTree.h"
#include "./tree/SingleLinkage.h"
#include "./tree/FastTree.h"
#include "./tree/MSTPrim.h"
#include "./tree/UPGMA.h"
#include "./tree/NeighborJoining.h"
#include "./tree/DistanceCalculator.h"
#include "./core/io_service.h"
#include "./utils/log.h"

#undef min
#undef max

#include <algorithm>
#include <set>
#include <random>
#include <map>
#include <list>
#include <stack>
#include <thread>
#include <memory>

#include <iostream>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <sstream>

using namespace std;

#define SHOW_PROGRESS

// *******************************************************************
CFAMSA::CFAMSA(CParams& _params) :
	params(_params), 
	instruction_set(params.instruction_set), 
	final_profile(nullptr),
	atp(_params.n_threads, _params.n_threads)
{
	initScoreMatrix();
}

// *******************************************************************
CFAMSA::~CFAMSA()
{
	delete final_profile;
}

// *******************************************************************
void CFAMSA::initScoreMatrix()
{
	score_matrix.resize(NO_AMINOACIDS);

	auto& sm_matrix = ScoringMatrices::get_matrix(params.matrix_type);

#ifdef HUGE_ALIGNMENTS
	for(int i = 0; i < NO_AMINOACIDS; ++i)
	{
		score_vector.emplace_back(sm_matrix[i][i]);
		for(int j = 0; j < NO_AMINOACIDS; ++j)
			score_matrix[i].emplace_back(sm_matrix[i][j]);
	}
#else
	for(int i = 0; i < NO_AMINOACIDS; ++i)
	{
		score_vector.emplace_back((score_t) round(sm_matrix[i][i] * cost_cast_factor));
		for(int j = 0; j < NO_AMINOACIDS; ++j)
			score_matrix[i].emplace_back((score_t) round(sm_matrix[i][j] * cost_cast_factor));
	}
#endif
}

// *******************************************************************
void CFAMSA::adjustParams(int n_seqs)
{

	if ((params.gt_heuristic != GT::None) && (n_seqs < params.sample_size)) {
		params.gt_heuristic = GT::None;
	}

	if(params.enable_gap_rescaling)
	{
		double gap_scaler = log2(n_seqs / (double) params.scaler_log);
		if(n_seqs < (int) params.scaler_log)
			gap_scaler = 1.0;
		else
			gap_scaler = 1.0 + (gap_scaler / params.scaler_div);

		params.gap_ext		 = (score_t) (params.gap_ext       * gap_scaler);
		params.gap_open	 	 = (score_t) (params.gap_open      * gap_scaler);
		params.gap_term_ext  = (score_t) (params.gap_term_ext  * gap_scaler);
		params.gap_term_open = (score_t) (params.gap_term_open * gap_scaler);
	}

	params.score_matrix = score_matrix;
	params.score_vector = score_vector;
}

// *******************************************************************
std::shared_ptr<AbstractTreeGenerator> CFAMSA::createTreeGenerator(const CParams& params) {
	
	std::shared_ptr<AbstractTreeGenerator> gen = nullptr;
	
	// distance only mode
	if (params.export_distances) {
		LOG_VERBOSE << "Calculating distances and storing in: " << params.output_file_name;

		std::shared_ptr<AbstractTreeGenerator> calculator;

		if (params.distance == Distance::indel_div_lcs) {
			gen = make_shared<DistanceCalculator<Distance::indel_div_lcs>>(
				params.n_threads, params.instruction_set, 
				params.output_file_name, params.generate_square_matrix, params.calculate_pid);

		}
		else if (params.distance == Distance::indel075_div_lcs) {
			gen = make_shared<DistanceCalculator<Distance::indel075_div_lcs>>(
				params.n_threads, params.instruction_set, 
				params.output_file_name, params.generate_square_matrix, params.calculate_pid);
		}
	
		return gen;
	}

	if (params.gt_method == GT::SLINK || (params.gt_heuristic != GT::None && params.gt_method == GT::MST_Prim)) {
		// single linkage (use also when MST used as a partial generator)
		if (params.distance == Distance::indel_div_lcs) {
			gen = make_shared<SingleLinkage<Distance::indel_div_lcs>>(params.n_threads, params.instruction_set);
		}
		else if (params.distance == Distance::indel075_div_lcs) {
			gen = make_shared<SingleLinkage<Distance::indel075_div_lcs>>(params.n_threads, params.instruction_set);
		}
	}
	else if (params.gt_method == GT::MST_Prim) {
		// Prim's minimum spanning tree 
		if (params.distance == Distance::indel_div_lcs) {
			gen = make_shared<MSTPrim<Distance::indel_div_lcs>>(params.n_threads, params.instruction_set);
		}
		else if (params.distance == Distance::indel075_div_lcs) {
			gen = make_shared<MSTPrim<Distance::indel075_div_lcs>>(params.n_threads, params.instruction_set);
		}
	}
	else if (params.gt_method == GT::UPGMA || params.gt_method == GT::UPGMA_modified) {
		// UPGMA
		if (params.distance == Distance::indel_div_lcs) {
			gen = make_shared<UPGMA<Distance::indel_div_lcs>>(params.n_threads, params.instruction_set,
				(params.gt_method == GT::UPGMA_modified));
		}
		else if (params.distance == Distance::indel075_div_lcs) {
			gen = make_shared<UPGMA<Distance::indel075_div_lcs>>(params.n_threads, params.instruction_set,
				(params.gt_method == GT::UPGMA_modified));
		}
	}
	else if (params.gt_method == GT::NJ) {
		// neighbour joining
		if (params.distance == Distance::indel_div_lcs) {
			gen = make_shared<NeighborJoining<Distance::indel_div_lcs>>(params.n_threads, params.instruction_set);
		}
		else if (params.distance == Distance::indel075_div_lcs) {
			gen = make_shared<NeighborJoining<Distance::indel075_div_lcs>>(params.n_threads, params.instruction_set);
		}
	} else {

		throw std::runtime_error("Error: Illegal guide tree method.");
	}

	// check if heuristic was specified
	if (params.gt_heuristic != GT::None) {

		// part-tree versus medoid-tree
		shared_ptr<IClustering> clustering = (params.gt_heuristic == GT::PartTree)
			? nullptr
			: make_shared<CLARANS>(params.cluster_fraction, params.cluster_iters);

		// local seed dumper
		class SeedDumper : public IFastTreeObserver {
		protected:
			std::ofstream ofs;
		public:
			SeedDumper(const std::string& fname) {
				ofs.open(fname);
			}

			void notifySeedsSelected(const std::vector<CSequence*>& seeds, int depth) override {
				if (depth == 0) {
					for (const auto s : seeds) {
						ofs << s->id.substr(1) << std::endl;
					}
				}
			}
		};

		// verify distance measure
		if (params.distance == Distance::indel_div_lcs) {
			auto ft = make_shared<FastTree<Distance::indel_div_lcs>>(
				params.n_threads,
				params.instruction_set,
				dynamic_pointer_cast<IPartialGenerator>(gen),
				params.subtree_size,
				clustering,
				params.sample_size);

			gen = ft;

			if (!params.seed_file_name.empty()) {
				ft->registerObserver(make_shared<SeedDumper>(params.seed_file_name));
			}
		}
		else if (params.distance == Distance::indel075_div_lcs) {
			auto ft = make_shared<FastTree<Distance::indel075_div_lcs>>(
				params.n_threads,
				params.instruction_set,
				dynamic_pointer_cast<IPartialGenerator>(gen),
				params.subtree_size,
				clustering,
				params.sample_size);

			gen = ft;

			if (!params.seed_file_name.empty()) {
				ft->registerObserver(make_shared<SeedDumper>(params.seed_file_name));
			}
		}

		
	}

	return gen;
}

// *******************************************************************
void CFAMSA::sortAndExtendSequences(std::vector<CSequence>& sequences) {

	// create vector of pointers just for sorting
	std::vector<CSequence*> seq_ptrs(sequences.size());
	std::transform(sequences.begin(), sequences.end(), seq_ptrs.begin(), [](CSequence& s) { return &s; });

	auto comparer = [](const CSequence* a, const CSequence* b)->bool {
		return a->length > b->length ||
			(a->length == b->length && std::lexicographical_compare(a->data, a->data + a->data_size, b->data, b->data + b->data_size));
	};

	std::stable_sort(seq_ptrs.begin(), seq_ptrs.end(), comparer);
	uint32_t max_seq_len = seq_ptrs[0]->length;

	std::vector<CSequence> output;
	output.reserve(sequences.size());

	// mark already allocated memory area to free
	auto mma = sequences.front().get_mma();
	if (mma) {
		mma->freeze();
	}

	// put sequences in the ordered collection with length extension
	for (int i = 0; i < seq_ptrs.size(); ++i) {
		output.emplace_back(std::move(sequences[seq_ptrs[i]->original_no]));
		output.back().DataResize(max_seq_len, UNKNOWN_SYMBOL);
	}
	output.swap(sequences);

	// free marked area
	if (mma) {
		mma->release_freezed();
	}
}

// *******************************************************************
void CFAMSA::extendSequences(std::vector<CSequence>& sequences) {
	
	// Temporarily resize the sequences by adding unknown symbols (just to simplify the pairwise comparison)
	auto& max_seq = *max_element(sequences.begin(), sequences.end(),
		[](const CSequence& x, const CSequence& y) { return x.length < y.length; });

	uint32_t max_seq_len = max_seq.length;

	auto mma = sequences.front().get_mma();

	if (mma) {
		mma->freeze();
	}

	int n_seqs = (int)sequences.size();
	for (int i = 0; i < n_seqs; ++i) {
		sequences[i].DataResize(max_seq_len, UNKNOWN_SYMBOL);
		
		/*
		
		if (i == representative.original_no) {
			// non-duplicate -> extend
			sequences[i].DataResize(max_seq_len, UNKNOWN_SYMBOL);
		}
		else {
			// duplicate -> just reallocate to free memory later
			sequences[i].DataResize(max_seq_len, UNKNOWN_SYMBOL);
			//sequences[i].DataResize(sequences[i].length, UNKNOWN_SYMBOL);
		}
		*/
	}

	if (mma) {
		mma->release_freezed();
	}
}

// *******************************************************************
void CFAMSA::shrinkSequences(std::vector<CSequence>& sequences) {
	auto mma = sequences.front().get_mma();

	if (mma) {
		mma->freeze();
	}

	int n_seqs = (int)sequences.size();
	for (int i = 0; i < n_seqs; ++i) {
		sequences[i].DataResize(sequences[i].length, UNKNOWN_SYMBOL);
	}

	if (mma) {
		mma->release_freezed();
	}
}

// *******************************************************************
void CFAMSA::removeDuplicates(std::vector<CSequence*>& sorted_seqs, std::vector<int>& original2sorted) {
	
	auto eq_comparer = [](const CSequence* a, const CSequence* b)->bool {
		return a->length == b->length && std::equal(a->data, a->data + a->length, b->data);
	};

	// update original2sorted mappings to take into account duplicates
	int cur_sorted_index = 0;
	for (int i = 1; i < (int)sorted_seqs.size(); ++i) {
		if (!eq_comparer(sorted_seqs[i], sorted_seqs[i - 1])) {
			++cur_sorted_index;
		}

		original2sorted[i] = cur_sorted_index;
	}

	auto newend = std::unique(sorted_seqs.begin(), sorted_seqs.end(), eq_comparer);
	sorted_seqs.erase(newend, sorted_seqs.end());
}

// *******************************************************************
// Compute Alignment according to guide tree
CProfile* CFAMSA::ComputeAlignment(std::vector<CGappedSequence*>& gapped_sequences, tree_structure& guide_tree)
{
	CProfile* profile = new CProfile(&params, &atp);

	profile->Clear();

	CProfileQueue pq(&gapped_sequences, &profiles, &guide_tree, params.n_threads);

	vector<thread *> workers(params.n_threads, nullptr);

	uint32_t computed_prof = 0;
	mutex mtx;

	size_t ref_thr = params.thr_internal_refinement;

	for(uint32_t i = 0; i < params.n_threads; ++i)
		workers[i] = new thread([&]{
			CGappedSequence *gs;
			CProfile *prof1;
			CProfile *prof2;
			CProfile *prof_sol;
			size_t prof_id;
			//bool only_task;
			uint32_t no_threads;
			uint32_t no_rows_per_box;

//			while(pq.GetTask(prof_id, gs, prof1, prof2, only_task))
			while(pq.GetTask(prof_id, gs, prof1, prof2, no_threads, no_rows_per_box))
			{
				if(gs != nullptr)
				{
					prof_sol = new CProfile(*gs, &params);
				}
				else
				{
					// Internal refinement
					if (prof1->Size() + prof2->Size() > ref_thr)
					{
						if (prof1->Size() <= ref_thr && prof1->Size() > 2)
							RefineAlignment(prof1, no_threads);
						if (prof2->Size() <= ref_thr && prof2->Size() > 2)
							RefineAlignment(prof2, no_threads);
					}

					prof_sol = new CProfile(prof1, prof2, &params, no_threads, no_rows_per_box, &atp);

					delete prof1;
					delete prof2;
				}

				pq.AddSolution(prof_id, prof_sol);

				if (params.very_verbose_mode)
				{
					lock_guard<mutex> lck(mtx);
					computed_prof++;

					if (computed_prof % 100 == 0 || 
						(computed_prof % 10 == 0 && (double) computed_prof / (2 * gapped_sequences.size() - 1) > 0.95))
					{
						LOG_NORMAL << "Computing alignment - " << fixed << setprecision(1) << 100.0 * computed_prof / (2 * gapped_sequences.size()-1) << 
							"%    (" << computed_prof << " of " << (2 * gapped_sequences.size() - 1) << ")\r";
						fflush(stdout);
					}
				}
			}
		});

	for(auto p : workers)
	{
		p->join();
		delete p;
	}
	workers.clear();

	profile = profiles.begin()->second;

	return profile;
}

// *******************************************************************
bool CFAMSA::GetAlignment(vector<CGappedSequence*> &result)
{
	if (!final_profile) {
		return false;
	}		

	result = final_profile->data;

	return !result.empty();
}

#ifdef DEVELOPER_MODE
// *******************************************************************
bool CFAMSA::LoadRefSequences()
{
	// ***** Read input file	
	LOG_DEBUG << "Processing reference sequence set: " << params.ref_file_name << "\n";

	if (IOService::loadFasta(params.ref_file_name, ref_sequences) < 2)
	{
		LOG_NORMAL << "Error: no (or incorrect) input file\n";
		return false;
	}

	return true;
}
#endif

// *******************************************************************
bool CFAMSA::ComputeMSA(vector<CSequence>& sequences)
{
	adjustParams((int)sequences.size());
	
	LOG_VERBOSE 
		<< "Params:\n"
		<< "  no. of threads: " << params.n_threads << "\n"
		<< "  guide tree method: " << GT::toString(params.gt_method) << "\n"
		<< "  guide tree speeding up heuristic: " << GT::toString(params.gt_heuristic) << "\n";

	if (params.gt_method == GT::imported) {
		LOG_VERBOSE << "  guide tree file: " << params.guide_tree_in_file << "\n";
	}
		
	LOG_VERBOSE
		<< "  distance measure: " << dist2str(params.distance) << "\n\n"

		<< "Advanced params:\n"
		<< "  no. of refinements: " << params.n_refinements << "\n"
		<< "  refinement threshold: " << params.thr_refinement << "\n"
		<< "  scroing matrix: " << ScoringMatrices::toString(params.matrix_type) << "\n"
		<< "  gap open cost (rescalled): " << params.gap_open << "\n"
		<< "  gap extension cost (rescalled): " << params.gap_ext << "\n"
		<< "  gap terminal open cost (rescalled): " << params.gap_term_open << "\n"
		<< "  gap terminal extension cost (rescalled): " << params.gap_term_ext << "\n"
		<< "  gap cost scaller log-term: " << params.scaler_log << "\n"
		<< "  gap cost scaller div-term: " << params.scaler_div << "\n"
		<< "  enable gap rescaling: " << params.enable_gap_rescaling << "\n"
		<< "  enable gap optimization: " << params.enable_gap_optimization << "\n"
		<< "  enable total score calculation: " << params.enable_total_score_calculation << "\n"
		<< "  refinement mode: " << Refinement::toString(params.refinement_mode) << "\n"
		<< "  guided alignment radius: " << params.guided_alignment_radius << "\n\n";
	

#ifdef DEVELOPER_MODE
	if (params.test_ref_sequences)
		LoadRefSequences();
#endif

	string instr_names[] = { "None", "SSE", "SSE2", "SSE3", "SSE3S", "SSE41", "SSE42", "AVX", "AVX2", "AVX512", "NEON"};
	LOG_VERBOSE << "Hardware configuration: " << endl
		<< " Number of threads: " << params.n_threads << endl
		<< " Instruction set: " << instr_names[(int)instruction_set] << endl << endl;

	GuideTree tree;
	std::vector<CSequence*> mapped_seqs(sequences.size());

	// store distance matrix
	if (params.export_distances) {
		LOG_VERBOSE << "Calculating distances and storing in: " << params.output_file_name;
		std::shared_ptr<AbstractTreeGenerator> calculator = createTreeGenerator(params);
		std::transform(sequences.begin(), sequences.end(), mapped_seqs.begin(), [](CSequence& s) { return &s; });
		extendSequences(sequences);
		(*calculator)(mapped_seqs, tree.raw());
		shrinkSequences(sequences);
		LOG_VERBOSE << " [OK]" << endl;
		
	} 
	else {

		timers[TIMER_SORTING].StartTimer();
		LOG_VERBOSE << "Sorting sequences...";
		sortAndExtendSequences(sequences);
		std::transform(sequences.begin(), sequences.end(), mapped_seqs.begin(), [](CSequence& s) { return &s; });
		timers[TIMER_SORTING].StopTimer();
		LOG_VERBOSE << " [OK]" << endl;

		std::vector<int> original2mapped(sequences.size());
		std::iota(original2mapped.begin(), original2mapped.end(), 0);

		// remove duplicates
		int dups = 0;
		if (!params.keepDuplicates) {
			LOG_VERBOSE << "Duplicate removal... ";
			removeDuplicates(mapped_seqs, original2mapped);
			dups = sequences.size() - mapped_seqs.size();
			LOG_VERBOSE << mapped_seqs.size() << "/" << sequences.size() << " sequences retained." << endl;
		}

		// only one unique sequence - move input to output end exit
		if (mapped_seqs.size() == 1) {
			final_profile = new CProfile(&params);
			for (int i = 0; i < (int)sequences.size(); ++i) {
				final_profile->AppendRawSequence(CGappedSequence(move(sequences[i])));
			}
			return true;
		}

		// store mappings and temporarily reset numerical identifiers (to make medoid trees work)
		for (int i = 0; i < (int)mapped_seqs.size(); ++i) {
			mapped_seqs[i]->sequence_no = i;
		}
		
		timers[TIMER_TREE_BUILD].StartTimer();
		if (params.gt_method == GT::imported) {
			LOG_VERBOSE << "Importing guide tree from: " << params.guide_tree_in_file;
			tree.loadNewick(params.guide_tree_in_file, sequences);
			tree.toUnique(original2mapped, (int)mapped_seqs.size());
		}
		else {
			std::shared_ptr<AbstractTreeGenerator> gen = createTreeGenerator(params);
			LOG_VERBOSE << "Computing guide tree...";
			(*gen)(mapped_seqs, tree.raw());
		}
		shrinkSequences(sequences);
		LOG_VERBOSE << " [OK]" << endl;
		timers[TIMER_TREE_BUILD].StopTimer();

		if (params.export_tree) {
			// store guide tree in Newick format...
			timers[TIMER_TREE_STORE].StartTimer();
			LOG_VERBOSE << "Storing guide tree in: " << params.output_file_name;
			tree.fromUnique(original2mapped);
			tree.saveNewick(params.output_file_name, sequences);
			LOG_VERBOSE << " [OK]" << endl;
			timers[TIMER_TREE_STORE].StopTimer();
		}
		else {
			// ... or perform an alignment
			std::vector<CGappedSequence> gapped_sequences;
			std::vector<CGappedSequence*> mapped_gapped_seqs(mapped_seqs.size(), nullptr);
			
			// convert sequences into gapped sequences
			gapped_sequences.reserve(sequences.size());
			for (int i = 0; i < (int)sequences.size(); ++i) {
				gapped_sequences.emplace_back(std::move(sequences[i]));
				if (mapped_gapped_seqs[original2mapped[i]] == nullptr) {
					mapped_gapped_seqs[original2mapped[i]] = &gapped_sequences[i];
				}
			}

			// clear input vectors
			std::vector<CSequence>().swap(sequences); 
			std::vector<CSequence*>().swap(mapped_seqs);

			timers[TIMER_ALIGNMENT].StartTimer();
			LOG_VERBOSE << "Computing alignment...";
			final_profile = ComputeAlignment(mapped_gapped_seqs, tree.raw());
			if (!final_profile) {
				return false;
			}
			LOG_VERBOSE << "[OK]" << endl;
			timers[TIMER_ALIGNMENT].StopTimer();

			timers[TIMER_REFINMENT].StartTimer();
			LOG_VERBOSE << "Computing refinement...";
			if (!RefineAlignment(final_profile, params.n_threads))
				return false;
			LOG_VERBOSE << "[OK]" << endl;
			timers[TIMER_REFINMENT].StopTimer();

			if (final_profile->Size() != mapped_gapped_seqs.size()) {
				throw std::runtime_error("Error: incomplete guide tree - report a bug");
			}

			// compose ordered alignment of unique sequences (without duplicates)
			std::vector<CGappedSequence*> ordered_unique_alignment(final_profile->data.size(), nullptr);
			for (int is = 0; is < (int)final_profile->data.size(); ++is) {
				ordered_unique_alignment[final_profile->data[is]->sequence_no] = final_profile->data[is];
			}

			// compose final ordered alignment 
			std::vector<CGappedSequence*> ordered_alignment(gapped_sequences.size(), nullptr);
			for (int i = 0; i < (int)gapped_sequences.size(); ++i) {
				CGappedSequence& current = gapped_sequences[i];
				CGappedSequence* representative = ordered_unique_alignment[original2mapped[i]];

				if (current.original_no == representative->original_no) {
					// unique sequence - put representative in profile
					ordered_alignment[current.original_no] = representative;
				}
				else {
					// duplicate - make a copy of a representative
					CGappedSequence* duplicate = new CGappedSequence(*representative);
					duplicate->id = std::move(current.id);
					duplicate->original_no = current.original_no;
					ordered_alignment[current.original_no] = duplicate;
				}
			}

			final_profile->data = std::move(ordered_alignment);
		}

		// store stats
		if (params.areStatsStored()) {

			// this is time consuming - only in very verbose
			if (params.very_verbose_mode) {
				int64_t sackin = tree.calculateSackinIndex();
				statistics.put("guide_tree.sackin", sackin);
				statistics.put("guide_tree.sackin_norm", sackin / (double)sequences.size());
			}
			statistics.put("input.n_duplicates", dups);
			statistics.put("time.sort", timers[TIMER_SORTING].GetElapsedTime());
			statistics.put("time.tree_build", timers[TIMER_TREE_BUILD].GetElapsedTime());
			statistics.put("time.tree_store", timers[TIMER_TREE_STORE].GetElapsedTime());
			statistics.put("time.alignment", timers[TIMER_ALIGNMENT].GetElapsedTime());
			statistics.put("time.refinement", timers[TIMER_REFINMENT].GetElapsedTime());
		}
	}

	return true;
}

// *******************************************************************
bool CFAMSA::alignProfiles(vector<CGappedSequence>& p1, vector<CGappedSequence>& p2) {

	CProfile prof_1 = CProfile(&params);
	CProfile prof_2 = CProfile(&params);

	timers[TIMER_ALIGNMENT].StartTimer();
	LOG_VERBOSE << "Computing alignment...";

	for (const auto& gs: p1) {
		prof_1.AppendRawSequence(gs);
	}
	for (const auto& gs : p2) {
		prof_2.AppendRawSequence(gs);
	}

	prof_1.CalculateCountersScores();
	prof_2.CalculateCountersScores();
	
	uint32_t no_rows_per_box = 4;

	refresh::active_thread_pool_v2 atp(1, 1);

	final_profile = new CProfile(&prof_1, &prof_2, &params, 1, no_rows_per_box, &atp);
	if (!final_profile) {
		return false;
	}
	LOG_VERBOSE << "[OK]" << endl;
	timers[TIMER_ALIGNMENT].StopTimer();


	timers[TIMER_REFINMENT].StartTimer();
	LOG_VERBOSE << "Computing refinement...";
	if (!RefineAlignment(final_profile, 1))
		return false;
	LOG_VERBOSE << "[OK]" << endl;
	timers[TIMER_REFINMENT].StopTimer();

	if (params.verbose_mode || params.very_verbose_mode) {
		statistics.put("time.alignment", timers[TIMER_ALIGNMENT].GetElapsedTime());
		statistics.put("time.refinement", timers[TIMER_REFINMENT].GetElapsedTime());
	}

	return 0;
}
