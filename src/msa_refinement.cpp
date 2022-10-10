#include "msa.h"

#include "./utils/log.h"
#include <iomanip>

// *******************************************************************
void CFAMSA::RefineRandom(CProfile* profile_to_refine, vector<size_t>& dest_prof_id)
{
	for (size_t i = 0; i < profile_to_refine->data.size(); ++i)
		dest_prof_id.emplace_back(rnd_rfn() % 2);

	if (count(dest_prof_id.begin(), dest_prof_id.end(), 0) == 0 ||
		count(dest_prof_id.begin(), dest_prof_id.end(), 1) == 0)		// Both profiles must contain at least 1 sequence
	{
		size_t id = rnd_rfn() % dest_prof_id.size();
		dest_prof_id[id] = !dest_prof_id[id];
	}
}

// *******************************************************************
void CFAMSA::RefineMostEmptyAndFullColumn(CProfile* profile_to_refine, vector<size_t>& dest_prof_id, vector<size_t>& gap_stats, bool valid_gap_stats)
{
	size_t size = profile_to_refine->data.front()->gapped_size;
	size_t card = profile_to_refine->data.size();

	dest_prof_id.clear();

	if (!valid_gap_stats)
		profile_to_refine->GetGapStats(gap_stats);

	vector<pair<size_t, size_t>> tmp;

	for (size_t i = 1; i <= size; ++i)
	{
		int x = (int)min(gap_stats[i], card - gap_stats[i]);
		if (x > 0)
			tmp.emplace_back(i, x);
	}

	stable_sort(tmp.begin(), tmp.end(), [](const pair<size_t, size_t>& x, const pair<size_t, size_t>& y) {
		if (x.second != y.second)
			return x.second < y.second;
		else
			return x.first < y.first;
		});

	if (tmp.empty())
	{
		RefineRandom(profile_to_refine, dest_prof_id);
		return;
	}

	size_t col_id = tmp[rnd_rfn() % tmp.size()].first;

	int first_prof_id = 0;
	int second_prof_id = 1;

	if (profile_to_refine->data[0]->GetSymbol(col_id) == GAP)
		swap(first_prof_id, second_prof_id);

	for (size_t j = 0; j < card; ++j)
		if (profile_to_refine->data[j]->GetSymbol(col_id) == GAP)
			dest_prof_id.emplace_back(first_prof_id);
		else
			dest_prof_id.emplace_back(second_prof_id);
}

// *******************************************************************
// Refine alignment
#ifdef DEBUG_MODE
bool CFAMSA::RefineAlignment(string output_file_name)
#else
bool CFAMSA::RefineAlignment(CProfile*& profile_to_refine)
#endif
{
	// Restart generator
	rnd_rfn.seed(5489u);

	if (params.refinement_mode == Refinement::OFF ||
		(params.refinement_mode == Refinement::AUTO && profile_to_refine->Size() > params.thr_refinement)) {
		return true;
	}

	size_t n_ref = params.n_refinements;
	size_t n_seq = profile_to_refine->Size();

	vector<size_t> gap_stats;

	if (n_ref > 2 * n_seq)
		n_ref = 2 * n_seq;
	if (n_ref > 0 && n_ref < 100 && n_seq < 100)
		n_ref = 100;

#ifdef DEBUG_MODE
	FILE* f_stat;

	if (output_file_name != "")
	{
		vector<CGappedSequence*> result;
		GetAlignment(result);

		COutputFile out_file;

		out_file.PutSequences(result);
		out_file.SaveFile(output_file_name + to_string(0));

		f_stat = fopen((output_file_name + "_stats").c_str(), "wt");

		fprintf(f_stat, "%d  %f\n", final_profile->width, final_profile->CalculateTotalScore());
	}
#endif

	int n_ref_succ = 0;
	score_t prev_total_score = profile_to_refine->CalculateTotalScore();

	sort(profile_to_refine->data.begin(), profile_to_refine->data.end(), [](CGappedSequence* p, CGappedSequence* q) {return p->id < q->id; });

	vector<size_t> dest_prof_id;
	vector<vector<size_t>> old_dest_prof_ids;

	vector<int> column_mapping1, column_mapping2;

	size_t i_ref;
	size_t i_succ_ref;
	bool valid_gap_stats = false;
#ifdef DEBUG_MODE
	int ref_upd[2] = { 0 };
	int hist_size[20] = { 0 };
#endif

	for (i_ref = i_succ_ref = 0; i_succ_ref < n_ref && i_ref < 20 * n_ref; ++i_ref)
	{
		LOG_DEBUG << "Computing refinement - " << fixed << setprecision(1) << 100.0 * (double)i_succ_ref / (double)n_ref << "%    (" << i_succ_ref << " of " << n_ref << ")  \r";

		CProfile profile1(&params), profile2(&params);

		RefineMostEmptyAndFullColumn(profile_to_refine, dest_prof_id, gap_stats, valid_gap_stats);
		valid_gap_stats = true;

		if (find(old_dest_prof_ids.begin(), old_dest_prof_ids.end(), dest_prof_id) == old_dest_prof_ids.end())
		{
			// Split into two profiles
			for (size_t i = 0; i < profile_to_refine->data.size(); ++i)
				if (dest_prof_id[i])
					profile1.AppendRawSequence(*profile_to_refine->data[i]);
				else
					profile2.AppendRawSequence(*profile_to_refine->data[i]);

			// Condense the profiles (remove empty columns)
			profile1.Condense(column_mapping1);
			profile2.Condense(column_mapping2);

			profile1.OptimizeGaps();
			profile2.OptimizeGaps();

			profile1.Size();
			profile2.Size();

#ifdef DEBUG_MODE
			int size_min = min(p1_size, p2_size);
			hist_size[min(9, size_min)]++;
#endif

			CProfile* prof = new CProfile(&params);

			// TODO: Enable parallelization here!
			prof->Align(&profile1, &profile2, 1, 0, &column_mapping1, &column_mapping2);
			sort(prof->data.begin(), prof->data.end(), [](CGappedSequence* p, CGappedSequence* q) {return p->id < q->id; });

			if (!(*prof == *profile_to_refine))		// if the new profile is the same as previous do not score it
			{
				prof->CalculateTotalScore();
#ifdef DEBUG_MODE
				ref_upd[0]++;
#endif

				if (prof->total_score >= prev_total_score)
				{
					prev_total_score = prof->total_score;
					swap(profile_to_refine, prof);
					++n_ref_succ;
					old_dest_prof_ids.clear();
					valid_gap_stats = false;
#ifdef DEBUG_MODE
					ref_upd[1]++;
					hist_size[10 + min(9, size_min)]++;
#endif
				}
			}

			delete prof;

			old_dest_prof_ids.emplace_back(dest_prof_id);
			i_succ_ref++;

#ifdef DEBUG_MODE
			if (output_file_name != "")
			{
				vector<CGappedSequence*> result;
				GetAlignment(result);

				COutputFile out_file;

				out_file.PutSequences(result);
				out_file.SaveFile(output_file_name + to_string(i_ref + 1));
				fprintf(f_stat, "%d  %f  p.sizes: %5d %5d\n", final_profile->width, final_profile->total_score, p1_size, p2_size);
			}
#endif
		}
	}

#ifdef DEBUG_MODE
	if (output_file_name != "")
		fclose(f_stat);
#endif

	return true;
}