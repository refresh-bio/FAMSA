#pragma once
#include <string>
#include <map>
#include <sstream>

#define GET_STATS

#ifdef GET_STATS
	#define STATS_WRITE(key,value) statistics[key] = value
	#define STATS_ADD(key,value) statistics[key] += value
#else
	#define STATS_WRITE(key,value) 
	#define STATS_ADD(key,value) 
#endif

class StatisticsProvider
{
public:
	std::map<std::string, double> statistics;

	virtual ~StatisticsProvider() {}

	void joinStats(const StatisticsProvider& other);

	void clearStats() { statistics.clear(); }

	std::string printStats();

	void saveStats(std::string file);
};