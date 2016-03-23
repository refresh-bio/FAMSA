#pragma once
#include <iostream>
#include <vector>
#include <memory>

#define LOG_NORMAL Log::getInstance(Log::LEVEL_NORMAL)
#define LOG_DEBUG Log::getInstance(Log::LEVEL_DEBUG)


class Log
{
public:
	static const int LEVEL_NORMAL;
	static const int LEVEL_DEBUG;

	void enable()	{ enabled = true; }
	void disable()	{ enabled = false; }
	
	static Log& getInstance(int level) {
		static std::vector<std::shared_ptr<Log>> logs;
		if (logs.size() == 0) {
			logs.push_back(std::shared_ptr<Log>(new Log()));
			logs.push_back(std::shared_ptr<Log>(new Log()));
		}

		return *logs[level];
	}

	template <class T>
	Log& operator<<(T v) {
		if (enabled) { *out << v; }
		return *this;
	}

	Log& operator<< (std::ostream& (*pf)(std::ostream&));
	Log& operator<< (std::ios& (*pf)(std::ios&));
	Log& operator<< (std::ios_base& (*pf)(std::ios_base&));

protected:
	bool enabled;
	std::ostream* out;


	Log();
};



