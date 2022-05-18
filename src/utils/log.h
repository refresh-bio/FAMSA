/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <cstdint>
#include <cstring>

#define LOG_NORMAL Log::getInstance(Log::LEVEL_NORMAL)
#define LOG_VERBOSE Log::getInstance(Log::LEVEL_VERBOSE)
#define LOG_DEBUG Log::getInstance(Log::LEVEL_DEBUG)

// *****************************************************************************************
//
class Log
{
public:
	static const int LEVEL_NORMAL;
	static const int LEVEL_VERBOSE;
	static const int LEVEL_DEBUG;

	void enable()	{ enabled = true; }
	void disable()	{ enabled = false; }
	bool isEnabled() const { return enabled; }
	
	// *****************************************************************************************
	//
	static Log& getInstance(int level) {
		static std::vector<std::shared_ptr<Log>> logs;
		if (logs.size() == 0) {
			logs.push_back(std::shared_ptr<Log>(new Log()));
			logs.push_back(std::shared_ptr<Log>(new Log()));
			logs.push_back(std::shared_ptr<Log>(new Log()));
		}

		return *logs[level];
	}

	// *****************************************************************************************
	//
	template <class T>
	Log& operator<<(T v) {
		if (enabled) { *out << v; }
		return *this;
	}

	Log& operator<< (std::ostream& (*pf)(std::ostream&));
	Log& operator<< (std::ios& (*pf)(std::ios&));
	Log& operator<< (std::ios_base& (*pf)(std::ios_base&));

	static std::string formatLargeNumber(uint64_t num, int minWidth = 0);

protected:
	bool enabled;
	std::ostream* out;

	Log();
};



// ************************************************************************************
class NumericConversions
{
public:
	static char digits[100000 * 5];
	static uint64_t powers10[15];
	struct _si {
		_si()
		{
			for (int i = 0; i < 100000; ++i)
			{
				int dig = i;

				digits[i * 5 + 4] = '0' + (dig % 10);
				dig /= 10;
				digits[i * 5 + 3] = '0' + (dig % 10);
				dig /= 10;
				digits[i * 5 + 2] = '0' + (dig % 10);
				dig /= 10;
				digits[i * 5 + 1] = '0' + (dig % 10);
				dig /= 10;
				digits[i * 5 + 0] = '0' + dig;
			}

			powers10[0] = 1;
			for (int i = 1; i < 15; ++i)
				powers10[i] = 10 * powers10[i - 1];
		}
	} static _init;

	static int NDigits(uint64_t v)
	{
		return (v < 10000)
			? (v < 100 ? (v < 10 ? 1 : 2) : (v < 1000 ? 3 : 4))
			: (v < 1000000 ? (v < 100000 ? 5 : 6) : (v < 10000000 ? 7 : 8));
	}

	static int Int2PChar(uint64_t val, char *str)
	{
		if (val >= 1000000000000000ull)
		{
			uint64_t dig1 = val / 1000000000000000ull;
			val -= dig1 * 1000000000000000ull;
			uint64_t dig2 = val / 10000000000ull;
			val -= dig2 * 10000000000ull;
			uint64_t dig3 = val / 100000ull;
			uint64_t dig4 = val - dig3 * 100000ull;

			int ndig = NDigits(dig1);

			std::memcpy(str, digits + dig1 * 5 + (5 - ndig), ndig);
			std::memcpy(str + ndig, digits + dig2 * 5, 5);
			std::memcpy(str + ndig + 5, digits + dig3 * 5, 5);
			std::memcpy(str + ndig + 10, digits + dig4 * 5, 5);

			return ndig + 15;
		}
		else if (val >= 10000000000ull)
		{
			uint64_t dig1 = val / 10000000000ull;
			val -= dig1 * 10000000000ull;
			uint64_t dig2 = val / 100000ull;
			uint64_t dig3 = val - dig2 * 100000ull;

			int ndig = NDigits(dig1);

			std::memcpy(str, digits + dig1 * 5 + (5 - ndig), ndig);
			std::memcpy(str + ndig, digits + dig2 * 5, 5);
			std::memcpy(str + ndig + 5, digits + dig3 * 5, 5);

			return ndig + 10;
		}
		else if (val >= 100000ull)
		{
			uint64_t dig1 = val / 100000ull;
			uint64_t dig2 = val - dig1 * 100000ull;

			int ndig = NDigits(dig1);

			memcpy(str, digits + dig1 * 5 + (5 - ndig), ndig);
			memcpy(str + ndig, digits + dig2 * 5, 5);

			return ndig + 5;
		}
		else
		{
			int ndig = NDigits(val);

			memcpy(str, digits + val * 5 + (5 - ndig), ndig);

			return ndig;
		}
	}

	static int Double2PChar(double val, uint32_t prec, char *str)
	{
		int64_t a = (int64_t)val;
		int64_t b = (int64_t)((1.0 + (val - (double)a)) * powers10[prec] + 0.5);

		int r1 = Int2PChar(a, str);
		int r2 = Int2PChar(b, str + r1);
		str[r1] = '.';

		return r1 + r2;
	}
};
