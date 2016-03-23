#pragma once
#include <algorithm>
#include <numeric>
#include <vector>

#ifdef _WIN32
	#include <windows.h>
#else
	#include <sys/time.h>
	#include <sys/resource.h>
	#include <unistd.h>
#endif


#define GET_TIMES

#ifdef GET_TIMES
	#define TIMER_CREATE(timer) Timer timer
	#define TIMER_START(timer) (timer).start()
	#define TIMER_MOVEON(timer) (timer).moveOn()
	#define TIMER_STOP(timer) (timer).stop()
	#define TIMER_SECONDS(timer) (timer).seconds()
	#define TIMER_MSECONDS(timer) (timer).seconds()
	#define TIMER_STOP_SAVE(timer, var) var=(timer).stop()
#else
	#define TIMER_CREATE(timer)
	#define TIMER_START(timer)
	#define TIMER_MOVEON(timer)
	#define TIMER_STOP(timer)
	#define TIMER_SECONDS(timer)
	#define TIMER_MSECONDS(timer)
	#define TIMER_STOP_SAVE(timer, var)
#endif


typedef struct 
{
	timeval start;
	timeval stop;
} stop_watch_t;


class Timer
{
public:
	/// <summary>
	/// Gets measured time in milliseconds.
	/// </summary>
	double mseconds() { return value; }

	/// <summary>
	/// Gets measured time in seconds.
	/// </summary>
	double seconds() { return value / 1000.0; }

	/// <summary>
	/// Initializes some stuff.
	/// </summary>
	Timer()
	{
	#ifdef _WIN32
		QueryPerformanceFrequency(&frequency);
	#endif
		reset();
	};

	/// <summary>
	/// Starts time measurement.
	/// </summary>
	void start()
	{
	#ifdef _WIN32
		QueryPerformanceCounter(&qp1);
	#else
		gettimeofday(&(timer.start),NULL);
	#endif
		value = 0.0;
	};

	/// <summary>
	/// Stops time measurement.
	/// </summary>
	/// <return>Time in seconds.</return>
	double stop()
	{
	#ifdef _WIN32
		QueryPerformanceCounter(&qp2);
		value += (1000000.0 * (double) (qp2.QuadPart - qp1.QuadPart)) / frequency.QuadPart / 1000.0;
	#else
		gettimeofday(&(timer.stop),NULL);
		timeval res;
		timersub(&(timer.stop),&(timer.start),&res);
		value += res.tv_sec * 1000.0 + res.tv_usec/1000.0; // 10^3 mSec per second

	#endif
		return value / 1000;
	};

	/// <summary>
	/// Continues time measurement.
	/// </summary>
	void moveOn()
	{
	#ifdef _WIN32
		QueryPerformanceCounter(&qp1);
	#else
		gettimeofday(&(timer.start),NULL);
	#endif
	}

	/// <summary>
	/// Resets the measured time.
	/// </summary>
	void reset()
	{
		value = 0.0;
	}

private:
	double value;

#ifdef _WIN32
	LARGE_INTEGER qp1;
	LARGE_INTEGER qp2;
	LARGE_INTEGER frequency;
#else
	stop_watch_t timer;
#endif

};
