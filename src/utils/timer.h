/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _TIMER_H
#define _TIMER_H

#ifdef WIN32
#include <windows.h>

// **********************************************************
typedef struct 
{
	LARGE_INTEGER start;
	LARGE_INTEGER stop;
} stop_watch_t;

// **********************************************************
class CStopWatch 
{
	stop_watch_t timer;
	LARGE_INTEGER frequency;
	double LIToSecs( LARGE_INTEGER & L);

public:
	CStopWatch();
	void StartTimer( );
	void StopTimer( );
	double GetElapsedTime();
};

// **********************************************************
typedef struct
{
	ULARGE_INTEGER start;
	ULARGE_INTEGER stop;
} thread_watch_t;



#else
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

typedef struct 
{
	timeval start;
	timeval stop;
} stop_watch_t;

class CStopWatch 
{
	stop_watch_t timer;

public:
	CStopWatch();
	void StartTimer( );
	void StopTimer( );
	double GetElapsedTime();
};

typedef timeval thread_watch_t;


#endif

#endif
// ***** EOF
