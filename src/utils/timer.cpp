/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifdef WIN32
#include <windows.h>
#endif

#include <cstdio> // NULL
#include "timer.h"


#ifdef WIN32
// **********************************************************
double CStopWatch::LIToSecs( LARGE_INTEGER & L) 
{
	return ((double)L.QuadPart /(double)frequency.QuadPart);
}

// **********************************************************
// CStopWatch
// **********************************************************
CStopWatch::CStopWatch()
{
	timer.start.QuadPart=0;
	timer.stop.QuadPart=0;	

	QueryPerformanceFrequency( &frequency );
}

// **********************************************************
void CStopWatch::StartTimer( ) 
{
    QueryPerformanceCounter(&timer.start);
}

// **********************************************************
void CStopWatch::StopTimer( ) 
{
    QueryPerformanceCounter(&timer.stop);
}

// **********************************************************
double CStopWatch::GetElapsedTime() 
{
	LARGE_INTEGER time;
	time.QuadPart = timer.stop.QuadPart - timer.start.QuadPart;
    return LIToSecs(time);
}


#else
// **********************************************************
CStopWatch::CStopWatch()
{
	gettimeofday(&(timer.start), NULL);
	timer.stop = timer.start;
}

// **********************************************************
void CStopWatch::StartTimer( ) 
{
	gettimeofday(&(timer.start),NULL);
}

// **********************************************************
void CStopWatch::StopTimer( ) 
{
	gettimeofday(&(timer.stop),NULL);
}

// **********************************************************
double CStopWatch::GetElapsedTime() 
{	
	timeval res;
	timersub(&(timer.stop),&(timer.start),&res);
	return res.tv_sec + res.tv_usec/1000000.0; // 10^6 uSec per second
}


#endif