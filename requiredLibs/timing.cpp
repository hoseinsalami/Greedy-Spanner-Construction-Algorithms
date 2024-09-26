/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2012  Bauke Conijn <email>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cassert>
#include <cstdio>
#include "timing.h"

#define USE_QFPC // escape to use old counter

#if defined _WIN32 || defined _WIN64
#ifdef USE_QFPC
// High resolution windows timer
#include <windows.h>
#include <winbase.h>
static LARGE_INTEGER begin, end, freq;

void Timing::reset()
{
	bool r = QueryPerformanceFrequency(&freq); // obtain frequency in seconds.
	//assert(r);
	QueryPerformanceCounter(&begin);
}

double Timing::elapsed()
{
	QueryPerformanceCounter(&end);
	return (end.QuadPart - begin.QuadPart)*1000./(double)freq.QuadPart;
}
#else
// Low resolution windows timer
#include <timer.h>
static unsigned int time = 0;

void Timing::reset()
{
	startTime();
	time = calculateElapsedTime();
}

double Timing::elapsed()
{
	return calculateElapsedTime() - time;
}
#endif
#else
#ifdef USE_QFPC
// High resolution linux timer
#include <time.h>
static timespec begin, end, freq;

void Timing::reset()
{
	int r = clock_getres(CLOCK_MONOTONIC, &freq);
	//assert(r==0); NOTE: commented this line because it gave compile errors
	//printf("Timer rsolution is %ldns\n",freq.tv_nsec);
	clock_gettime(CLOCK_MONOTONIC, &begin);
}

double Timing::elapsed()
{
	clock_gettime(CLOCK_MONOTONIC, &end);
	return (end.tv_sec-begin.tv_sec)*1000.0 + (end.tv_nsec-begin.tv_nsec)/1000000.0;
}
#else
// Low resolution linux timer
#include <time.h>
static time_t begin, end;

void Timing::reset()
{
	time(&begin);
}

double Timing::elapsed()
{
	time(&end);
	return (end - begin)* 1000.0 / CLOCKS_PER_SEC;
}
#endif
#endif

