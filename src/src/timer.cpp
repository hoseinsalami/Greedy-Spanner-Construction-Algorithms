#include "timer.h"

#include <iostream>

using namespace std;

LARGE_INTEGER timerFreq_;
LARGE_INTEGER counterAtStart_;

void startTime()
{
  QueryPerformanceFrequency(&timerFreq_);
  QueryPerformanceCounter(&counterAtStart_);
  TIMECAPS ptc;
  UINT cbtc = 8;
  MMRESULT result = timeGetDevCaps(&ptc, cbtc);
}

unsigned int calculateElapsedTime()
{
  if (timerFreq_.QuadPart == 0)
  {
    return -1;
  }
  else
  {
    LARGE_INTEGER c;
    QueryPerformanceCounter(&c);
    return static_cast<unsigned int>( (c.QuadPart - counterAtStart_.QuadPart) * 1000 / timerFreq_.QuadPart );
  }
}
