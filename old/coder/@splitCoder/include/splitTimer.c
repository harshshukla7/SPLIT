#include <stdio.h>
#include <stdlib.h>
#include "splitTimer.h"

splitTimer solvertimer;

void split_tic()
{
  /* read current clock cycles */
  solvertimer.tic = mach_absolute_time();
}

double split_toc()
{
  uint64_t duration; /* elapsed time in clock cycles*/
  solvertimer.toc = mach_absolute_time();
  duration = solvertimer.toc - solvertimer.tic;
  
  /*conversion from clock cycles to nanoseconds*/
  mach_timebase_info(&(solvertimer.tinfo));
  duration *= solvertimer.tinfo.numer;
  duration /= solvertimer.tinfo.denom;
  
  return (double)duration / 1000000000;
}
