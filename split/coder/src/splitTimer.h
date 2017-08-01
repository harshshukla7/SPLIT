#ifndef __SPLIT_TIMER
#define __SPLIT_TIMER

#include <mach/mach_time.h>

/* Use MAC OSX  mach_time for timing */
typedef struct splitTimer{
  uint64_t tic;
  uint64_t toc;
  mach_timebase_info_data_t tinfo;
  
} splitTimer;

double  split_tic();
double split_toc(double start);

#endif