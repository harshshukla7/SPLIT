#include <mach/mach_time.h>

/* Use MAC OSX  mach_time for timing */
typedef struct splitTimer{
  uint64_t tic;
  uint64_t toc;
  mach_timebase_info_data_t tinfo;
  
} splitTimer;

void   split_tic();
double split_toc();
