#include <assert.h>
#include <CoreServices/CoreServices.h>
#include <mach/mach.h>
#include <mach/mach_time.h>
#include <unistd.h>
#include "splitTimer.h"

static uint64_t start;

double split_tic()
{
  /* read current clock cycles */
  return (double)mach_absolute_time();
}

double split_toc(double start)
{
  double          elapsed;
  uint64_t        elapsedNano;
  static mach_timebase_info_data_t    sTimebaseInfo;

  uint64_t end = mach_absolute_time();

    // Calculate the duration.

  elapsed = (double)end - start;

    // Convert to nanoseconds.

    // If this is the first time we've run, get the timebase.
    // We can use denom == 0 to indicate that sTimebaseInfo is 
    // uninitialised because it makes no sense to have a zero 
    // denominator is a fraction.

  if ( sTimebaseInfo.denom == 0 ) {
    (void) mach_timebase_info(&sTimebaseInfo);
  }

    // Do the maths. We hope that the multiplication doesn't 
    // overflow; the price you pay for working in fixed point.

  return elapsed * (double)(sTimebaseInfo.numer) / (double)(sTimebaseInfo.denom);
}
