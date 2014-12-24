#include "fast_math.h"
#include <stdio.h>
#include <sys/time.h>
#include <math.h>

// call this function to start a microsecond-resolution timer
struct timeval timer_start(){
  struct timeval start_time;
  struct timezone tz;
  gettimeofday(&start_time, &tz);
  return start_time;
}

// call this function to end a timer, returning microseconds elapsed as a long
long timer_end(struct timeval start_time){
  struct timeval end_time;
  struct timezone tz;
  gettimeofday(&end_time, &tz);
  long elapsed = (end_time.tv_sec - start_time.tv_sec)*1e6 + 
    end_time.tv_usec - start_time.tv_usec;
  return elapsed;
}

int main (int argc, char *argv[]) {
  double x,y,v;
  struct timeval tstart;
  double elapsed;
  
  tstart = timer_start();
  for (y=-5e2; y<5e2; y++)
    for (x=-5e2; x<5e2; x++)
      v = poly_atan2(y, x);
  elapsed = timer_end(tstart)/1.e3;
  printf("x=%f, y=%f, v=%f", x, y, v);
  printf("Elapsed time: %.3g ms\n",elapsed);

  tstart = timer_start();
  for (y=-5e2; y<5e2; y++)
    for (x=-5e2; x<5e2; x++)
      v = atan2(y, x);
  elapsed = timer_end(tstart)/1.e3;
  printf("x=%f, y=%f, v=%f", x, y, v);
  printf("Elapsed time: %.3g ms\n",elapsed);

  return 0;
}
