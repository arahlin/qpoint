#include "fast_math.h"
#include <stdio.h>
#include <sys/time.h>
#include <math.h>

// call this function to start a microsecond-resolution timer
struct timeval timer_start(){
  struct timeval start_time;
  gettimeofday(&start_time, NULL);
  return start_time;
}

// call this function to end a timer, returning microseconds elapsed as a long
long timer_end(struct timeval start_time){
  struct timeval end_time;
  gettimeofday(&end_time, NULL);
  long elapsed = (end_time.tv_sec - start_time.tv_sec)*1e6 + 
    end_time.tv_usec - start_time.tv_usec;
  return elapsed;
}

int main (int argc, char *argv[]) {
  double x,y,v;
  struct timeval tstart;
  double elapsed;
  int len = 1e3;
  
  tstart = timer_start();
  for (y=-len; y<len; y++)
    for (x=-len; x<len; x++)
      v = poly_atan2(y, x);
  elapsed = timer_end(tstart)/1.e3;
  printf("Polynomial x=%f, y=%f, v=%f\n", x, y, v);
  printf("Polynomial elapsed time: %.3g ms for %d operations\n",elapsed,4*len*len);

  tstart = timer_start();
  for (y=-len; y<len; y++)
    for (x=-len; x<len; x++)
      v = atan2(y, x);
  elapsed = timer_end(tstart)/1.e3;
  printf("Built-in x=%f, y=%f, v=%f\n", x, y, v);
  printf("Built-in elapsed time: %.3g ms for %d operations\n",elapsed,4*len*len);

  return 0;
}
