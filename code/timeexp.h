// 
// File: timeexp.h -- function prototypes for StartTimer and 
// ElapsedTime, functions for timing experiments
// Source: Shiflet, Angela B. Data Structure in C++
//

#ifndef TIMEEXP_H
#define TIMEEXP_H

#include <time.h>

void   StartTimer(void);
  // Usage: StartTimer();
  // NOTE:  Call this function just before the beginning of the code which
  // 	    you would like to measure the running time for.

double ElapsedTime(void);
  // Usage: cout << "The running time is: "; 
  //        cout << ElapsedTime() << " milliseconds."<< endl;
  // NOTE:  1) The ElapsedTime() function returns the elapsed time 
  //	       (in milliseconds) since the StartTimer() function was 
  //           called.
  //        2) A millisecond is equal to one thousandth of a second.
#endif
