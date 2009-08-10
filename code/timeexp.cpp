// 
// File: timeexp.cpp -- function definitions for StartTimer and 
// ElapsedTime, functions for timing experiments
// Source: Shiflet, Angela B. Data Structures in C++a
//
#include "timeexp.h"

static time_t StartingTime;   // starting time of clock 

// 
// Function to initialize static global variable StartingTime
// to the number of ticks of the clock
// Pre:  none
// Post: StartingTime has been initialized to clock()
//

void StartTimer(void)
{
   StartingTime = clock();
}

//
// Function to return the elapsed time in seconds since 
// static global variable StartingTime was initialized
// Pre:  StartingTime has been initialized to a value of clock.
// Post: The elapsed time since the initialization of StartingTime 
//       has been returned.
//

double ElapsedTime(void)
{
   //return (double(clock() - StartingTime));
   return (double(clock() - StartingTime)/CLOCKS_PER_SEC);
}
