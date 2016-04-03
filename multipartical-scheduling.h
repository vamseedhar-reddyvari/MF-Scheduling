//
// Created by Vamseedhar on 3/27/16.
//

#ifndef MFE_SCHEDULING_MULTIPARTICAL_SCHEDULING_H_H
#define MFE_SCHEDULING_MULTIPARTICAL_SCHEDULING_H_H

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>

using namespace std;

//Type defs
typedef vector<uint> uarray;
typedef vector<double> darray;
typedef unsigned int uint;
typedef vector<int> iarray;

// Global Constants
//unsigned int No_Users = 10000;
double PRECISION_STATE      = 0.01;              /* sampling period for state space */
double PRECISION_BID  = 0.15;            /* sampling period for bid space */
int MAXARRIVAL        =  100;                /* Arrival is a random procees, Maximum arrival in a time slot */
int SERVICE           =  500;                /* serivice in each time slot  */
double BETA           =  0.9;             /* 1 - regenration probability */
unsigned int MAXREGEN =  100;         /* regeneration yields a random queue length, Maximum value of regenerated queue */
unsigned int MINREGEN =    0;
unsigned int M        =   10;              /* Number of users in each cell */
unsigned int MAXBID   = 3500;          /*  Bid space trunction, max bid < 2000 */
unsigned int QFINAL   = 2000;          /*  State trunction, max queuelength < 15 :  MAXARRIVAL * (1/(1-beta)) * a large constant, Queue exceeds this value with small probability  */
unsigned int MAXQUEUE = 2500;
double epsilon        = 0.008;

darray generate_linspace(double start,double end, int No_Values);
iarray generate_intspace(double start,double end, int No_Values);


uint Total_Users= 100000;
uint No_Users_Per_Cell = 10;
uint No_Cells = Total_Users/No_Users_Per_Cell;
iarray Cells_Space = generate_intspace(0,No_Cells,No_Cells);


#endif //MFE_SCHEDULING_MULTIPARTICAL_SCHEDULING_H_H
