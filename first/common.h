//
// Created by kostas on 8/27/18.
//

#ifndef NTDS_COMMON_H
#define NTDS_COMMON_H
#include <algorithm>
#define BILLION  1000000000LL
timespec startT, endT;

void startTimer(){
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startT);
}

double endTimer(){
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endT);
    return endT.tv_sec * BILLION + endT.tv_nsec - (startT.tv_sec * BILLION + startT.tv_nsec);
}
#endif //NTDS_COMMON_H
